/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Wikki GmbH
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Field-level helpers for foamCoherentConvert:
    reading typed fields, permuting patch values into the coherent
    face order and writing them with an explicit stream format.

\*---------------------------------------------------------------------------*/

#ifndef Foam_coherentConvertFields_H
#define Foam_coherentConvertFields_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOobjectList.H"
#include "HashSet.H"
#include "HashTable.H"
#include "directFvPatchFieldMapper.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Per-patch permutation from the legacy face order into the coherent
//- face order, keyed by patch name.
//  The coherent mesh writer stable-sorts all faces by their owner cell
//  and the reader keeps that order within each patch, so the coherent
//  order of a patch is Foam::sortedOrder(faceCells) of the legacy patch.
//  Only non-identity permutations are stored.
struct patchPermutations
{
    //- Permutation: coherent index -> legacy index
    HashTable<labelList> order;

    //- Legacy face centres in coherent order (for verification)
    HashTable<pointField> centres;

    //- Extract from a legacy mesh
    void reset(const polyMesh& mesh)
    {
        order.clear();
        centres.clear();

        for (const polyPatch& pp : mesh.boundaryMesh())
        {
            labelList perm(Foam::sortedOrder(pp.faceCells()));

            centres.set(pp.name(), pointField(pp.faceCentres(), perm));

            bool identity = true;
            forAll(perm, i)
            {
                if (perm[i] != i)
                {
                    identity = false;
                    break;
                }
            }

            if (!identity)
            {
                order.set(pp.name(), std::move(perm));
            }
        }
    }

    //- Verify the ordering assumption against the coherent mesh.
    //  Neighbour-side processor patches are not checked (their values
    //  are never written) and are removed from the permutation table.
    void verify(const polyMesh& mesh)
    {
        const scalar tol = 1e-8*mesh.bounds().mag();

        for (const polyPatch& pp : mesh.boundaryMesh())
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && !ppp->owner())
            {
                order.erase(pp.name());
                centres.erase(pp.name());
                continue;
            }

            const auto iter = centres.cfind(pp.name());
            if (!iter.good())
            {
                continue;
            }

            const pointField& legacy = iter.val();
            const vectorField coherent(pp.faceCentres());

            if (legacy.size() != coherent.size())
            {
                FatalErrorInFunction
                    << "Patch " << pp.name() << " has " << coherent.size()
                    << " faces on the coherent mesh but " << legacy.size()
                    << " on the legacy mesh" << nl
                    << exit(FatalError);
            }

            forAll(legacy, facei)
            {
                if (mag(legacy[facei] - coherent[facei]) > tol)
                {
                    FatalErrorInFunction
                        << "Face ordering of patch " << pp.name()
                        << " on the coherent mesh does not match the"
                        << " expected order (face " << facei << ": "
                        << coherent[facei] << " vs " << legacy[facei]
                        << ")" << nl
                        << exit(FatalError);
                }
            }
        }
    }
};


//- Permute the patch values of a GeometricField into coherent order
template<class Type, template<class> class PatchField, class GeoMesh>
void permutePatchFields
(
    GeometricField<Type, PatchField, GeoMesh>& fld,
    const patchPermutations& perms
)
{
    auto& bfld = fld.boundaryFieldRef();

    forAll(bfld, patchi)
    {
        const auto iter = perms.order.cfind(bfld[patchi].patch().name());

        if (iter.good() && bfld[patchi].size() == iter.val().size())
        {
            bfld[patchi].autoMap(directFvPatchFieldMapper(iter.val()));
        }
    }
}


//- No patch values for a DimensionedField
template<class Type, class GeoMesh>
void permutePatchFields
(
    DimensionedField<Type, GeoMesh>&,
    const patchPermutations&
)
{}


//- Construct a GeometricField without reading its old-time field,
//- which is converted like any other field
template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>* newField
(
    const IOobject& io,
    const fvMesh& mesh,
    GeometricField<Type, PatchField, GeoMesh>*
)
{
    return new GeometricField<Type, PatchField, GeoMesh>(io, mesh, false);
}


//- Construct a DimensionedField
template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>* newField
(
    const IOobject& io,
    const fvMesh& mesh,
    DimensionedField<Type, GeoMesh>*
)
{
    return new DimensionedField<Type, GeoMesh>(io, mesh);
}


//- Read the named fields of the given type, optionally permute the
//- patch values, and write with the given stream format.
//  Fields remain registered (in storedObjects) so that later fields
//  can look them up. Returns the number of fields converted.
template<class FieldType>
label convertFields
(
    const fvMesh& mesh,
    const wordList& fieldNames,
    const patchPermutations* perms,
    const IOstreamOption streamOpt,
    const bool writeNow,
    DynamicList<regIOobject*>& storedObjects,
    wordHashSet& converted
)
{
    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        if (!nFields)
        {
            Info<< "    " << FieldType::typeName << ':';
        }
        Info<< ' ' << fieldName;

        FieldType* fldPtr = newField
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh.thisDb(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh,
            static_cast<FieldType*>(nullptr)
        );

        if (!fldPtr->store())
        {
            FatalErrorInFunction
                << "Field " << fieldName << " is already registered"
                << exit(FatalError);
        }
        storedObjects.push_back(fldPtr);

        if (perms)
        {
            permutePatchFields(*fldPtr, *perms);
        }

        if (writeNow)
        {
            if (!fldPtr->writeObject(streamOpt, true))
            {
                FatalErrorInFunction
                    << "Failed writing field " << fieldName
                    << exit(FatalError);
            }
        }

        converted.insert(fieldName);
        ++nFields;
    }

    if (nFields) Info<< nl;

    return nFields;
}


//- Select the names with the given header class
inline wordList selectNames
(
    const word& className,
    const wordList& names,
    const wordList& classes,
    const boolList& mask
)
{
    DynamicList<word> selected(names.size());

    forAll(names, i)
    {
        if (mask[i] && classes[i] == className)
        {
            selected.push_back(names[i]);
        }
    }

    Foam::sort(selected);

    return wordList(std::move(selected));
}

} // End namespace Foam

#endif

// ************************************************************************* //
