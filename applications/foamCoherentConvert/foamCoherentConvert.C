/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Sergey Lesnik
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

Application
    foamCoherentConvert

Group
    grpMiscUtilities

Author
    Sergey Lesnik, Wikki GmbH

Description
    Converts the mesh and the fields of a case between the legacy formats
    (ascii/binary) and the coherent (ADIOS2) format, in place.

    Supported are the polyMesh (without zones) and the volume, surface and
    volume-internal fields, i.e. everything the coherent format can hold.
    Other objects (point fields, lagrangian clouds, sets, zones, ...) are
    reported and left untouched.

    The conversion does not depend on the 'writeFormat' entry of the
    controlDict, but the entry has to be adjusted afterwards to match the
    new format of the case.

    Legacy to coherent (target coherent):
    - The current partitioning is retained: a serial run yields a single
      partition, a parallel run on N ranks yields N partitions.
      Use 'renumberMesh -decompose' for repartitioning.
    - Since the coherent format orders the faces of each patch by owner
      cell, non-uniform patch values are permuted accordingly.
    - For serial runs the field files are replaced by their coherent
      counterparts (same location).

    Coherent to legacy (target ascii or binary):
    - A parallel run writes a decomposed case into the processor
      directories (without procAddressing).
    Removal of the converted source files:
    - By default nothing is removed, the source files are reported and
      left in place. The 'writeFormat' entry then decides which ones
      are used.
    - With -clean the converted source files are removed: the legacy mesh
      files (points, faces, owner, neighbour, boundary) after a conversion
      to coherent, and the coherent files (data.bp, polyMesh/coherent, the
      field stubs of parallel runs) after a conversion to ascii/binary.
      Only the selected times are affected, other times keep their files.
    - Zones, sets etc. are never removed, but note that the coherent
      format changes the face and point numbering.

    For parallel runs the legacy files are handled with the uncollated
    fileHandler unless -fileHandler is specified, regardless of the
    automatic selection for 'writeFormat coherent'.

    Ascii output uses the writePrecision and writeCompression settings of
    the controlDict.

Usage
    \b foamCoherentConvert coherent|ascii|binary [OPTION]

    Arguments:
      - \par format
        Target format: coherent, ascii or binary

    Options:
      - \par -no-mesh
        Do not convert the mesh

      - \par -no-fields
        Do not convert the fields

      - \par -clean
        Remove the converted source files (mesh and fields of the
        selected times)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "IFstream.H"
#include "OSspecific.H"
#include "cloud.H"
#include "coupledPolyPatch.H"
#include "regionProperties.H"
#include "polyMeshCoherentMapper.H"
#include "SliceStreamRepo.H"
#include "collatedFileOperation.H"
#include "uncollatedFileOperation.H"

#include "coherentConvertFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Field types supported by the coherent format
#define forAllCoherentFieldTypes(Macro)                                       \
    Macro(volScalarField);                                                    \
    Macro(volVectorField);                                                    \
    Macro(volSphericalTensorField);                                           \
    Macro(volSymmTensorField);                                                \
    Macro(volTensorField);                                                    \
    Macro(volScalarField::Internal);                                          \
    Macro(volVectorField::Internal);                                          \
    Macro(volSphericalTensorField::Internal);                                 \
    Macro(volSymmTensorField::Internal);                                      \
    Macro(volTensorField::Internal);                                          \
    Macro(surfaceScalarField);                                                \
    Macro(surfaceVectorField);                                                \
    Macro(surfaceSphericalTensorField);                                       \
    Macro(surfaceSymmTensorField);                                            \
    Macro(surfaceTensorField)


//- Objects in the global (non-processor) time directory of a region.
//  Read on the master and broadcast.
void listGlobalObjects
(
    const Time& runTime,
    const word& regionName,
    wordList& names,
    wordList& classes,
    boolList& isCoherent
)
{
    names.clear();
    classes.clear();
    isCoherent.clear();

    if (UPstream::master())
    {
        const fileName dir
        (
            runTime.globalTimePath()/polyMesh::regionName(regionName)
        );

        for (const fileName& f : Foam::readDir(dir, fileName::FILE))
        {
            IFstream is(dir/f);

            if (!is.good() || word::validate(f) != f)
            {
                continue;
            }

            IOobject io
            (
                word(f),
                runTime.timeName(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            );

            if (io.readHeader(is))
            {
                names.push_back(io.name());
                classes.push_back(io.headerClassName());
                isCoherent.push_back
                (
                    is.format() == IOstreamOption::COHERENT
                );
            }
        }
    }

    Pstream::broadcasts(UPstream::worldComm, names, classes, isCoherent);
}


//- Instance of a legacy mesh (faces file) at or before the current time,
//- empty if none
word legacyMeshInstance(const Time& runTime, const word& regionName)
{
    word inst
    (
        runTime.findInstance
        (
            polyMesh::meshDir(regionName),
            "faces",
            IOobject::READ_IF_PRESENT,
            word::null,
            false  // no constant fallback
        )
    );

    if (!returnReduceAnd(!inst.empty()))
    {
        inst.clear();
    }

    return inst;
}


//- Instance of a coherent mesh, empty if none
word coherentMeshInstance(const Time& runTime)
{
    // The search is done on the master only, which leaves the cached time
    // information of the master-based fileHandlers (collated etc.)
    // inconsistent between the ranks
    word inst(polyMeshCoherentMapper::findMeshInstance(runTime).first());
    fileHandler().flush();

    return inst;
}


//- Path of a legacy mesh file via the fileHandler (collated etc.),
//- empty if not found. Must be called on all ranks.
fileName legacyMeshFilePath
(
    const Time& runTime,
    const fileName& instance,
    const word& regionName,
    const word& name
)
{
    IOobject io
    (
        name,
        instance,
        polyMesh::meshDir(regionName),
        runTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    return fileHandler().filePath(false, io, word::null, false);
}


//- Directory holding the legacy mesh files, as seen by the fileHandler
fileName legacyMeshDirPath
(
    const Time& runTime,
    const fileName& instance,
    const word& regionName
)
{
    const fileName faces
    (
        legacyMeshFilePath(runTime, instance, regionName, "faces")
    );

    return
    (
        faces.empty()
      ? runTime.path()/instance/polyMesh::meshDir(regionName)
      : faces.path()
    );
}


//- The legacy mesh files that are converted
static const wordHashSet convertedMeshFiles
({
    "points", "faces", "owner", "neighbour", "boundary"
});


//- Remove the converted legacy mesh files. Other files (zones, sets, ...)
//- are not restored by a conversion back and are therefore kept
void removeLegacyMeshFiles
(
    const Time& runTime,
    const fileName& instance,
    const word& regionName
)
{
    for (const word& f : convertedMeshFiles)
    {
        const fileName path
        (
            legacyMeshFilePath(runTime, instance, regionName, f)
        );

        if (!path.empty())
        {
            Foam::rm(path);
        }
    }
}


//- Report files in the mesh directory that are not converted and may
//- not correspond to the converted mesh (master only)
void reportOtherMeshFiles(const Time& runTime, const fileName& dir)
{
    if (!UPstream::master())
    {
        return;
    }

    wordHashSet remaining;

    for (const fileName& f : Foam::readDir(dir, fileName::FILE))
    {
        if (!convertedMeshFiles.contains(f) && f != "coherent")
        {
            remaining.insert(f);
        }
    }

    for (const fileName& d : Foam::readDir(dir, fileName::DIRECTORY))
    {
        if (d != "data.bp")
        {
            remaining.insert(d + "/");
        }
    }

    if (remaining.size())
    {
        Warning
            << "Files in " << runTime.relativePath(dir)
            << " are not converted and may not correspond to the"
            << " face and point numbering of the converted mesh: "
            << flatOutput(remaining.sortedToc()) << nl;
    }
}


//- Warn if a legacy and a coherent mesh coexist in the same directory
void warnDuplicateMesh
(
    const word& legacyInst,
    const word& coherentInst,
    const word& regionName,
    const bool toCoherent
)
{
    if (!legacyInst.empty() && legacyInst == coherentInst)
    {
        Warning
            << "Both a legacy and a coherent mesh exist in "
            << legacyInst/polyMesh::meshDir(regionName) << nl
            << "    The " << (toCoherent ? "legacy" : "coherent")
            << " mesh is used for the conversion"
            << (toCoherent ? " and overwrites the coherent one" : "")
            << nl << endl;
    }
}


//- Warn about mesh information that the coherent format cannot hold
void warnUnsupportedMeshInfo(const polyMesh& mesh)
{
    if
    (
        returnReduceOr
        (
            mesh.pointZones().size()
         || mesh.faceZones().size()
         || mesh.cellZones().size()
        )
    )
    {
        Warning
            << "Mesh zones are not supported by the coherent format."
            << " The zone files are kept, but refer to the original"
            << " face and point numbering" << nl;
    }

    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        if (isA<coupledPolyPatch>(pp) && !isA<processorPolyPatch>(pp))
        {
            Warning
                << "Coupled patch " << pp.name() << " (type " << pp.type()
                << ") is not supported by the coherent format" << nl;
        }
    }
}


void removeStoredObjects(DynamicList<regIOobject*>& storedObjects)
{
    while (!storedObjects.empty())
    {
        storedObjects.back()->checkOut();
        storedObjects.pop_back();
    }
}


void reportSkipped(const wordHashSet& skipped, const char* reason)
{
    if (skipped.size())
    {
        Info<< "    skipped (" << reason << "): "
            << flatOutput(skipped.sortedToc()) << nl;
    }
}


void reportLagrangian(const Time& runTime, const word& regionName)
{
    const fileName dir
    (
        runTime.timePath()/polyMesh::regionName(regionName)/cloud::prefix
    );

    if (returnReduceOr(Foam::isDir(dir)))
    {
        Info<< "    skipped (not supported by coherent format): "
            << cloud::prefix << '/' << nl;
    }
}


//- Remove names that are already coherent
wordList notAlreadyCoherent
(
    const wordList& names,
    const wordHashSet& alreadyCoherent
)
{
    DynamicList<word> selected(names.size());

    for (const word& name : names)
    {
        if (!alreadyCoherent.contains(name))
        {
            selected.push_back(name);
        }
    }

    return wordList(std::move(selected));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Legacy -> coherent
void convertToCoherent
(
    Time& runTime,
    const word& regionName,
    const instantList& timeDirs,
    const IOstreamOption::streamFormat legacyFmt,
    const bool noMesh,
    const bool noFields,
    const bool clean
)
{
    const IOstreamOption coherentOpt(IOstreamOption::COHERENT);

    // Reading the legacy mesh requires a legacy format
    runTime.writeFormat(legacyFmt);

    const word legacyInst(legacyMeshInstance(runTime, regionName));
    const word coherentInst(coherentMeshInstance(runTime));

    if (!noMesh)
    {
        warnDuplicateMesh(legacyInst, coherentInst, regionName, true);
    }

    const bool legacyMesh = !legacyInst.empty();

    if (!legacyMesh)
    {
        if (!noMesh)
        {
            FatalErrorInFunction
                << "No legacy mesh found for region "
                << polyMesh::regionName(regionName) << nl
                << "Use -no-mesh if the mesh is already in coherent format."
                << nl << exit(FatalError);
        }

        Warning
            << "No legacy mesh found. Non-uniform patch values are assumed"
            << " to be already in the coherent face order" << nl;
    }

    patchPermutations perms;
    label nCellsLegacy = -1;
    fileName legacyMeshInst;

    if (legacyMesh)
    {
        Info<< "Reading legacy mesh" << nl;

        autoPtr<fvMesh> meshPtr
        (
            new fvMesh
            (
                IOobject
                (
                    regionName,
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                ),
                false
            )
        );
        meshPtr->init(true);

        fvMesh& mesh = meshPtr();

        nCellsLegacy = mesh.nCells();
        warnUnsupportedMeshInfo(mesh);
        perms.reset(mesh);

        if (!noMesh)
        {
            Info<< "Writing coherent mesh to " << mesh.facesInstance()
                << nl;

            // Re-enable writing of the mesh (disabled after reading)
            mesh.setInstance(mesh.facesInstance());

            if (!mesh.writeObject(coherentOpt, true))
            {
                FatalErrorInFunction
                    << "Failed writing coherent mesh"
                    << exit(FatalError);
            }

            legacyMeshInst = mesh.facesInstance();
        }

        // The coherent mesh is read afresh below
        meshPtr.clear();
    }

    // Remove (with -clean) or keep the converted legacy mesh
    const auto removeLegacyMesh = [&]()
    {
        if (legacyMeshInst.empty())
        {
            return;
        }

        const fileName dir
        (
            legacyMeshDirPath(runTime, legacyMeshInst, regionName)
        );

        if (clean)
        {
            Info<< nl << "Removing the converted legacy mesh files from "
                << runTime.relativePath(dir) << ": "
                << flatOutput(convertedMeshFiles.sortedToc()) << nl;

            removeLegacyMeshFiles(runTime, legacyMeshInst, regionName);
        }
        else
        {
            Info<< nl << "Keeping the converted legacy mesh files in "
                << runTime.relativePath(dir) << ": "
                << flatOutput(convertedMeshFiles.sortedToc())
                << " (use -clean to remove)" << nl;
        }

        reportOtherMeshFiles(runTime, dir);
    };

    if (noMesh && coherentInst.empty())
    {
        FatalErrorInFunction
            << "No coherent mesh found for region "
            << polyMesh::regionName(regionName) << nl
            << exit(FatalError);
    }

    if (noFields)
    {
        removeLegacyMesh();
        return;
    }

    // Reading the coherent mesh registers the CoherentMesh object
    // needed for writing coherent fields
    runTime.writeFormat(IOstreamOption::COHERENT);

    Info<< "Reading coherent mesh" << nl;

    fvMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        ),
        false
    );
    mesh.init(true);

    if (nCellsLegacy >= 0 && nCellsLegacy != mesh.nCells())
    {
        FatalErrorInFunction
            << "Coherent mesh has " << mesh.nCells()
            << " cells, but the legacy mesh has " << nCellsLegacy << nl
            << exit(FatalError);
    }

    if (legacyMesh)
    {
        perms.verify(mesh);
    }

    Info<< nl;

    forAll(timeDirs, timei)
    {
        if (timeDirs[timei].name() == runTime.constant())
        {
            continue;
        }

        runTime.setTime(timeDirs[timei], timei);
        Info<< "Time = " << runTime.timeName() << nl;

        // Fields that are already coherent (stubs in the global directory)
        wordHashSet alreadyCoherent;
        {
            wordList names, classes;
            boolList isCoherent;
            listGlobalObjects(runTime, regionName, names, classes, isCoherent);

            forAll(names, i)
            {
                if (isCoherent[i])
                {
                    alreadyCoherent.insert(names[i]);
                }
            }
        }

        // Legacy fields (in the processor directories for parallel runs)
        IOobjectList objects(mesh, runTime.timeName());

        DynamicList<regIOobject*> storedObjects;
        wordHashSet converted;

        #undef  doLocalCode
        #define doLocalCode(FieldType)                                        \
        {                                                                     \
            wordList names                                                    \
            (                                                                 \
                notAlreadyCoherent                                            \
                (                                                             \
                    objects.sortedNames<FieldType>(true),                     \
                    alreadyCoherent                                           \
                )                                                             \
            );                                                                \
            convertFields<FieldType>                                          \
            (                                                                 \
                mesh,                                                         \
                names,                                                        \
                (legacyMesh ? &perms : nullptr),                              \
                coherentOpt,                                                  \
                true,                                                         \
                storedObjects,                                                \
                converted                                                     \
            );                                                                \
        }

        forAllCoherentFieldTypes(doLocalCode);
        #undef doLocalCode

        wordHashSet skipped(objects.sortedNames(true));
        skipped -= converted;
        skipped -= alreadyCoherent;

        reportSkipped(alreadyCoherent, "already coherent");
        reportSkipped(skipped, "not supported by coherent format");
        reportLagrangian(runTime, regionName);

        removeStoredObjects(storedObjects);
    }

    SliceStreamRepo::closeInstance();

    removeLegacyMesh();
}


//- Coherent -> legacy
void convertToLegacy
(
    Time& runTime,
    const word& regionName,
    const instantList& timeDirs,
    const IOstreamOption outOpt,
    const bool noMesh,
    const bool noFields,
    const bool clean
)
{
    // Reading the coherent mesh and fields requires the coherent format
    runTime.writeFormat(IOstreamOption::COHERENT);

    const word legacyInst(legacyMeshInstance(runTime, regionName));
    const word coherentInst(coherentMeshInstance(runTime));

    warnDuplicateMesh(legacyInst, coherentInst, regionName, false);

    if (coherentInst.empty())
    {
        FatalErrorInFunction
            << "No coherent mesh found for region "
            << polyMesh::regionName(regionName) << nl
            << exit(FatalError);
    }

    Info<< "Reading coherent mesh" << nl;

    fvMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        ),
        false
    );
    mesh.init(true);

    Info<< nl;

    // Converted field names per time (for cleanup)
    List<wordHashSet> convertedPerTime(timeDirs.size());

    if (!noFields)
    {
        forAll(timeDirs, timei)
        {
            if (timeDirs[timei].name() == runTime.constant())
            {
                continue;
            }

            runTime.setTime(timeDirs[timei], timei);
            Info<< "Time = " << runTime.timeName() << nl;

            wordList names, classes;
            boolList isCoherent;
            listGlobalObjects(runTime, regionName, names, classes, isCoherent);

            DynamicList<regIOobject*> storedObjects;
            wordHashSet& converted = convertedPerTime[timei];

            // Read all fields first (some patch fields look up others)
            #undef  doLocalCode
            #define doLocalCode(FieldType)                                    \
            convertFields<FieldType>                                          \
            (                                                                 \
                mesh,                                                         \
                selectNames(FieldType::typeName, names, classes, isCoherent), \
                nullptr,                                                      \
                outOpt,                                                       \
                false,                                                        \
                storedObjects,                                                \
                converted                                                     \
            )

            forAllCoherentFieldTypes(doLocalCode);
            #undef doLocalCode

            // Close the readers before writing
            SliceStreamRepo::closeInstance();

            for (const regIOobject* obj : storedObjects)
            {
                if (!obj->writeObject(outOpt, true))
                {
                    FatalErrorInFunction
                        << "Failed writing field " << obj->name()
                        << exit(FatalError);
                }
            }

            wordHashSet skipped(names);
            skipped -= converted;
            reportSkipped(skipped, "not coherent or not supported");

            removeStoredObjects(storedObjects);
        }
    }

    if (!noMesh)
    {
        Info<< "Writing legacy mesh to " << mesh.facesInstance() << nl;

        // Re-enable writing of the mesh (disabled after reading).
        // Zones are not read from the coherent format: avoid empty files
        mesh.setInstance(mesh.facesInstance());
        mesh.pointZones().writeOpt(IOobject::NO_WRITE);
        mesh.faceZones().writeOpt(IOobject::NO_WRITE);
        mesh.cellZones().writeOpt(IOobject::NO_WRITE);

        // Existing legacy mesh files (including zones and sets) do not
        // correspond to the coherent mesh, which has a different
        // point and face ordering
        // Stale legacy mesh files (they are overwritten anyway)
        removeLegacyMeshFiles(runTime, mesh.facesInstance(), regionName);

        if (!mesh.writeObject(outOpt, true))
        {
            FatalErrorInFunction
                << "Failed writing legacy mesh"
                << exit(FatalError);
        }

        reportOtherMeshFiles
        (
            runTime,
            legacyMeshDirPath(runTime, mesh.facesInstance(), regionName)
        );

        if (UPstream::parRun())
        {
            Info<< "    Note: no procAddressing is written." << nl;
        }
    }

    // Remove (with -clean) or keep the converted coherent files
    SliceStreamRepo::closeInstance();

    // Paths of the converted coherent files (global directory)
    DynamicList<fileName> paths;

    const fileName regionDir(polyMesh::regionName(regionName));

    forAll(timeDirs, timei)
    {
        if (convertedPerTime[timei].empty())
        {
            continue;
        }

        const fileName dir
        (
            runTime.globalPath()/timeDirs[timei].name()/regionDir
        );

        if (Foam::isDir(dir/"data.bp"))
        {
            paths.push_back(dir/"data.bp");
        }

        // The stubs: overwritten by the legacy files in serial,
        // stale in parallel
        if (UPstream::parRun())
        {
            for (const word& name : convertedPerTime[timei].sortedToc())
            {
                if (Foam::isFile(dir/name))
                {
                    paths.push_back(dir/name);
                }
            }
        }
    }

    if (!noMesh)
    {
        const fileName dir
        (
            runTime.globalPath()/mesh.facesInstance()/mesh.meshDir()
        );

        for (const char* f : {"data.bp", "coherent"})
        {
            if (Foam::isDir(dir/f) || Foam::isFile(dir/f))
            {
                paths.push_back(dir/f);
            }
        }
    }

    if (paths.empty())
    {
        return;
    }

    fileNameList relative(paths.size());
    forAll(paths, i)
    {
        relative[i] = runTime.relativePath(paths[i]);
    }

    if (clean)
    {
        Info<< nl << "Removing the converted coherent files: "
            << flatOutput(relative) << nl;

        if (UPstream::master())
        {
            for (const fileName& path : paths)
            {
                if (Foam::isDir(path))
                {
                    Foam::rmDir(path);
                }
                else
                {
                    Foam::rm(path);
                }
            }
        }
    }
    else
    {
        Info<< nl << "Keeping the converted coherent files: "
            << flatOutput(relative) << " (use -clean to remove)" << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert mesh and fields between legacy and coherent formats"
    );

    timeSelector::addOptions();
    argList::addArgument
    (
        "format",
        "Target format: coherent | ascii | binary"
    );
    argList::addBoolOption
    (
        "no-mesh",
        "Do not convert the mesh"
    );
    argList::addBoolOption
    (
        "no-fields",
        "Do not convert the fields"
    );
    argList::addBoolOption
    (
        "clean",
        "Remove the converted source files (mesh and fields of the"
        " selected times)"
    );
    #include "addAllRegionOptions.H"

    // The presence of processor directories depends on the source format,
    // not on the writeFormat entry of the controlDict
    argList::noCheckProcessorDirectories();

    #include "setRootCase.H"

    const word target(args.get<word>(1));
    const bool toCoherent = (target == "coherent");

    if (!toCoherent && target != "ascii" && target != "binary")
    {
        FatalError
            << "Unknown target format '" << target
            << "', expected 'coherent', 'ascii' or 'binary'" << nl
            << exit(FatalError);
    }

    const bool noMesh = args.found("no-mesh");
    const bool noFields = args.found("no-fields");
    const bool clean = args.found("clean");

    #include "createTime.H"
    #include "getAllRegionOptions.H"

    // A parallel run with 'writeFormat coherent' automatically switches to
    // the collated fileHandler (argList). The legacy files are nevertheless
    // handled uncollated, unless the fileHandler was specified explicitly
    if
    (
        UPstream::parRun()
     && !args.found("fileHandler")
     && runTime.writeFormat() == IOstreamOption::COHERENT
     && fileHandler().type()
     == fileOperations::collatedFileOperation::typeName
    )
    {
        Info<< "Using uncollated fileHandler for the legacy files"
            << " (specify -fileHandler to override)" << nl << nl;

        fileHandler
        (
            autoPtr<fileOperation>
            (
                new fileOperations::uncollatedFileOperation(false)
            )
        );

        // Remove the (empty) processorsN directory created by argList
        if (UPstream::master())
        {
            const fileName dir
            (
                runTime.globalPath()
              / ("processors" + Foam::name(UPstream::nProcs()))
            );

            if
            (
                Foam::isDir(dir)
             && Foam::readDir(dir, fileName::FILE).empty()
             && Foam::readDir(dir, fileName::DIRECTORY).empty()
            )
            {
                Foam::rmDir(dir);
            }
        }
    }

    // Legacy format: the target, or (for reading the legacy input when
    // converting to coherent) the controlDict entry unless coherent
    IOstreamOption::streamFormat legacyFmt = runTime.writeFormat();

    if (!toCoherent)
    {
        legacyFmt = IOstreamOption::formatNames.get(target);
    }
    else if (legacyFmt == IOstreamOption::COHERENT)
    {
        legacyFmt = IOstreamOption::BINARY;
    }

    const IOstreamOption legacyOpt
    (
        legacyFmt,
        (
            legacyFmt == IOstreamOption::ASCII
          ? runTime.writeCompression()
          : IOstreamOption::UNCOMPRESSED
        )
    );


    // Time selection

    instantList timeDirs;

    if (toCoherent)
    {
        // Legacy case: normal (processor-local) time directories
        timeDirs = timeSelector::select0(runTime, args);
    }
    else
    {
        // Coherent case: time directories in the global directory only
        instantList times;

        if (UPstream::master())
        {
            const bool oldParRun = UPstream::parRun(false);
            times = TimePaths::findTimes(runTime.globalPath(), runTime.constant());
            UPstream::parRun(oldParRun);
        }
        Pstream::broadcast(times);

        timeDirs = timeSelector::select(times, args, runTime.constant());
    }

    if (timeDirs.empty())
    {
        FatalErrorInFunction
            << "No times selected" << nl
            << exit(FatalError);
    }

    // The mesh is searched backwards from the last selected time
    runTime.setTime(timeDirs.last(), timeDirs.size() - 1);

    Info<< "Converting to " << target << " format" << nl << nl;


    for (const word& regionName : regionNames)
    {
        if (regionNames.size() > 1)
        {
            Info<< "Region " << regionName << nl;
        }

        if (toCoherent)
        {
            convertToCoherent
            (
                runTime,
                regionName,
                timeDirs,
                legacyFmt,
                noMesh,
                noFields,
                clean
            );
        }
        else
        {
            convertToLegacy
            (
                runTime,
                regionName,
                timeDirs,
                legacyOpt,
                noMesh,
                noFields,
                clean
            );
        }

        Info<< nl;
    }

    Info<< "To use the converted case, set in system/controlDict:" << nl
        << "    writeFormat     " << target << ';' << nl;

    if (!toCoherent && UPstream::parRun())
    {
        Info<< "and in system/decomposeParDict:" << nl
            << "    numberOfSubdomains " << UPstream::nProcs() << ';' << nl;
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
