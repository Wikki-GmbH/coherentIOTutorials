# Coherent I/O format for OpenFOAM
This repository provides a short introduction and a tutorial case for usage of the coherent format in the current state.

A revised I/O method, the coherent format, has been developed within the scope of the exaFOAM project. Initial implementation in foam-extend has been recently integrated into OpenCFD’s version of OpenFOAM. The new format has numerous advantages for parallel execution: no need of reconstruction, easy manipulation of boundary conditions of a field within a single file, substantial reduction of number of files. The implementation uses ADIOS2 software package as transport layer that enables user-defined data aggregation and subfiling.

Implementation is available here (ensure that the `coherent-preview.v2406` branch is selected):
https://gitlab.com/openfoam/community/exafoam/io/-/tree/coherent-preview.v2406

Coherent I/O training from the 18th OpenFOAM Workshop:
https://www.youtube.com/watch?v=P12VmgSTd44

Paper explaining the coherent mesh format:
Weiß, R.G., Lesnik, S., Galeazzo, F.C.C., Ruopp, A., Rusche, H., 2024. Coherent mesh representation for parallel I/O of unstructured polyhedral meshes. J Supercomput. https://doi.org/10.1007/s11227-024-06051-7

## Case running
- Specify “writeFormat coherent” in the controlDict dictionary.
- Create mesh (e.g. blockMesh or snappyHexMesh).
- For a serial run, the fields may be in a legacy format. User needs to make sure that cell/face order of the fields corresponds to the coherent mesh.
- For a parallel run with the coherent format, the coherent fields need to be present in the time step folders within the root directory. Fields that specify a legacy format (e.g. ascii) in the header are recognized as fields from a serial run and will cause an error. Thus, make sure that all field files indicate “coherent” format in the file headers (even if the fields consist only of uniform fields that has the same representation in all the formats). Alternatively, fields in legacy format may be supplied. These need to be located as usual in processorXX (or processorsXX for collated) directories. Again, user needs to make sure that cell/face order of the fields corresponds to the coherent mesh.

## General notes
- Mesh, geometric and dimensioned fields are supported by the new format. Other data structures will try to use a fallback, which utilizes the collated file handler.
- To some extent, format mixing is supported.
- If “writeFormat coherent” is not set in controlDict, reading with the coherent format is not attempted.

### Not supported at the moment
- Mesh conversion
- Pre-processing tools not tested
- Field conversion from the coherent to a legacy format

### Convert from a legacy (ascii, binary, collated) to the coherent format
By now, only fields may be converted. This means that you will need to create a mesh in the coherent format that corresponds to the original mesh. See dedicated Allrun scripts in the tutorials provided here.

## Useful commands
| Command | Description |
|-|-|
| `bpls` | A tool to view .bp files (ADIOS2) |
| `bpls constant/polyMesh/data.bp` | View contents of the mesh |
| `bpls 0.505/data.bp` | View contents from time step 0.505 |
| `bpls 0.505/data.bp -d U/internalField -c "3,3,1" -n 3` | View first 3 elements of internal field U from time step 0.505 |
| `renumberMesh -decompose` | First, renumber coherent mesh using the decomposition settings from system/decomposeParDict. Then, renumber each resulting partition with the default renumbering algorithm (Cuthill-Mckee). |
| `renumberMesh -decompose -renumber-method none` |  Renumber coherent mesh using the decomposition settings from system/decomposeParDict |
| `renumberMesh -decompose -dry-run -write-maps` | Create a VTU file (.vtu extension readable by paraview) with the coloring resulting from the decomposition and renumbering algorithms without the actual renumbering |