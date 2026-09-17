# Coherent I/O format for OpenFOAM
This repository provides a short introduction and a tutorial case for usage of the coherent format in the current state.

A revised I/O method, the coherent format, has been developed within the scope of the exaFOAM project. Initial implementation in foam-extend has been recently integrated into OpenCFD’s version of OpenFOAM. The new format has numerous advantages for parallel execution:
- no need for reconstruction
- easy manipulation of boundary conditions of a field within a single file
- substantial reduction of number of files.

The implementation uses ADIOS2 software package as transport layer that enables user-defined data aggregation and subfiling.

The current state is available here (ensure that the `coherent-preview.v2412` branch is selected):
https://gitlab.com/openfoam/community/exafoam/io/-/tree/coherent-preview.v2412

Coherent I/O trainings:
- 18th OpenFOAM Workshop: https://www.youtube.com/watch?v=P12VmgSTd44
- 19th OpenFOAM Workshop: https://www.youtube.com/watch?v=gRsWkkap4os

Paper explaining the coherent mesh format:
Weiß, R.G., Lesnik, S., Galeazzo, F.C.C., Ruopp, A., Rusche, H., 2024. Coherent mesh representation for parallel I/O of unstructured polyhedral meshes. J Supercomput. https://doi.org/10.1007/s11227-024-06051-7

## Case running
- Specify “writeFormat coherent” in the controlDict dictionary.
- Create mesh, e.g. with blockMesh or snappyHexMesh or convert an existing mesh from a legacy format.
- For a serial run, the fields may be in a legacy format. User needs to make sure that cell/face order of the fields corresponds to the coherent mesh.
- For a parallel run with the coherent format, the coherent fields need to be present in the time step folders within the root directory. Fields that specify a legacy format (e.g. ascii) in the header are recognized as fields from a serial run and will cause an error. Thus, make sure that all field files indicate “coherent” format in the file headers (even if the fields consist only of uniform fields that has the same representation in all the formats). Alternatively, fields in legacy format may be supplied. These need to be located as usual in processorXX (or processorsXX for collated) directories. Again, user needs to make sure that cell/face order of the fields corresponds to the coherent mesh.

## Tutorials
Few tutorials showcasing possible workflows with the coherent format are provided. The functionality is demonstrated via Allrun scripts.
- cavity: standard lid-driven cavity case from the tutorials folder
- occDrivAerStaticMesh: a DrivAer case modified for the 1st OpenFOAM HPC Challenge (source: https://develop.openfoam.com/committees/hpc/-/tree/develop/incompressible/simpleFoam/occDrivAerStaticMesh)

The scripts are named `Allrun.<serial|parallel>.<tool>.<scenario>`; `Allrun` and `Allrun-parallel` are the plain coherent runs.
- They demonstrate different approaches for the conversion between the legacy and the coherent format:
  - `foamCoherentConvert` is the preferred approach; the `foamCoherentConvert` tool is located in this repository; simple to use, less chances of doing things wrong
  - `renumberMesh` is a hacky way and does not support conversion of coherent to legacy but is already available from the coherent-preview branch; the two-step procedure with `renumberMesh`: mesh first, then the fields, with `writeFormat` adjusted in between.

## Format converter: foamCoherentConvert
`applications/foamCoherentConvert` converts mesh and fields between the legacy formats (ascii/binary, uncollated or collated) and the coherent format, in place and in both directions. Build with the coherent OpenFOAM environment sourced: `applications/Allwmake`.
```
foamCoherentConvert coherent|ascii|binary [-no-mesh] [-no-fields] [-clean] [-time ...] [-region ...] [-parallel] [-fileHandler collated]
```
- Independent of the `writeFormat` entry of the controlDict, which has to be adjusted afterwards (a reminder is printed).
- Legacy to coherent: mesh and all selected times in one go, partitioning retained (N ranks yield N partitions; repartition with `renumberMesh -decompose`). Non-uniform patch values are permuted to the face order of the coherent format.
- Coherent to legacy (`ascii` or `binary`): a parallel run writes a decomposed case into the processor directories (without procAddressing).
- Collated cases: pass `-fileHandler collated` in both directions, for reading collated input and for writing collated output.
- Nothing is removed by default; `-clean` removes the converted source files of the selected times (never zones and sets, which do not match the renumbered faces and points of a mesh converted back).
- Converted are the mesh (without zones) and the volume, surface and volume-internal fields; other objects are reported and left untouched.

See `cavity/Allrun.*.foamCoherentConvert.*` for round trips in serial, parallel and with collated files.

## General notes
- Mesh, geometric and dimensioned fields are supported by the new format.
- All I/O data that is not supported by the coherent format yet is attempted to be handled via collated file handler. Thus, processorsXX folder may be created during runs.
- If “writeFormat coherent” is not set in controlDict, reading with the coherent format is not attempted.
- Mesh and field conversion between the legacy and coherent formats is possible using foamCoherentConvert (see above); mesh conversion and repartitioning also using renumberMesh.
- Starting a run from the fields in legacy format is supported both in serial and parallel with a condition that the mesh in coherent format is available and the cell numbering matches those of fields.

### Not supported at the moment
- Converting both the mesh and fields in one go with renumberMesh (use foamCoherentConvert).
- Running `snappyHexMesh` in parallel (probably, a path lookup issue).
- Cell-, face-, point-zones.
- Cyclic BC.
- Non-core functionality: dynamic mesh, lagrangian, finite area etc.
- Paraview reader.

## Useful commands
| Command | Description |
|-|-|
| `bpls` | A tool to view .bp files (ADIOS2) |
| `bpls constant/polyMesh/data.bp` | View contents of the mesh |
| `bpls 0.505/data.bp` | View contents from time step 0.505 |
| `bpls 0.505/data.bp -d U/internalField -c "3,3,1" -n 3` | View first 3 elements of internal field U from time step 0.505 |
| `foamCoherentConvert coherent` | Convert mesh and fields from legacy to coherent format (works also in parallel) |
| `foamCoherentConvert ascii` | Convert mesh and fields from coherent to legacy format |
| `renumberMesh -no-fields -renumber-method none -write-coherent -overwrite` | Convert mesh (works also in parallel); prerequisite: writeFormat is NOT set to coherent |
| `renumberMesh -renumber-method none` | Convert fields (works also in parallel); prerequisite: mesh available in coherent format and writeFormat is set to coherent  |
| `renumberMesh -decompose -overwrite` | First, renumber coherent mesh using the decomposition settings from system/decomposeParDict. Then, renumber each resulting partition with the default renumbering algorithm (Cuthill-McKee). Without `-overwrite` the mesh is written to the next time step and not found by a run starting from the current time. |
| `renumberMesh -decompose -renumber-method none` |  Partition coherent mesh using the decomposition settings from system/decomposeParDict without CM renumbering |
| `renumberMesh -decompose -dry-run -write-maps` | Create a VTU file (.vtu extension readable by paraview) with the coloring resulting from the decomposition and renumbering algorithms without the actual renumbering |
