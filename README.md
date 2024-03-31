# SOP-MULTI
SOP-MULTI is a LAMMPS package for simulating IDP and Multidomain Proteins


Usage instructions:
Place all the files and folders from lammps-SOP-MULTI-src in the src folder of the LAMMPS source code and use <br/>
$ make yes-SOP-MULTI <br/>
$ make yes-OPENMP (optional, openmp acceleration. Note: No significant enhancement in the runtime was observed with single molecule systems)

## In Serial environments
$ make serial
## In Parallel environments
$ make mpi

Once successfully compiled you will notice the executable lmp_serial, lmp_mpi or lmp_omp. This executable can be used along with the input files
to run the simulations.

This Package has been tested with lammps-2Aug2023 version
