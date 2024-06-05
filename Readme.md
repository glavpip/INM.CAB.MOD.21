INMOM version with modified advection scheme based on CABARET

This repository contains the code of the OGCM, which is a subject of the paper "Gusev, A.V. A computational complex for numerical simulation of the North Atlantic dynamics in eddy-resolving mode" submitted to "Computers & Geosciences".

Content of the repository.

!----------------------------------------------------------------------------------------------------------------
PETSC

The subfolder petsc-3.18.5 contains an external library for solving linear equation sets (in our case, shallow water equations). To build the library for user's platform, one needs to enter the directory petsc-3.18.5 anr execute ./configure . In the case of failure, the compilers and other keys should be inserted manually. For instance,

./configure PETSC_ARCH=centos-intelmpi --prefix=petsc-lib --with-cc=mpiicc --with-cxx=mpiicpc --with-fc=mpiifort --with-blas-lapack-dir=/opt/intel/compilers_and_libraries_2019.5.281/linux/mkl

Then follow installer's appointments and store the paths to $HOME/.bashrc

!----------------------------------------------------------------------------------------------------------------
OGCM

The subfolder OGCM contains source files of OGCM itself without PETSC. To build the OGCM, one needs to build PETSC first.
Then enter the folder OGCM and run 'make all'. The executable inmsom.exe is to be created.

OGCM/Control - folder with subroutines controlling the OGCM operation.
OGCM/Function - folder with subroutines containing solution of individual physical subproblems.
OGCM/Inc - folder with include-files containing domain configuration, physical and technical parameters.
OGCM/INMSOM - folder with project files for Microsoft Visual Studio with PGI Fortran.
OGCM/Modules - folder with declaration of arrays within fortran modules.
OGCM/Service - folder with service procedures.
OGCM/inmsom_head.f90 - main head files of the OGCM.
OGCM/Makefile.inc - file with parameters to be used within make procedure.
OGCM/makefike - makefile itself.
OGCM/ocean_run.par - file with task parameters for run.
OGCM/phys_proc.par - file with physical parameters to run.
OGCM/run_petsc.qs - an example of a script for running the task in the framework of the Slurm resource manager.

!---------------------------------------------------------------------------------------------------------------
DATA

The full datasets needed for performing 1 year run takes more than 30Gbyte, what is out of Github functionality. In the case of necessity, these data can be provided on request by using external services.
