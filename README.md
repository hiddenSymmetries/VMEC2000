## Overview

This VMEC repo is to build VMEC as isolated piece of software. It also comes with support for building a python extension,
so VMEC can be run directly from python, with all variables passed directly in memory.

## Downloading
Download the package from git. And change to the root directory of VMEC source code by running
```bash
cd <VMEC_ROOT>
```

## Compiling

Compiling VMEC requires MPI, parallel version of HDF5, netcdf-c, netcdf-fortran, numerical libraries such as BLAS, LAPACK, BLACS, and SCALAPACK.

If you want to build the python extension, skip down to the "Python extension compiling" section below. Building the python extension will also build the VMEC executable.

### PPPL Cluster
On PPPL cluster, the dependencies can be enabled using `module` command as
```bash
module load intel openmpi blacs scalapack netcdf-c netcdf-fortran
```
The location of the libraries associated with the packages can be found from
```bash
echo $LD_LIBRARY_PATH
```

(Optional) If interested in using cmake, load cmake also using module command.

The one could use either `make` or `cmake` to build VMEC. 
1. Make
-------
If using make, edit the `Makefile` to correctly reflect the location of the packages and the compilers. Then run 
```bash
make all
```

2. CMake
--------
The following command was used to configure cmake build setup for VMEC on PPPL cluster
```bash
cmake -S. -Bbuild -GNinja -DCMAKE_Fortran_COMPILER=mpifort -DNETCDF_INC_PATH=/usr/pppl/intel/2019-pkgs/netcdf-fortran-4.5.2/include/ -DNETCDF_LIB_PATH=/usr/pppl/intel/2019-pkgs/netcdf-fortran-4.5.2/lib -DSCALAPACK_LIB_NAME=scalapack -DBLACS_LIB_NAME=mpiblacs -DBLACS_CINIT_NAME=mpiblacsCinit -DBLACS_F77INIT_NAME=mpiblacsF77init --trace-source=CMakeLists.txt 2>&1 | tee log
```
There are few points to note on the above command
  - The above command gives a verbose output and also store the output in `log` file. 
  - Ninja build system is used. If your system doesn't have ninja installed, remove the -G option. The default is Makefile setup. Alternatively you can install ninja in a python virtual environment using pip.

After successful completion of cmake configuration step, go to `build` directory and run `ninja` command as
```bash
cd build; ninja
````

## Python Extension Compiling
In the VMEC root folder, edit the `cmake_config_file.json` as necessary for your system. Several example `.json` files are provided in the `cmake/machines` folder. The default one works for CentOS 8 with netcdf and openmpi installed from CentOS repos and Intel MKL.

There are two ways to build the python extension. 
1. You could simply  run
```bash
pip install .
```
Sometimes you may get error associated with f90wrap installation. In that case, install numpy before running the above command:
```bash
pip install numpy
pip install .
```
If you are using conda virtual environment, the numpy version used by pip to build the VMEC extension could be different from the numpy version installed in the conda virtaul environment. The numpy version mismatch could pose some issues. For such cases, use the second approach given below.

2. In this approach, we directly run the setup.py script present in the VMEC root folder. Ensure that the following modules are installed: `cmake`, `scikit-build`, `ninja`, `numpy`, and `f90wrap`. These packages could be installed by running the following commands:
```bash
pip install numpy
pip install cmake scikit-build ninja f90wrap
```
Some times you may have to install `numpy` before trying to install `f90wrap`.

Then, run 
```bash
python setup.py build_ext
python setup.py install
``` 
in succession. At this point, you should be able to import the `vmec` module in python. To test this, you can try the following command from the shell:
```bash
python -c "import vmec; print('success')"
```

