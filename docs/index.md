# Overview

This repo contains both Fortran source code and Python interface for VMEC. 


# Installation of python interface (Unix like OS only)

0. Prerequisities: 
    a. Install all dependencies such as SZIP, HDF5, NETCDF (C and Fortran), and MPI 
and populate **LD_LIBRARY_PATH** with
the locations of the libraries. For BLAS and LAPACK, users are urged to use Intel MKL even for
gcc suite of compilers.

    b. If you are installating VMEC on supercomputing centers, where libraries are made avaialable 
through modules command, load all the required modules listed in prerequisites.

1. Create a python virtual environment and install pip into the virtual environment

2. Download the source code and go to root folder named `VMEC2000`. 

4. Edit the cmake_config_file.json file to guide the CMake build system to find the required libraries.
Example files are available in `python/machines` folder

5. If using Intel MKL compilers,
Define MKLROOT environment variable as
```bash
export MKLROOT=<path/to/MKL directory>
```

5. Install the python interface by running
```bash
 pip install -e .
```

6. Check if your installation worked by loading python and then
```python
import vmec
```
If the installation was successfull, the above command should not throw an error.

