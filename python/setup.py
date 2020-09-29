import sys

import setuptools
from skbuild import setup

setup(
    name="vmec",
    version="0.0.1",
    packages=["vmec"],
    install_requires=['f90wrap == v0.2.3'],
    cmake_args=[
           '-DCMAKE_Fortran_COMPILER=mpifort',
           '-DNETCDF_INC_PATH=/usr/pppl/intel/2019-pkgs/netcdf-fortran-4.5.2/include',
           '-DNETCDF_LIB_PATH=/usr/pppl/intel/2019-pkgs/netcdf-fortran-4.5.2/lib',
           '-DSCALAPACK_LIB_NAME=scalapack',
           '-DBLACS_LIB_NAME=mpiblacs',
           '-DBLACS_CINIT_NAME=mpiblacsCinit',
           '-DBLACS_F77INIT_NAME=mpiblacsF77init'],
    cmake_source_dir=".."
)
