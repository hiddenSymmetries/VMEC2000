import sys

import setuptools
from skbuild import setup

setup(
    name="vmec",
    version="0.0.1",
    packages=["vmec"],
    install_requires=['f90wrap == v0.2.3'],
    cmake_args=[
	       '-DCMAKE_C_COMPILER=mpicc',
	       '-DCMAKE_CXX_COMPILER=mpicxx',
           '-DCMAKE_Fortran_COMPILER=mpifort',
           #'-DNETCDF_INC_PATH=/usr/pppl/intel/2019-pkgs/netcdf-fortran-4.5.2/include',
           #'-DNETCDF_LIB_PATH=/usr/pppl/intel/2019-pkgs/netcdf-fortran-4.5.2/lib',
           '-DNETCDF_INC_PATH=/usr/pppl/gcc/9.3-pkgs/netcdf-fortran-4.5.2/include',
           '-DNETCDF_LIB_PATH=/usr/pppl/gcc/9.3-pkgs/netcdf-fortran-4.5.2/lib',
           '-DBLAS_LIBRARIES=/usr/pppl/gcc/9.3-pkgs/lapack-3.9.0/lib/librefblas.a',
           '-DLAPACK_LIBRARIES=/usr/pppl/gcc/9.3-pkgs/lapack-3.9.0/lib/liblapack.a',
           '-DSCALAPACK_LIB_NAME=scalapack',
           '-DBLACS_LIB_NAME=mpiblacs',
           '-DBLACS_CINIT_NAME=mpiblacsCinit',
           '-DBLACS_F77INIT_NAME=mpiblacsF77init'],
    cmake_source_dir=".."
)
