#!/usr/bin/env python3

import json
import argparse

import setuptools
from skbuild import setup

#parser = argparse.ArgumentParser(
#        description="Process arguments to building VMEC python interface")
#parser.add_argument('--cmake_config_file', 
#                    default="pppl_gcc.json")
#args = parser.parse_known_args()
#print(args)

#with open('pppl_intel.json') as fp:
with open('cmake_config_file.json') as fp:
    d = json.load(fp)

setup(
    name="vmec",
    version="0.0.1",
    packages=["vmec"],
    install_requires=['f90wrap == v0.2.3'],
    cmake_args=d['cmake_args'],
           # [
           #'-DCMAKE_C_COMPILER=mpicc',
           #'-DCMAKE_CXX_COMPILER=mpicxx',
           #'-DCMAKE_Fortran_COMPILER=mpifort',
           ##'-DNETCDF_INC_PATH=/usr/pppl/intel/2019-pkgs/netcdf-fortran-4.5.2/include',
           ##'-DNETCDF_LIB_PATH=/usr/pppl/intel/2019-pkgs/netcdf-fortran-4.5.2/lib',
           #'-DNETCDF_INC_PATH=/usr/pppl/gcc/9.3-pkgs/netcdf-fortran-4.5.2/include',
           #'-DNETCDF_LIB_PATH=/usr/pppl/gcc/9.3-pkgs/netcdf-fortran-4.5.2/lib',
           ##'-DBLAS_LIBRARIES=/usr/pppl/gcc/9.3-pkgs/lapack-3.9.0/lib/librefblas.a',
           ##'-DLAPACK_LIBRARIES=/usr/pppl/gcc/9.3-pkgs/lapack-3.9.0/lib/liblapack.a',
           #'-DBLA_VENDOR=Intel10_64lp',
           #'-DSCALAPACK_LIB_NAME=scalapack',
           #'-DBLACS_LIB_NAME=mpiblacs',
           #'-DBLACS_CINIT_NAME=mpiblacsCinit',
           #'-DBLACS_F77INIT_NAME=mpiblacsF77init'],
    cmake_source_dir=".."
)
