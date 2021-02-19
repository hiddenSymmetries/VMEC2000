#!/usr/bin/env python3
import sys
import json
import argparse
import setuptools
import pathlib
import os

from os.path import basename, splitext
from glob import glob

print("system.platform is {}".format(sys.platform))
if (sys.platform == "darwin"):
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-dynamiclib')

from skbuild import setup

fldr_path = pathlib.Path(__file__).parent.absolute()

with open(os.path.join(fldr_path, 'cmake_config_file.json')) as fp:
    d = json.load(fp)

class EmptyListWithLength(list):
    def __len__(self):
        return 1

setup(
    name="vmec",
    version="0.0.2",
    license="MIT",
    url="https://gitlab.com/mbkumar/VMEC2000",
    packages=['vmec'],
    package_dir={'': 'src'},
    #py_modules=[splitext(basename(path))[0] for path in glob('src/vmec/*.py')],
    install_requires=['f90wrap == v0.2.3'],
    python_requires=">=3.7",
    ext_modules=EmptyListWithLength(),
    description="Python wrapper for VMEC2000",
    maintainer="Bharat Medasani",
    classifiers=[
        "Development Status :: 1 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Information Technology",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Nuclear Fusion Community",
        "Operating System :: MacOS",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Fortran",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Nuclear Fusion",
        "Topic :: Scientific/Engineering :: Nuclear Fusion :: Stellarator Research",
        "Topic :: Scientific/Engineering :: Optimization",
        "Topic :: Software Development",],
    cmake_args=d['cmake_args'],
    cmake_source_dir=os.path.abspath("..")
)
