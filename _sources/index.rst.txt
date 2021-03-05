.. VMEC documentation master file, created by
   sphinx-quickstart on Wed Mar  3 22:18:43 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

**********************
Python Wrapper to VMEC
**********************

`VMEC <https://github.com/ORNL-Fusion/PARVMEC>`_ is an MHD equilibrium solver. 
VMEC is primarily written in Fortran and called from command line. The `python
wrapper to VMEC <https://github.com/hiddenSymmetries/VMEC2000>`_ enables users 
to call `VMEC` from python. 

Installation
============
Prerequisites
-------------
Create a virtual environment (preferably with conda) and install *numpy*.

Build & Install
---------------
At this stage, you must get the source from https://github.com/hiddenSymmetries/VMEC2000
and then do ``pip install``. 
Note: If you are using *pip*, install *numpy* also using *pip* at the prerequisites step.

Alternatively, you could install the dependencies listed in *pyproject.toml* manually and
run ``python setup.py install``.

At a later date we will provide support for ``pip install vmec``.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation

Features
========
- `VMEC python wrapper` uses a minimal VMEC source code decoupled from bigger packages such as
  `STELLOPT` and `V3FIT`. The users could compile just the Fortran code if desired.
- Provides python wrappers to parts of VMEC, particularly those needed by end-users.

Authors
=======
Caoxiang Zhu, Bharat Medasani, and Matt Landremann.

.. toctree::
   :maxdepth: 2
   :caption: API

   vmec

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
