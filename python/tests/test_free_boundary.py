#!/usr/bin/env python

import os
import shutil
import tempfile
import unittest
import logging
import numpy as np
from scipy.io import netcdf_file
from scipy.special import ellipk, ellipe
from mpi4py import MPI
import vmec

logger = logging.getLogger('[{}]'.format(MPI.COMM_WORLD.Get_rank()) + __name__)
logging.basicConfig(level=logging.INFO)

ictrl = np.zeros(5, dtype=np.int32)
verbose = False
reset_file = ''
fcomm = MPI.COMM_WORLD.py2f()

MU0 = 4.0e-7 * np.pi

# A small free-boundary tokamak held by a toroidal field and two
# poloidal-field loops. The vacuum field is analytic, so the only grid
# dependence of the result is the interpolation of the mgrid tables.
INPUT = """&INDATA
  LFREEB = T  MGRID_FILE = '{mgrid}'
  EXTCUR(1) = 1.7E+07  EXTCUR(2) = -1.7E+05
  DELT = 0.9  TCON0 = 2.0  NFP = 1  MPOL = 6  NTOR = 0  NTHETA = 24
  NS_ARRAY = 16 31  FTOL_ARRAY = 1.0E-30 1.0E-12  NITER_ARRAY = 1000 4000
  NSTEP = 200  NVACSKIP = 3  GAMMA = 0.0
  PHIEDGE = 1.2  NCURR = 1  CURTOR = 4.0E+05
  AC = 1.0 -1.0  AM = 2.0E+04 -2.0E+04
  RAXIS_CC = 1.75  ZAXIS_CS = 0.0
  RBC(0,0) = 1.7  ZBS(0,0) = 0.0
  RBC(0,1) = 0.45  ZBS(0,1) = 0.5
/
"""


def loop_field(R, Z, a, zc):
    """Field (B_R, B_Z) of a circular loop of radius a at height zc, per ampere."""
    dz = Z - zc
    den = (a + R)**2 + dz**2
    m = 4 * a * R / den
    K, E = ellipk(m), ellipe(m)
    c = MU0 / (2 * np.pi) / np.sqrt(den)
    br = c * dz / R * (-K + (a**2 + R**2 + dz**2) / ((a - R)**2 + dz**2) * E)
    bz = c * (K + (a**2 - R**2 - dz**2) / ((a - R)**2 + dz**2) * E)
    return br, bz


def write_mgrid(path, nr, nz, rmin=1.0, rmax=2.6, zmin=-0.9, zmax=0.9):
    """Scaled-mode mgrid with a 1/R toroidal field and a pair of PF loops."""
    R, Z = np.meshgrid(np.linspace(rmin, rmax, nr), np.linspace(zmin, zmax, nz))
    zero = np.zeros_like(R)
    br1, bz1 = loop_field(R, Z, 2.8, 0.8)
    br2, bz2 = loop_field(R, Z, 2.8, -0.8)
    groups = [(zero, MU0 / (2 * np.pi * R), zero), (br1 + br2, zero, bz1 + bz2)]
    f = netcdf_file(path, 'w')
    f.createDimension('stringsize', 30)
    f.createDimension('external_coil_groups', len(groups))
    f.createDimension('dim_00001', 1)
    f.createDimension('external_coils', len(groups))
    f.createDimension('rad', nr)
    f.createDimension('zee', nz)
    f.createDimension('phi', 1)
    for name, value in (('ir', nr), ('jz', nz), ('kp', 1), ('nfp', 1),
                        ('nextcur', len(groups))):
        f.createVariable(name, 'i', ())[()] = value
    for name, value in (('rmin', rmin), ('zmin', zmin), ('rmax', rmax),
                        ('zmax', zmax)):
        f.createVariable(name, 'd', ())[()] = value
    labels = f.createVariable('coil_group', 'c',
                              ('external_coil_groups', 'stringsize'))
    for i, name in enumerate(('tf', 'pf')):
        labels[i, :] = np.array(list(name.ljust(30)), dtype='S1')
    f.createVariable('mgrid_mode', 'c', ('dim_00001',))[:] = np.array([b'S'])
    f.createVariable('raw_coil_cur', 'd', ('external_coils',))[:] = 1.0
    for i, fields in enumerate(groups):
        for comp, value in zip(('br', 'bp', 'bz'), fields):
            var = f.createVariable('{}_{:03d}'.format(comp, i + 1), 'd',
                                   ('phi', 'zee', 'rad'))
            var[:] = value[None]
    f.close()


class FreeBoundaryMgridTests(unittest.TestCase):

    def run_free_boundary(self, ext, nr, nz):
        comm = MPI.COMM_WORLD
        if comm.Get_rank() == 0:
            write_mgrid('mgrid_{}.nc'.format(ext), nr, nz)
            with open('input.' + ext, 'w') as f:
                f.write(INPUT.format(mgrid='mgrid_{}.nc'.format(ext)))
        comm.Barrier()
        ictrl[:] = 0
        ictrl[0] = 1 + 2 + 4 + 8
        vmec.runvmec(ictrl, 'input.' + ext, verbose, fcomm, reset_file)
        self.assertEqual(ictrl[1], 11)
        comm.Barrier()
        f = netcdf_file('wout_{}.nc'.format(ext), mmap=False)
        raxis = np.atleast_1d(f.variables['raxis_cc'][()])[0]
        volume = f.variables['volume_p'][()]
        f.close()
        vmec.cleanup(True)
        logger.info('{}: raxis = {}, volume = {}'.format(ext, raxis, volume))
        return raxis, volume

    def test_mgrid_resolution(self):
        """
        The free-boundary solution on a 50 mm mgrid must agree with the one
        on a 12.5 mm mgrid (129 x 145 nodes), which gave raxis = 1.600552 m
        and volume = 5.358354 m^3 with cubic interpolation. Bilinear
        interpolation of the 50 mm grid moves the axis by 49 mm, cubic
        interpolation by less than 0.3 mm.
        """
        comm = MPI.COMM_WORLD
        cwd = os.getcwd()
        tmp = comm.bcast(tempfile.mkdtemp() if comm.Get_rank() == 0 else None)
        os.chdir(tmp)
        try:
            raxis, volume = self.run_free_boundary('mgrid_coarse', 33, 37)
        finally:
            os.chdir(cwd)
            comm.Barrier()
            if comm.Get_rank() == 0:
                shutil.rmtree(tmp)
        self.assertLess(abs(raxis - 1.600552), 1.0e-3)
        self.assertLess(abs(volume / 5.358354 - 1), 1.0e-3)

    def test_second_run(self):
        """A second free-boundary run in the same process gives the same result."""
        comm = MPI.COMM_WORLD
        cwd = os.getcwd()
        tmp = comm.bcast(tempfile.mkdtemp() if comm.Get_rank() == 0 else None)
        os.chdir(tmp)
        try:
            first = self.run_free_boundary('first', 33, 37)
            second = self.run_free_boundary('second', 33, 37)
        finally:
            os.chdir(cwd)
            comm.Barrier()
            if comm.Get_rank() == 0:
                shutil.rmtree(tmp)
        np.testing.assert_allclose(second, first, rtol=1.0e-12)


if __name__ == "__main__":
    unittest.main()
