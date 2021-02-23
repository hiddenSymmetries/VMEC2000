#!/usr/bin/env python

import unittest
import os
import logging
import numpy as np
from mpi4py import MPI
import vmec

logger = logging.getLogger('[{}]'.format(MPI.COMM_WORLD.Get_rank()) + __name__)
logging.basicConfig(level=logging.INFO)

ictrl = np.zeros(5, dtype=np.int32)
verbose = True
reset_file = ''

# Flags used by runvmec():
restart_flag = 1
readin_flag = 2
timestep_flag = 4
output_flag = 8
cleanup_flag = 16
reset_jacdt_flag = 32

fcomm = MPI.COMM_WORLD.py2f()

filename = os.path.join(os.path.dirname(__file__), 'input.li383_low_res')

class VmecTests(unittest.TestCase):

    def test_read_input(self):
        """
        Try reading a VMEC input file.
        """
        ictrl[:] = 0
        ictrl[0] = restart_flag + readin_flag
        logger.info("Calling runvmec. ictrl={} comm={}".format(ictrl, fcomm))
        vmec.runvmec(ictrl, filename, verbose, fcomm, reset_file)

        self.assertEqual(ictrl[1], 0)

        self.assertEqual(vmec.vmec_input.nfp, 3)
        self.assertEqual(vmec.vmec_input.mpol, 4)
        self.assertEqual(vmec.vmec_input.ntor, 3)
        print('rbc.shape:', vmec.vmec_input.rbc.shape)
        print('rbc:',vmec.vmec_input.rbc[101:103, 0:4])

        # n = 0, m = 0:
        self.assertAlmostEqual(vmec.vmec_input.rbc[101,0], 1.3782)

        # n = 0, m = 1:
        self.assertAlmostEqual(vmec.vmec_input.zbs[101,1], 4.6465E-01)

        # n = 1, m = 1:
        self.assertAlmostEqual(vmec.vmec_input.zbs[102,1], 1.6516E-01)

        logger.info("About to cleanup(False)")
        vmec.cleanup(False)

    def test_run_read(self):
        """
        Try running VMEC, then reading in the output.
        """
        ictrl[:] = 0
        ictrl[0] = 1 + 2 + 4 + 8
        vmec.runvmec(ictrl, filename, verbose, fcomm, reset_file)

        self.assertEqual(ictrl[1], 11)

        self.assertEqual(vmec.vmec_input.nfp, 3)
        self.assertEqual(vmec.vmec_input.mpol, 4)
        self.assertEqual(vmec.vmec_input.ntor, 3)
        print('rbc.shape:', vmec.vmec_input.rbc.shape)
        print('rbc:',vmec.vmec_input.rbc[101:103, 0:4])

        # n = 0, m = 0:
        self.assertAlmostEqual(vmec.vmec_input.rbc[101,0], 1.3782)

        # n = 0, m = 1:
        self.assertAlmostEqual(vmec.vmec_input.zbs[101,1], 4.6465E-01)

        # n = 1, m = 1:
        self.assertAlmostEqual(vmec.vmec_input.zbs[102,1], 1.6516E-01)

        # Now try reading in the output
        wout_file = os.path.join(os.path.dirname(__file__), 'wout_li383_low_res.nc')
        ierr = 0
        vmec.read_wout_mod.read_wout_file(wout_file, ierr)
        self.assertEqual(ierr, 0)
        self.assertAlmostEqual(vmec.read_wout_mod.betatot, \
                                   0.0426215030653306, places=4)

        print('iotaf.shape:',vmec.read_wout_mod.iotaf.shape)
        print('rmnc.shape:',vmec.read_wout_mod.rmnc.shape)

        self.assertAlmostEqual(vmec.read_wout_mod.iotaf[-1], \
                                   0.654868168783638, places=4)

        self.assertAlmostEqual(vmec.read_wout_mod.rmnc[0, 0], \
                                   1.4773028173065, places=4)
        
        logger.info("About to cleanup(True)")
        vmec.cleanup(True)

if __name__ == "__main__":
    unittest.main()
