import unittest
import numpy as np
from mpi4py import MPI
import os
import vmec

reset_file = ''

class F90wrapVmecTests(unittest.TestCase):

    def test_read_input(self):
        """
        Try reading a VMEC input file.
        """
        self.ictrl[0] = run_modes['input']
        vmec_f90wrap.runvmec(self.ictrl, self.filename, self.verbose, \
                                 self.fcomm, reset_file)

        self.assertTrue(self.ictrl[1] in success_codes)

        self.assertEqual(vmec_f90wrap.vmec_input.nfp, 3)
        self.assertEqual(vmec_f90wrap.vmec_input.mpol, 4)
        self.assertEqual(vmec_f90wrap.vmec_input.ntor, 3)
        print('rbc.shape:', vmec_f90wrap.vmec_input.rbc.shape)
        print('rbc:',vmec_f90wrap.vmec_input.rbc[101:103, 0:4])

        # n = 0, m = 0:
        self.assertAlmostEqual(vmec_f90wrap.vmec_input.rbc[101,0], 1.3782)

        # n = 0, m = 1:
        self.assertAlmostEqual(vmec_f90wrap.vmec_input.zbs[101,1], 4.6465E-01)

        # n = 1, m = 1:
        self.assertAlmostEqual(vmec_f90wrap.vmec_input.zbs[102,1], 1.6516E-01)



    def test_run_read(self):
        """
        Try running VMEC, then reading in results from the wout file.
        """

        self.ictrl[0] = 1 + 2 + 4 + 8
        vmec_f90wrap.runvmec(self.ictrl, self.filename, self.verbose, \
                                 self.fcomm, reset_file)

        self.assertTrue(self.ictrl[1] in success_codes)

        self.assertEqual(vmec_f90wrap.vmec_input.nfp, 3)
        self.assertEqual(vmec_f90wrap.vmec_input.mpol, 4)
        self.assertEqual(vmec_f90wrap.vmec_input.ntor, 3)
        print('rbc.shape:', vmec_f90wrap.vmec_input.rbc.shape)
        print('rbc:',vmec_f90wrap.vmec_input.rbc[101:103, 0:4])

        # n = 0, m = 0:
        self.assertAlmostEqual(vmec_f90wrap.vmec_input.rbc[101,0], 1.3782)

        # n = 0, m = 1:
        self.assertAlmostEqual(vmec_f90wrap.vmec_input.zbs[101,1], 4.6465E-01)

        # n = 1, m = 1:
        self.assertAlmostEqual(vmec_f90wrap.vmec_input.zbs[102,1], 1.6516E-01)

        # Now try reading in the output
        wout_file = os.path.join(os.path.dirname(__file__), 'wout_li383_low_res.nc')
        ierr = 0
        vmec_f90wrap.read_wout_mod.read_wout_file(wout_file, ierr)
        self.assertEqual(ierr, 0)
        self.assertAlmostEqual(vmec_f90wrap.read_wout_mod.betatot, \
                                   0.0426215030653306, places=4)

        print('iotaf.shape:',vmec_f90wrap.read_wout_mod.iotaf.shape)
        print('rmnc.shape:',vmec_f90wrap.read_wout_mod.rmnc.shape)

        self.assertAlmostEqual(vmec_f90wrap.read_wout_mod.iotaf[-1], \
                                   0.654868168783638, places=4)

        self.assertAlmostEqual(vmec_f90wrap.read_wout_mod.rmnc[0, 0], \
                                   1.4773028173065, places=4)


if __name__ == "__main__":
    unittest.main()
