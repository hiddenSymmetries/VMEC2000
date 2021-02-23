#!/usr/bin/env python

import unittest
import logging
import os
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

class RerunTests(unittest.TestCase):

    def test_rerun(self):
        """
        Read in multiple input files, running VMEC various numbers of times in between.
        """
        files = ['circular_tokamak',
                 'up_down_asymmetric_tokamak',
                 'li383_low_res',
                 'LandremanSenguptaPlunk_section5p3_low_res']

        #file_indices = list(range(len(files)))

        for nrun in range(3):
            for jfile in range(len(files)):
                filename = os.path.join(os.path.dirname(__file__), 'input.' + files[jfile])

                # Read in the input file:
                ictrl[:] = 0
                ictrl[0] = restart_flag + readin_flag
                logger.info("Calling runvmec to read input file {}. ictrl={} comm={}" \
                            .format(filename, ictrl, fcomm))
                vmec.runvmec(ictrl, filename, verbose, fcomm, reset_file)
                self.assertEqual(ictrl[1], 0)

                # Deallocate arrays:
                logger.info("About to cleanup(False)")
                vmec.cleanup(False) # False because we have not time-stepped.

                for jrun in range(nrun):
                    # Set axis to 0 to force VMEC to recompute its
                    # initial guess for the axis. This way each run of
                    # the code should produce independent results, no
                    # matter how many times the code is run.
                    vmec.vmec_input.raxis_cc = 0
                    vmec.vmec_input.raxis_cs = 0
                    vmec.vmec_input.zaxis_cc = 0
                    vmec.vmec_input.zaxis_cs = 0
                    
                    logger.info("About to reinit. nrun={} jfile={} jrun={}" \
                                .format(nrun, jfile, jrun))
                    vmec.reinit()

                    # Time-step:
                    ictrl[:] = 0
                    ictrl[0] = restart_flag + reset_jacdt_flag + timestep_flag + output_flag
                    logger.info("Calling runvmec to timestep. ictrl={} comm={}".format(ictrl, fcomm))
                    vmec.runvmec(ictrl, filename, verbose, fcomm, reset_file)
                    self.assertTrue(ictrl[1] in [0, 11]) # 11 if converged, 0 if not.

                    # Deallocate arrays:
                    logger.info("About to cleanup(True). nrun={} jfile={} jrun={}" \
                                .format(nrun, jfile, jrun))
                    vmec.cleanup(True)
                    

if __name__ == "__main__":
    unittest.main()
