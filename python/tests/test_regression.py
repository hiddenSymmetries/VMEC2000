#!/usr/bin/env python

import unittest
import logging
import os
import numpy as np
from scipy.io import netcdf_file
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

def compare(f1, f2, field):
    """
    Assert that the difference between a field of the new output file
    and the reference file is sufficiently small. Also print the
    differences to the log.
    """
    x1 = f1.variables[field][()]
    x2 = f2.variables[field][()]
    logger.info('Diff in field {} is {}'.format(field, np.max(np.abs(x2 - x1))))
    np.testing.assert_allclose(x1, x2, atol=1e-10, rtol=1e-6)

    
class RegressionTests(unittest.TestCase):

    def test_regression(self):
        """
        Read in multiple input files, running VMEC various numbers of
        times in between.  Several output quantities are compared to
        known reference values.  This test case covers both cases in
        which VMEC converges and cases in which it does not, due to an
        unachievably small ftol.  This test case includes several
        geometries, covering both axisymmetry and non-axisymmetry, and
        it covers both stellarator-symmetry and
        non-stellarator-symmetry.
        """
        files = ['circular_tokamak',
                 'up_down_asymmetric_tokamak',
                 'li383_low_res',
                 'LandremanSenguptaPlunk_section5p3_low_res']

        #file_indices = list(range(len(files)))

        for nrun in range(4):
            for jfile in range(len(files)):
                filename = os.path.join(os.path.dirname(__file__), 'input.' + files[jfile])

                # Read in the input file:
                ictrl[:] = 0
                ictrl[0] = restart_flag + readin_flag
                logger.info("Calling runvmec to read input file {}. ictrl={} comm={}" \
                            .format(filename, ictrl, fcomm))
                vmec.runvmec(ictrl, filename, verbose, fcomm, reset_file)
                self.assertEqual(ictrl[1], 0)
                ftol_array_orig = np.copy(vmec.vmec_input.ftol_array)

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

                    # Cover both cases in which VMEC converges and
                    # cases in which it does not. To do this, for
                    # every other run we set ftol to a value so small
                    # that it is unattainable.
                    should_converge = (np.mod(jrun, 2) == 0)
                    vmec.vmec_input.ftol_array = ftol_array_orig
                    if not should_converge:
                        vmec.vmec_input.ftol_array[:3] = 1.0e-30
                    
                    logger.info("About to reinit. nrun={} jfile={} jrun={} should_converge={}" \
                                .format(nrun, jfile, jrun, should_converge))
                    logger.info("ftol_array_orig: {}".format(ftol_array_orig[:5]))
                    logger.info("ftol_array: {}".format(vmec.vmec_input.ftol_array[:5]))
                    vmec.reinit()

                    # Time-step:
                    ictrl[:] = 0
                    ictrl[0] = restart_flag + reset_jacdt_flag + timestep_flag + output_flag
                    logger.info("Calling runvmec to timestep. ictrl={} comm={}".format(ictrl, fcomm))
                    vmec.runvmec(ictrl, filename, verbose, fcomm, reset_file)

                    # Expected return codes for ictrl[1], from vmec_params.f:
                    # 2 = more_iter_flag
                    # 11 = successful_term_flag
                    if should_converge:
                        self.assertEqual(ictrl[1], 11)
                    else:
                        self.assertEqual(ictrl[1], 2)

                    # Regression tests: compare output to reference files.
                    if should_converge:
                        # Before trying to read the wout file, make sure master is
                        # done writing it.
                        MPI.COMM_WORLD.Barrier()

                        # New wout output file is in the current
                        # working directory, whereas the reference
                        # wout file is in a fixed directory.
                        wout_file = 'wout_' + files[jfile] + '.nc'
                        reference_file = os.path.join(os.path.dirname(__file__), 'wout_' + files[jfile] + '_reference.nc')

                        ierr = 0
                        logger.info('About to read output file {}'.format(wout_file))
                        f1 = netcdf_file(wout_file, mmap=False)
                        #vmec.read_wout_mod.read_wout_file(wout_file, ierr)
                        #self.assertEqual(ierr, 0)

                        logger.info('About to read reference output file {}'.format(reference_file))
                        f2 = netcdf_file(reference_file, mmap=False)

                        compare(f1, f2, 'iotaf')
                        compare(f1, f2, 'rmnc')
                        compare(f1, f2, 'zmns')
                        compare(f1, f2, 'lmns')
                        compare(f1, f2, 'bmnc')
                        #np.testing.assert_allclose(vmec.read_wout_mod.iotaf, f2.variables['iotaf'][()])
                        
                        f1.close()
                        f2.close()
                        
                        
                    # Deallocate arrays:
                    logger.info("About to cleanup(True). nrun={} jfile={} jrun={}" \
                                .format(nrun, jfile, jrun))
                    vmec.cleanup(True)
                    

if __name__ == "__main__":
    unittest.main()
