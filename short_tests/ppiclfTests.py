#!/usr/bin/env python
from lib.ppiclfTestCase import *
from unittest import skip
from shutil import copyfile

import re

###############################################################################

class stokes_2d(ppiclfTestCase):
    example_subdir  =  'stokes_2d'
    case_name       =  'stokes_2d'

#    def setUp(self):
#        # nothing for now


    @Parallel
    def test_Parallel(self):
#        # self.config_size()
         self.build_ppiclf()
         self.run_ppiclf()
#
#        # en = self.get_value_from_log('proj error: ', column=-1)
#        # self.assertAlmostEqualDelayed(en, target_val=3.434015E-08, delta=1e-7, label='proj error')
#
#        # self.assertDelayedFailures()

###############################################################

if __name__ == '__main__':
    import argparse, os

    # Get arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--f77", default='mpif77', help="The Fortran 77 compiler to use [default: mpif77]")
    parser.add_argument("--cc", default='mpicc',  help="The C compiler to use [default: mpicc]")
    parser.add_argument("--ifmpi", default='true', choices=['true', 'false'], help="Enable/disable parallel tests with MPI [default: true]")
    parser.add_argument("--nprocs", default='4', help="Number of processes to use for MPI tests [default: 4]")
    parser.add_argument("-v", "--verbose", action='store_true', help="Enable verbose output")
 
    args = parser.parse_args()

    # Set environment
    os.environ['CC'] = args.cc
    os.environ['FC'] = args.f77
    os.environ['IFMPI'] = args.ifmpi
    os.environ['PARALLEL_PROCS'] = args.nprocs
    if args.verbose:
        os.environ['VERBOSE_TESTS'] = 'true'
        ut_verbose = 2
    else:
        os.environ['VERBOSE_TESTS'] = 'false'
        ut_verbose = 1

    testList = (
               stokes_2d
               ) 

    suite = unittest.TestSuite([unittest.TestLoader().loadTestsFromTestCase(t) for t in testList])
    unittest.TextTestRunner(verbosity=ut_verbose, buffer=True).run(suite)
