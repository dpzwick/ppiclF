#!/usr/bin/env python
from lib.ppiclfTestCase import *
from unittest import skip
from shutil import copyfile

import re

###############################################################################
# python -m 'unittest' ppiclfTests
#              -or-
# python -m 'unittest' ppiclfTests.testnamehere.test_Parallel
###############################################################################

class test1(ppiclfTestCase):
    example_subdir  =  'test1'
    case_name       =  'test1'

    def npart_nghost_nbin(self,cname,npt,npmn,npmx,ngt,ngmn,ngmx,nbx,nby,nbz,nbt,delt=1e-7):
          prefix_name = '(' + cname + ')'
          my_name = prefix_name + ' npart total: '
          en = self.get_value_from_log(my_name, column=-1)
          self.assertAlmostEqualDelayed(en, target_val=npt, delta=delt, label=my_name)
          my_name = prefix_name + ' npart min: '
          en = self.get_value_from_log(my_name, column=-1)
          self.assertAlmostEqualDelayed(en, target_val=npmn, delta=delt, label=my_name)
          my_name = prefix_name + ' npart max: '
          en = self.get_value_from_log(my_name, column=-1)
          self.assertAlmostEqualDelayed(en, target_val=npmx, delta=delt, label=my_name)
          my_name = prefix_name + ' nghost total: '
          en = self.get_value_from_log(my_name, column=-1)
          self.assertAlmostEqualDelayed(en, target_val=ngt, delta=delt, label=my_name)
          my_name = prefix_name + ' nghost min: '
          en = self.get_value_from_log(my_name, column=-1)
          self.assertAlmostEqualDelayed(en, target_val=ngmn, delta=delt, label=my_name)
          my_name = prefix_name + ' nghost max: '
          en = self.get_value_from_log(my_name, column=-1)
          self.assertAlmostEqualDelayed(en, target_val=ngmx, delta=delt, label=my_name)
          my_name = prefix_name + ' nbin x: '
          en = self.get_value_from_log(my_name, column=-1)
          self.assertAlmostEqualDelayed(en, target_val=nbx, delta=delt, label=my_name)
          my_name = prefix_name + ' nbin y: '
          en = self.get_value_from_log(my_name, column=-1)
          self.assertAlmostEqualDelayed(en, target_val=nby, delta=delt, label=my_name)
          my_name = prefix_name + ' nbin z: '
          en = self.get_value_from_log(my_name, column=-1)
          self.assertAlmostEqualDelayed(en, target_val=nbz, delta=delt, label=my_name)
          my_name = prefix_name + ' nbin total: '
          en = self.get_value_from_log(my_name, column=-1)
          self.assertAlmostEqualDelayed(en, target_val=nbt, delta=delt, label=my_name)

    @Parallel
    def test_Parallel(self):
         self.build_ppiclf()
         self.run_ppiclf()

         # Test A
	 self.npart_nghost_nbin(cname='A',npt=1.0,npmn=0.0,npmx=1.0,ngt=0.0,ngmn=0.0,ngmx=0.0,nbx=2.0,nby=1.0,nbz=1.0,nbt=2.0)
         my_name = '(A) TimeI error: '
         en = self.get_value_from_log(my_name, column=-1)
         self.assertAlmostEqualDelayed(en, target_val=0.0, delta=1e-7, label=my_name)
	 
         # Test B
	 self.npart_nghost_nbin(cname='B',npt=2.0,npmn=1.0,npmx=1.0,ngt=0.0,ngmn=0.0,ngmx=0.0,nbx=2.0,nby=1.0,nbz=1.0,nbt=2.0)

         # Test C
	 self.npart_nghost_nbin(cname='C',npt=2.0,npmn=1.0,npmx=1.0,ngt=2.0,ngmn=1.0,ngmx=1.0,nbx=2.0,nby=1.0,nbz=1.0,nbt=2.0)

         # Test D
	 self.npart_nghost_nbin(cname='D',npt=2.0,npmn=1.0,npmx=1.0,ngt=2.0,ngmn=1.0,ngmx=1.0,nbx=2.0,nby=1.0,nbz=1.0,nbt=2.0)

         # Test E
	 self.npart_nghost_nbin(cname='E',npt=3.0,npmn=1.0,npmx=2.0,ngt=3.0,ngmn=1.0,ngmx=2.0,nbx=2.0,nby=1.0,nbz=1.0,nbt=2.0)
         my_name = '(E) Prjct error: '
         en = self.get_value_from_log(my_name, column=-1)
         self.assertAlmostEqualDelayed(en, target_val=0.0, delta=1e-7, label=my_name)

         self.assertDelayedFailures()

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
               test1
               ) 

    suite = unittest.TestSuite([unittest.TestLoader().loadTestsFromTestCase(t) for t in testList])
    unittest.TextTestRunner(verbosity=ut_verbose, buffer=True).run(suite)
