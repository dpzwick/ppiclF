import unittest
import inspect
import os
from functools import wraps

###############################################################################
#  DECORATORS
###############################################################################

def Parallel(method):
    @wraps(method)
    def wrapper(self, *args, **kwargs):
        self.mpi_procs = self.parallel_procs
        if not self.ifmpi:
            self.skipTest("Skipping \"{0}\"; MPI is not enabled.".format(self.id()))
        else:
            self.log_suffix = '.general'
            if self.ifmpi:
                self.log_suffix += '.parallel'
            else:
                self.log_suffix += '.serial'
            method(self, *args, **kwargs)
    return wrapper

###############################################################################
#  BASE TEST CASE
###############################################################################

class ppiclfTestCase(unittest.TestCase):
    """ Base class for ppiclf unittests

    This defines a setUpClass method to:
        (a) get the relevant environment variables for compilers, directories
    All subclassed TestCases will need to do these things.

    Class attributes:
        f77 (str):            The Fortran 77 compiler to use     [default: 'gfortran']
        cc (str):             The C compiler to use              [default: 'gcc']
        ifmpi (bool):         Perform compilation/tests with MPI [default: False]
        source_root (str):    Path to Nek source directory;overridden by $NEK_SOURCE_ROOT env variable
                              [default: '$HOME/nek5_svn/trunk/nek']
        tools_root (str):     Path to Nek tools directory; overridden by $TOOLS_ROOT env variable
                              [default: '$HOME/nek5_svn/trunk/tools']
        examples_root (str):  Path to Nek examples directory; overridden by $EXAMPLES_ROOT env variable
                              [default: '$HOME/nek5_svn/examples']
        makenek (str):        Path to makenek                    [default: source_root/makenek]
        tools_bin (str):      Directory to place compiled tools  [default: tools_root/bin]

    Subclass attributes:
        These aren't meaningful in the base class.  They're intended for a subclass that represents
        a particular example problem.
        example_subdir (str): The subdirectory for the subclass' example.  Assumed that it's in example_root
        rea_file (str):       The .rea file for the subclass' example, minus the '.rea' extension.  Assumed
                              that it's in example_root/example_dir
        size_file (str):      The SIZE file for the subclass' example.  Assuemed that it's in
                              example_root/example_subdir
    """
    # Defined in subclasses only; declared here to make syntax checker happy
    example_subdir      = ""
    case_name           = ""

    def __init__(self, *args, **kwargs):
        # These can be overridden by self.get_opts
        self.f77            = "mpif77"
        self.cc             = "mpicc"
        self.pplist         = ""
        self.usr_lflags     = ""
        self.ifmpi          = True

        self.source_root    = os.path.dirname(os.path.dirname(inspect.getabsfile(self.__class__)))
        self.examples_root  = os.path.dirname(inspect.getabsfile(self.__class__))
        self.make           = os.path.join(self.source_root, 'Makefile')
        self.log_root       = ""
        self.verbose        = True
        self.serial_procs   = 1
        self.parallel_procs = 2
        self.size_params    = {}

        # These are overridden by method decorators (Parallel, ..)
        self.log_suffix = ""
        self.mpi_procs  = None

        # Empy list of delayed fails
        self._delayed_failures = []

        self.get_opts()

        unittest.TestCase.__init__(self, *args, **kwargs)


    def assertAlmostEqualDelayed(self, test_val, target_val, delta, label):
        if abs(test_val-target_val) <= delta:
            msg = '    SUCCESS: {0}: Test value {1} equals target value {2} +/- {3}'.format(label, test_val, target_val, delta)
        else:
            msg = '    FAILURE: {0}: Test value {1} exceeds target value {2} +/- {3}'.format(label, test_val, target_val, delta)
            self._delayed_failures.append(msg)
        print(msg)


    def assertIsNotNullDelayed(self, test_val, label):
        if test_val:
            msg = 'SUCCESS: Found phrase "{0}" in logfile.'.format(label)
        else:
            msg = 'FAILURE: Unexpectedly did not find phrase "{0}" in logfile'.format(label)
            self._delayed_failures.append(msg)
        print(msg)

    def assertIsNullDelayed(self, test_val, label):
        if test_val:
            msg = 'FAILURE: Found phrase "{0}" in logfile.'.format(label)
            self._delayed_failures.append(msg)
        else:
            msg = 'SUCCESS: Did not find phrase "{0}" in logfile'.format(label)
        print(msg)

    def assertDelayedFailures(self):
        if self._delayed_failures:
            report = [
                '\n\nFailed assertions:{0}\n'.format(len(self._delayed_failures))
            ]
            for i,failure in enumerate(self._delayed_failures, start=1):
                report.append('{0}: {1}'.format(i, failure))
            #self._delayed_failures = []
            self.fail('\n'.join(report))


    def get_opts(self):

        print("Getting setup options...")

        # Get compiler options from env
        self.f77            = os.environ.get('FC', self.f77)
        self.cc             = os.environ.get('CC', self.cc)
        self.pplist         = os.environ.get('PPLIST', self.pplist)
        self.usr_lflags     = os.environ.get('USR_LFLAGS', self.usr_lflags)
        self.ifmpi          = os.environ.get('MPI', self.ifmpi)

        # Get paths from env
        try:
            self.source_root = os.path.abspath(os.environ['SOURCE_ROOT'])
        except KeyError:
            pass
        else:
            self.make           = os.path.join(self.source_root, 'Makefile')

        self.examples_root = os.path.abspath(os.environ.get('EXAMPLES_ROOT', self.examples_root))

        try:
            self.log_root = os.path.abspath(os.environ['LOG_ROOT'])
        except KeyError:
            pass

        self.verbose        = str(os.environ.get('VERBOSE_TESTS', self.verbose)).lower() == 'true'
        self.parallel_procs = int(os.environ.get('PARALLEL_PROCS', self.parallel_procs))

        # Print everything out
        for varname, varval in (
                ('FC', self.f77),
                ('CC', self.cc),
                ('PPLIST', self.pplist),
                ('USR_LFLAGS', self.usr_lflags),
                ('IFMPI', self.ifmpi),
                ('SOURCE_ROOT', self.source_root),
                ('EXAMPLES_ROOT', self.examples_root),
                ('LOG_ROOT', self.log_root),
                ('VERBOSE_TESTS', self.verbose),
                ('PARALLEL_PROCS', self.parallel_procs)
        ):
            if varval:
                print('    Using {0:14} = "{1}"'.format(varname, varval))

        # Verify that pathnames are valid
        for varname, varval in (
                ('SOURCE_ROOT', self.source_root),
                ('EXAMPLES_ROOT', self.examples_root),
                ('LOG_ROOT', self.log_root),
        ):
            if varval and not os.path.isdir(varval):
                raise OSError('The {0} directory "{1}" does not exist. Please the env variable ${0} to a valid directory.'.format(varname, varval))

        print("Finished getting setup options!")

    def build_ppiclf(self, opts=None):
        from lib.ppiclfBinBuild import build_ppiclf
        cls = self.__class__

        all_opts = dict(
            FC = self.f77,
            CC = self.cc,
            PPLIST = self.pplist,
            USR_LFLAGS = self.usr_lflags,
            MPI = int(self.ifmpi),
        )
        if opts:
            all_opts.update(opts)

        build_ppiclf(
            source_root = self.source_root,
            cwd         = os.path.join(self.examples_root, cls.example_subdir),
            opts        = all_opts,
            verbose     = self.verbose,
        )

    def run_ppiclf(self, rea_file=None):
        from lib.ppiclfBinRun import run_ppiclf
        cls = self.__class__

        run_ppiclf(
            cwd        = os.path.join(self.examples_root, cls.example_subdir),
            rea_file   = cls.case_name if not rea_file else rea_file,
            ifmpi      = self.ifmpi,
            log_suffix = self.log_suffix,
            n_procs    = self.mpi_procs,
            verbose    = self.verbose
        )

    def get_value_from_log(self, label, column, row=0, logfile=None):
        cls = self.__class__
        if not logfile:
            logfile = os.path.join(
                self.examples_root,
                cls.example_subdir,
                '{0}.log.{1}{2}'.format(cls.case_name, self.mpi_procs, self.log_suffix)
            )
        # Get all lines with label
        with open(logfile, 'r') as f:
            line_list = [l for l in f if label in l]
        if not line_list:
            raise ValueError("Could not find label \"{0}\" in logfile \"{1}\".  The run may have failed.".format(label, logfile))
        try:
            value = float(line_list[row].split()[column])
        except ValueError:
            raise ValueError("Attempted to parse non-numerical value in logfile, \"{0}\".  Logfile may be malformatted or run may have failed".format(logfile))
        except IndexError:
            raise IndexError("Fewer rows and/or columns than expected in logfile, \"{0}\".  Logfile may be malformmated or run may have failed.".format(logfile))
        else:
            return value

    def get_phrase_from_log(self, label, logfile=None, row=0):
        cls = self.__class__
        if not logfile:
            logfile = os.path.join(
                self.examples_root,
                cls.example_subdir,
                '{0}.log.{1}{2}'.format(cls.case_name, self.mpi_procs, self.log_suffix)
            )

        with open(logfile, 'r') as f:
            line_list = [l for l in f if label in l]

        try:
            line = line_list[row]
        except IndexError:
            return None
        else:
            return line
