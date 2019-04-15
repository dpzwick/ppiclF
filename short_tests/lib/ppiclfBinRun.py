import os
import sys
from subprocess import call, check_call, PIPE, STDOUT, Popen, CalledProcessError

def run_ppiclf(cwd, rea_file, ifmpi, log_suffix='', n_procs=1, verbose=False):
    # Paths to executables, files
    test         = os.path.join(cwd, 'test.out')
    logfile      = os.path.join(cwd, '{0}.log.{1}{2}'.format(rea_file, n_procs, log_suffix))
    if ifmpi:
        command = ['mpiexec', '-np', str(n_procs), test]
    else:
        command = [test]

    print("Running test...")
    print('    Using command "{0}"'.format(' '.join(command)))
    print('    Using working directory "{0}"'.format(cwd))

    # Any error here is unexepected
    try:
        if verbose:
            with open(logfile, 'w') as f:
                proc =Popen(command, cwd=cwd, stderr=STDOUT, stdout=PIPE)
                for line in proc.stdout:
                    sys.stdout.write(line)
                    f.write(line)
        else:
            with open(logfile, 'w') as f:
                call(command, cwd=cwd, stdout=f)

    except Exception as E:
        # TODO: Change to warnings.warn()
        print('Could not successfully run test! Caught error: {0}'.format(E))
    else:
        print('Finished running test!')
