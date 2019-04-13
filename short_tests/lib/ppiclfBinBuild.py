import os
import sys
from subprocess import call, check_call, Popen, PIPE, STDOUT

def build_ppiclf(source_root, cwd=None, opts=None, verbose=False):

    if not opts:
        _opts = {}
    else:
        _opts = opts.copy()
    _opts.update(SOURCE_ROOT=source_root)

    print('Compiling ppiclf...')
    print('    Using working directory "{0}"'.format(cwd))
    for key, val in _opts.iteritems():
        print('    Using {0}="{1}"'.format(key, val))

    my_env = os.environ.copy()
    if source_root         : my_env["SOURCE_ROOT"] = source_root
    if _opts.get('F77')    : my_env["FC"] = _opts.get('F77') 
    if _opts.get('CC')     : my_env["CC"] = _opts.get('CC')
    if _opts.get('PPLIST') : my_env["PPLIST"] = _opts.get('PPLIST') 

    logfile     = os.path.join(cwd, 'build.log')

    # copy F-File
    proc = Popen(
        ['cp',os.path.join(format(cwd),'user_routines/ppiclf_user.f'),os.path.join(source_root, 'source', 'ppiclf_user.f')], 
        cwd=cwd,
        env=my_env,
        stdin=PIPE)
    proc.wait()

    # copy H-File
    proc = Popen(
        ['cp',os.path.join(format(cwd),'user_routines/PPICLF_USER.h'),os.path.join(source_root, 'source', 'PPICLF_USER.h')], 
        cwd=cwd,
        env=my_env,
        stdin=PIPE)
    proc.wait()

    # Clean ppiclF library
    proc = Popen(
        ['make','clean'], 
        cwd=source_root,
        env=my_env,
        stdin=PIPE)
    proc.wait()

    # Make ppiclF library
    proc = Popen(
        'make', 
        cwd=source_root,
        env=my_env,
        stdin=PIPE, 
        stderr=STDOUT) 
    proc.wait()

    # Clean example case
    proc = Popen(
        ['make','clean'], 
        cwd=cwd,
        env=my_env,
        stdin=PIPE, 
        stderr=STDOUT) 
    proc.wait()

    # Make example case
    proc = Popen(
        'make', 
        cwd=cwd,
        env=my_env,
        stdin=PIPE, 
        stderr=STDOUT) 
    proc.wait()

    if proc.returncode != 0:
       f = open(logfile, "r")
       text = f.read()
       print text
       f.close()
       exit(-1)
