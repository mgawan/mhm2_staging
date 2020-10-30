#!/usr/bin/env python

 # HipMer v 2.0, Copyright (c) 2020, The Regents of the University of California,
 # through Lawrence Berkeley National Laboratory (subject to receipt of any required
 # approvals from the U.S. Dept. of Energy).  All rights reserved."

 # Redistribution and use in source and binary forms, with or without modification,
 # are permitted provided that the following conditions are met:

 # (1) Redistributions of source code must retain the above copyright notice, this
 # list of conditions and the following disclaimer.

 # (2) Redistributions in binary form must reproduce the above copyright notice,
 # this list of conditions and the following disclaimer in the documentation and/or
 # other materials provided with the distribution.

 # (3) Neither the name of the University of California, Lawrence Berkeley National
 # Laboratory, U.S. Dept. of Energy nor the names of its contributors may be used to
 # endorse or promote products derived from this software without specific prior
 # written permission.

 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 # EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 # OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 # SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 # INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 # TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 # BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 # CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 # ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 # DAMAGE.

 # You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades
 # to the features, functionality or performance of the source code ("Enhancements") to
 # anyone; however, if you choose to make your Enhancements available either publicly,
 # or directly to Lawrence Berkeley National Laboratory, without imposing a separate
 # written license agreement for such Enhancements, then you hereby grant the following
 # license: a  non-exclusive, royalty-free perpetual license to install, use, modify,
 # prepare derivative works, incorporate into other computer software, distribute, and
 # sublicense such enhancements or derivative works thereof, in binary and source code
 # form.

from __future__ import print_function

import signal
import subprocess
import sys
import os
import datetime
import time
import traceback
import argparse
import threading
import io
import string
import multiprocessing
import collections

SIGNAMES = ['SIGHUP', 'SIGINT', 'SIGQUIT', 'SIGILL', 'SIGTRAP', 'SIGABRT', 'SIGBUS', 'SIGFPE', 'SIGKILL', 'SIGUSR1',
            'SIGSEGV', 'SIGUSR2', 'SIGPIPE', 'SIGALRM', 'SIGTERM', 'SIGSTKFLT', 'SIGCHLD', 'SIGCONT', 'SIGSTOP', 'SIGTSTP',
            'SIGTTIN', 'SIGTTOU', 'SIGURG', 'SIGXCPU', 'SIGXFSZ', 'SIGVTALRM', 'SIGPROF', 'SIGWINCH', 'SIGIO', 'SIGPWR', 'SIGSYS']

_orig_sighdlr = None
_proc = None
_output_dir = ''
_err_thread = None
_stop_thread = False

def print_red(*args):
    print("\033[91m", *args, sep='',  end='', file=sys.stderr)
    print("\033[00m", file=sys.stderr)

_defaultCores = None
def get_hdw_cores_per_node():
    """Query the hardware for physical cores"""
    global _defaultCores
    if _defaultCores is not None:
        return _defaultCores
    try:
        import psutil
        cores = psutil.cpu_count(logical=False)
        print("Found %d cpus from psutil" % cores)
    except (NameError, ImportError):
        #print("Could not get cpus from psutil")
        pass
    # always trust lscpu, not psutil
    # NOTE some versions of psutil has bugs and comes up with the *WRONG* physical cores
    if True:
        import platform
        cpus = multiprocessing.cpu_count()
        hyperthreads = 1
        if platform.system() == 'Darwin':
            for line in os.popen('sysctl -n hw.physicalcpu').readlines():
                hyperthreads = cpus / int(line)
            print("Found %d cpus and %d hyperthreads from sysctl" % (cpus, hyperthreads))
        else:
            for line in os.popen('lscpu').readlines():
                if line.startswith('Thread(s) per core'):
                    hyperthreads = int(line.split()[3])
            print("Found %d cpus and %d hyperthreads from lscpu" % (cpus, hyperthreads))
        cores = int(cpus / hyperthreads)
    _defaultCores = cores
    return cores

def get_job_id():
    """Query the environment for a job"""
    for key in ['PBS_JOBID', 'SLURM_JOBID', 'LSB_JOBID', 'JOB_ID', 'COBALT_JOBID', 'LOAD_STEP_ID', 'LBS_JOBID']:
        if key in os.environ:
            return os.environ.get(key)
    return str(os.getpid())

def get_job_name():
    """Query the env for the name of a job"""
    for key in ['PBS_JOBNAME', 'JOBNAME', 'SLURM_JOB_NAME', 'LSB_JOBNAME', 'JOB_NAME', 'LSB_JOBNAME']:
        if key in os.environ:
            return os.environ.get(key)
    return ""

def is_cobalt_job():
    return os.environ.get('COBALT_JOBID') is not None

def is_slurm_job():
    return os.environ.get('SLURM_JOB_ID') is not None

def is_pbs_job():
    return os.environ.get('PBS_JOBID') is not None

def is_lsb_job():
    return os.environ.get('LSB_JOBID') is not None

def is_ge_job():
    return os.environ.get('JOB_ID') is not None

def is_ll_job():
    return os.environ.get('LOAD_STEP_ID') is not None


def get_slurm_cores_per_node(defaultCores = 0):
    # Only trust this environment variable from slurm, otherwise trust the hardware
    if defaultCores == 0:
        defaultCores = get_hdw_cores_per_node()
    ntasks_per_node = os.environ.get('SLURM_NTASKS_PER_NODE')
    if ntasks_per_node:
        print("Found tasks per node from SLURM_NTASKS_PER_NODE=", ntasks_per_node)
        return int(ntasks_per_node)
    # This SLURM variable defaults to all the hyperthreads if not overriden by the sbatch option --ntasks-per-node
    ntasks_per_node = os.environ.get('SLURM_TASKS_PER_NODE') # SLURM_TASKS_PER_NODE=32(x4)
    if ntasks_per_node:
        if ntasks_per_node.find('(') > 0:
            ntasks_per_node = int(ntasks_per_node[:ntasks_per_node.find('(')])
        else:
            ntasks_per_node = int(ntasks_per_node)
        if ntasks_per_node <= defaultCores:
            print("Detected slurm job restricts cores to ", ntasks_per_node, " because of SLURM_TASKS_PER_NODE=", os.environ.get('SLURM_TASKS_PER_NODE'))
            return ntasks_per_node
        print("Using default cores of ", defaultCores, ". Ignoring tasks per node ", ntasks_per_node, " from SLURM_TASKS_PER_NODE=", os.environ.get('SLURM_TASKS_PER_NODE'))
    return defaultCores

def get_cobalt_cores_per_node():
    ntasks_per_node = os.environ.get('COBALT_PARTCORES')
    return int(ntasks_per_node)

def get_lsb_cores_per_node():
    # LSB_MCPU_HOSTS=batch2 1 h22n07 42 h22n08 42
    lsb_mcpu = os.environ.get('LSB_MCPU_HOSTS')
    host_core = lsb_mcpu.split()
    return int(host_core[-1])


def get_job_cores_per_node(defaultCores = 0):
    """Query the job environment for the number of cores per node to use, if available"""
    if defaultCores == 0:
        defaultCores = get_hdw_cores_per_node()
    if 'GASNET_PSHM_NODES' in os.environ:
        print("Detected procs_per_node from GASNET_PSHM_NODES=",os.getenv('GASNET_PSHM_NODES'))
        return int(os.getenv('GASNET_PSHM_NODES'))
    ntasks_per_node = None
    if is_slurm_job():
        ntasks_per_node = get_slurm_cores_per_node(defaultCores)
    if is_lsb_job():
        ntasks_per_node = get_lsb_cores_per_node()
    if is_cobalt_job():
        ntasks_per_node = get_cobalt_cores_per_node()
    if ntasks_per_node is not None:
        return ntasks_per_node
    return defaultCores

def get_slurm_job_nodes():
    """Query the SLURM job environment for the number of nodes"""
    nodes = os.environ.get('SLURM_JOB_NUM_NODES')
    if nodes is None:
        nodes = os.environ.get('SLURM_NNODES')
    if nodes:
        return int(nodes)
    print("Warning: could not determine the number of nodes in this SLURM job (%d). Only using 1" % (get_job_id()))
    return 1

def get_lsb_job_nodes():
    """Query the LFS job environment for the number of nodes"""
    # LSB_MCPU_HOSTS=batch2 1 h22n07 42 h22n08 42
    nodes = os.environ.get('LSB_MCPU_HOSTS')
    if nodes:
        return int( (len(nodes.split()) - 2) / 2)
    print("Warning: could not determine the number of nodes in this LSF job (%s). Only using 1" % (get_job_id()))
    return 1

def get_pbs_job_nodes():
    """Query the PBS job environment for the number of nodes"""
    nodesfile = os.environ.get('PBS_NODEFILE')
    if nodesfile is not None:
        nodes = 0
        with open(nodesfile, 'r') as f:
            for line in f:
                nodes += 1
        return nodes
    print("Warning: could not determine the number of nodes in this PBS job (%d). Only using 1" % (get_job_id()))
    return 1

def get_ge_job_nodes():
    """Query the Grid Engine job environment for the number of nodes"""
    nodes = os.environ.get("NHOSTS")
    if nodes is not None:
        return int(nodes)
    print("Warning: could not determine the number of nodes in this SGE job (%d). Only using 1" % (get_job_id()))
    return 1

def get_cobalt_job_nodes():
    """Query the COBALT job environment for the number of nodes"""
    nodes = os.environ.get("COBALT_JOBSIZE")
    if nodes is not None:
        return int(nodes)
    print("Warning: could not determine the number of nodes in this COBALT job (%s). Only using 1" % (get_job_id()))
    return 1

def get_job_nodes():
    """Query the job environment for the number of nodes"""
    if is_slurm_job():
        return get_slurm_job_nodes()
    if is_lsb_job():
        return get_lsb_job_nodes()
    if is_pbs_job():
        return get_pbs_job_nodes()
    if is_ge_job():
        return get_ge_job_nodes()
    if is_cobalt_job():
        return get_cobalt_job_nodes()
    print("Warning: could not determine the number of nodes in this unsupported scheduler job (%s). Only using 1" % (get_job_id()))
    return 1

def get_job_desc():
    job_name = get_job_name()
    if job_name == "":
        return "PID " + get_job_id()
    return "job " + get_job_id() + " (" + job_name + ")"


def which(file_name):
    if os.path.exists(file_name) and os.access(file_name, os.X_OK):
        return file_name
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None

def handle_interrupt(signum, frame):
    global _orig_sighdlr
    global _stop_thread
    print_red('\n\nInterrupt received, signal', signum)
    _stop_thread = True
    signal.signal(signal.SIGINT, _orig_sighdlr)
    exit_all(1)


def exit_all(status):
    global _proc
    if _proc:
        #os.kill(_proc.pid, signal.SIGINT)
        try:
            _proc.terminate()
        except OSError:
            pass
            #print("Process ", _proc, " is already terminated\n")
    sys.exit(status)


def die(*args):
    print_red('\nFATAL ERROR: ', *args)
    sys.stdout.flush()
    sys.stderr.flush()
    exit_all(1)


def check_exec(cmd, args, expected):
    test_exec = which(cmd)
    if not test_exec:
        die('Cannot find ', cmd)
    try:
        result = subprocess.check_output([test_exec, args]).decode()
        if expected not in result:
            die(test_exec, ' failed to execute')
    except subprocess.CalledProcessError as err:
        die('Could not execute ', test_exec +': ', err)

def capture_err(err_msgs):
    global _proc
    global _stop_thread
    for line in iter(_proc.stderr.readline, b''):
        line = line.decode()
        # filter out all but warnings
        # errors causing crashes will come to light later
        if 'WARNING' in line:
            if 'GASNet was configured without multi-rail support' not in line and 'GASNET_AM_CREDITS_SLACK reduced to GASNET_AM_CREDITS_PP' not in line:
                sys.stderr.write(line)
                sys.stderr.flush()
        # FIXME: check for messages about memory failures
        if 'UPC++ could not allocate' in line:
            print_red('ERROR: UPC++ memory allocation failure')
        err_msgs.append(line)
        if _stop_thread:
            return
    _proc.wait()

def check_err_msgs(err_msgs):
    warnings = []
    errors = []
    for msg in err_msgs:
        msg = msg.strip()
        if 'WARNING' in msg:
            warnings.append(msg)
        elif msg[:2] == '+ ':
            # this is 'set -x' console echo of a command
            pass
        elif 'GASNet reporting enabled' in msg:
            # this is just info
            pass
        elif msg != '':
            errors.append(msg)
    if len(warnings) > 0:
        print('There were', len(warnings), 'warnings:', file=sys.stderr)
        for warning in list(collections.OrderedDict.fromkeys(warnings)):
            sys.stderr.write(warning + '\n')
    if len(errors) > 0:
        print('There were', len(errors), 'errors:', file=sys.stderr)
        print_sigint = False
        for err in errors:
            if 'SIGINT(2)' in err and print_sigint:
                continue
            print_sigint = True
            sys.stderr.write(err + '\n')
    return len(errors) + len(warnings)

def print_err_msgs(err_msgs, return_status):
    global _output_dir
    num_problems = check_err_msgs(err_msgs)
    err_msgs.append('==============================================')
    if len(_output_dir) == 0:
        _output_dir = os.getcwd() + "/"
        # we have not yet entered the output directory, so this is a failure of the command line
        # and we need to dump all the error messages to the console
        print_red("No output dir was created yet")
        for msg in err_msgs:
            print(msg)
            sys.exit(return_status)
    else:
        if _output_dir[0] != '/':
            _output_dir = os.getcwd() + "/" + _output_dir
        suspect_oom = None
        if return_status != 0:
            if return_status == 9: # SIGKILL
                suspect_oom = "Got SIGKILLed"
            err_msgs.append("Return status: %d\n" % (return_status))
            print_red("MHM2 failed")
        # keep track of all msg copies so we don't print duplicates
        seen_msgs = {}
        per_rank_dir = _output_dir + 'per_rank/'
        if not os.path.exists(per_rank_dir):
            per_rank_dir = _output_dir
        err_log = per_rank_dir + 'err.log'
        with open(err_log, 'a') as f:
            for msg in err_msgs:
                clean_msg = msg.strip()
                #clean_msg = re.sub('\(proc .+\)', '(proc XX)', msg.strip())
                if clean_msg not in seen_msgs:
                    f.write(clean_msg + '\n')
                    f.flush()
                    seen_msgs[clean_msg] = True
                    if 'SIGBUS' in clean_msg or 'bound CqGetEvent GNI_RC_TRANSACTION_ERROR' in clean_msg or 'oom-kill' in clean_msg or 'bad_alloc' in clean_msg or 'SIGKILL' in clean_msg \
                        or 'Cannot allocate memory' in clean_msg or 'mmap failed' in clean_msg:
                        suspect_oom = clean_msg
            if suspect_oom is not None:
                f.write("Out of memory is suspected because of: %s\n" %(suspect_oom))
                print_red("Out of memory is suspected based on the errors in err.log such as: ", suspect_oom, "\n")
        if num_problems > 0:
            print_red("Check " + err_log + " for details")

def main():
    global _orig_sighdlr
    global _proc
    global _output_dir
    global _err_thread

    start_time = time.time()
    _orig_sighdlr = signal.getsignal(signal.SIGINT)
    signal.signal(signal.SIGINT, handle_interrupt)

    argparser = argparse.ArgumentParser(add_help=False)
    argparser.add_argument("--auto-resume", action="store_true", help="Automatically resume after a failure")
    argparser.add_argument("--shared-heap", default="10%", help="Shared heap as a percentage of memory")
    #argparser.add_argument("--procs-per-node", default=0, help="Processes to spawn per node (default auto-detect cores)")
    argparser.add_argument("--procs", default=0, type=int, help="Total numer of processes")
    argparser.add_argument("--trace-dir", default=None, help="Output directory for stacktrace")
    argparser.add_argument("--stats-dir", default=None, help="Output directory for stacktrace")

    options, unknown_options = argparser.parse_known_args()

    if options.auto_resume:
        print("--auto-resume is enabled: will try to restart if run fails")

    check_exec('upcxx-run', '-h', 'UPC++')
    # expect mhm2 to be in same directory as mhm2.py
    mhm2_binary_path = os.path.split(sys.argv[0])[0] + '/mhm2'
    if not (os.path.exists(mhm2_binary_path) or which(mhm2_binary_path)):
        die("Cannot find binary mhm2 in '", mhm2_binary_path, "'")

    #cores_per_node = int(options.procs_per_node)
    #if cores_per_node == 0:
    #    cores_per_node = get_job_cores_per_node()
    num_nodes = get_job_nodes()
    if options.procs == 0:
        options.procs = num_nodes * get_job_cores_per_node()

    cmd = ['upcxx-run', '-n', str(options.procs), '-N', str(num_nodes)]

    # special spawner for summit -- executes jsrun and picks up job size from the environment!
    if 'LMOD_SYSTEM_NAME' in os.environ and os.environ['LMOD_SYSTEM_NAME'] == "summit":
        print("This is Summit - executing custom script mhm2-upcxx-run-summit to spawn the job")
        # expect mhm2-upcxx-run-summit to be in same directory as mhm2.py too
        cmd = [mhm2_binary_path + "-upcxx-run-summit"]
        if 'UPCXX_RUN_SUMMIT_OPTS' in os.environ:
            cmd.extend(os.environ['UPCXX_RUN_SUMMIT_OPTS'].split())

    if 'UPCXX_SHARED_HEAP_SIZE' not in os.environ:
        cmd.extend(['-shared-heap', options.shared_heap]) # both upcxx-run and upcxx-run-summit support this

    cmd.extend(['--', mhm2_binary_path])
    cmd.extend(unknown_options)

    print("Executing mhm2 with " + get_job_desc() + " on " + str(num_nodes) + " nodes.")
    print("Executing as: " + " ".join(sys.argv))

    cores = get_job_cores_per_node()
    noderanks = '0'
    halfnoderanks = '0,%d' % (cores/2)
    for n in range(1, num_nodes):
        noderanks += ',' + str(n*cores)
        halfnoderanks += ',' + str(n*cores) + ',' + str(n*cores+cores/2)

    # set extra GASNET environments from build and/or options to mhm2.py
    runtime_vars = """@MHM2PY_RUNTIME_ENV@"""
    if runtime_vars == '@MHM2PY_RUNTIME' + '_ENV@':
        runtime_vars = ''
    runtime_output_vars = ''

    if options.stats_dir is not None:
        if not os.path.isdir(options.stats_dir):
            os.mkdir(options.stats_dir)
        runtime_vars += ' GASNET_STATSFILE="%s/stats.%%", ' % (os.path.realpath(options.stats_dir))
        runtime_vars += runtime_output_vars
        if os.environ.get("GASNET_STATSNODES") is None:
            runtime_vars += ' GASNET_STATSNODES="%s", ' % noderanks
    if options.trace_dir is not None:
        if not os.path.isdir(options.trace_dir):
            os.mkdir(options.trace_dir)
        runtime_vars += ' GASNET_TRACEFILE="%s/trace.%%", ' % (os.path.realpath(options.trace_dir))
        if os.environ.get("GASNET_TRACENODES") is None:
            runtime_vars += ' GASNET_TRACENODES="%s", ' % halfnoderanks
        if os.environ.get("GASNET_TRACEMASK") is None:
            runtime_vars += ' GASNET_TRACEMASK="GPWBNIH", ' # some of the more useful and less verbose trace options

    # it appears that this GASNET_COLL_SCRATCH_SIZE  is still needed
    print("Setting GASNET_COLL_SCRATCH_SIZE=4M", runtime_vars)
    runenv = eval('dict(os.environ, GASNET_COLL_SCRATCH_SIZE="4M", %s MHM2_RUNTIME_PLACEHOLDER="")' % (runtime_vars))
    #print("Runtime environment: ", runenv)

    mhm2_lib_path = os.path.split(sys.argv[0])[0] + '/../lib'
    if not os.path.exists(mhm2_lib_path):
        die("Cannot find mhm2 lib install in '", mhm2_lib_path, "'")

    # This should no longer be necessary with the static GPU build fixes of df8cc23, leaving here in case problems reoccur
    #if which('nvcc'):
    #    # FIXME: this ugly hack is because we need to load a shared library on Cori GPU nodes,
    #    # which can't be done with the craype environment. Not needed anywhere else :(
    #    # The intel library path is only needed for the intel compiler. Sigh.
    #    runenv['LD_LIBRARY_PATH'] = mhm2_lib_path + ':/usr/lib64/slurmpmi/:/opt/intel/compilers_and_libraries_2019.3.199/linux/compiler/lib/intel64_lin/'
    #    print('Setting LD_LIBRARY_PATH=' + runenv['LD_LIBRARY_PATH'])

    restarting = False
    err_msgs = []
    while True:
        print(str(datetime.datetime.now()) + ' ' + 'executing:\n', ' '.join(cmd))
        started_exec = False
        completed_round = False
        try:
            _proc = subprocess.Popen(cmd, env=runenv, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # thread captures the error stream
            _err_thread = threading.Thread(target=capture_err, args=(err_msgs,))
            _err_thread.start()
            for line in iter(_proc.stdout.readline, b''):
                if not started_exec:
                    print('Started executing at ' + str(datetime.datetime.now()) + ' with PID ' + str(_proc.pid))
                    started_exec = True

                line = line.decode()
                sys.stdout.write(line)
                sys.stdout.flush()
                if '  output = ' in line:
                    _output_dir = line.split()[3]
                    onlyascii = ''.join([s for s in _output_dir if ord(s) < 127 and ord(s) >= 32])
                    _output_dir = onlyascii
                    if _output_dir.endswith('[0m'):
                        _output_dir = _output_dir[:-3]
                    if _output_dir[-1] != '/':
                        _output_dir += '/'
                    # get rid of any leftover error logs if not restarting
                    try:
                        if not restarting:
                            os.rename(_output_dir + '/per_rank/err.log', _output_dir + '/per_rank/err.log' + str(datetime.datetime.now().isoformat()))
                    except:
                        pass

                if 'Completed ' in line and 'initialization' not in line:
                    completed_round = True

            _err_thread.join()
            if _proc.returncode < 0:
                _proc.returncode *= -1
            if _proc.returncode not in [0, 15] or not status:
                signame = ''
                if _proc.returncode <= len(SIGNAMES) and _proc.returncode > 0:
                    signame = ' (' + SIGNAMES[_proc.returncode - 1] + ')'
                if _proc.returncode == 127:
                    # 127 is the return code from the CLI parser, so we don't want to print this
                    found_error = False
                    for msg in err_msgs:
                        if msg.startswith('Error'):
                            print(msg, end='')
                            found_error = True
                        elif msg.startswith('INI was not able to parse'):
                            print("  Unable to parse entry '" + msg.split()[-1] + "' in the config file")
                        else:
                            print(" ", msg, end='')
                    if found_error:
                        return 127
                    else:
                        return 0
                print_red("\nERROR: subprocess terminated with return code ", _proc.returncode)
                signals_found = {}
                for err_msg in err_msgs:
                    for signame in SIGNAMES:
                        if signame in err_msg:
                            if not signame in signals_found:
                                signals_found[signame] = 0
                            signals_found[signame] += 1
                for signame in SIGNAMES:
                    if signame in signals_found:
                        print_red("  Found ", signals_found[signame], " occurences of ", signame)
                        got_signal = signame
                #err_msgs.append("ERROR: subprocess terminated with return code " + str(_proc.returncode) + " " + signame)
                print_err_msgs(err_msgs, _proc.returncode)
                if completed_round and options.auto_resume:
                    print_red('Trying to restart with output directory ', _output_dir)
                    restarting = True
                    err_msgs = []
                    cmd.append('--restart')
                    if _output_dir[:-1] not in cmd:
                        cmd.extend(['-o', _output_dir])
                    time.sleep(5)
                else:
                    if options.auto_resume:
                        print_red("No additional completed round. Could not restart, exiting...")
                    return signal.SIGABRT
            else:
                final_assembly = _output_dir + "final_assembly.fasta"
                if os.path.exists(final_assembly):
                    print("Final assembly can be found at ", final_assembly)
                else:
                    err_msgs.append("Could not find the final assembly!  It should be at %s\n" % (final_assembly))
                print_err_msgs(err_msgs, _proc.returncode)
                print('Overall time taken (including any restarts): %.2f s' % (time.time() - start_time))
                break
        except:
            print_red("Got an exception")
            traceback.print_tb(sys.exc_info()[2], limit=100)
            print_err_msgs(err_msgs, -1)
            if _proc:
                try:
                    print_red("\nTerminating subprocess after exception: ", sys.exc_info(), "\n")
                    traceback.print_tb(sys.exc_info()[2], limit=100)
                    _proc.terminate()
                except OSError:
                    pass
                except:
                    print_red("\nUnexpected error in forced termination of subprocess: ", sys.exc_info())
                    traceback.print_tb(sys.exc_info()[2], limit=100)
                    raise
            raise

    return 0

if __name__ == "__main__":
    status = 1
    try:
        status = main()
    except SystemExit:
        if status != 127:
            raise
    except:
        e = sys.exc_info()[0]
        print_red("\n", "\nCaught an exception %s in mhm2.py!\n\n" % e)
        traceback.print_exc(file=sys.stderr)
    finally:
        exit_all(status)
