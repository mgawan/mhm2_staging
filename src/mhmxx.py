#!/usr/bin/env python

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
#import re

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

def get_hwd_cores_per_node():
    """Query the hardware for physical cores"""
    try:
        import psutil
        cores = psutil.cpu_count(logical=False)
    except (NameError, ImportError):
        import platform
        import multiprocessing
        cpus = multiprocessing.cpu_count()
        hyperthreads = 1
        if platform.system() == 'Darwin':
            for line in os.popen('sysctl -n hw.physicalcpu').readlines():
                 hyperthreads = cpus / int(line)
        else:
            for line in os.popen('lscpu').readlines():
                if line.startswith('Thread(s) per core'):
                    hyperthreads = int(line.split()[3])
        cores = cpus / hyperthreads
    return cores


def get_job_id():
    """Query the environment for a job"""
    for key in ['PBS_JOBID', 'SLURM_JOBID', 'LSB_JOBID', 'JOB_ID', 'LOAD_STEP_ID']:
      if key in os.environ:
        return os.environ.get(key)
    return None

def get_job_name():
    """Query the env for the name of a job"""
    for key in ['PBS_JOBNAME', 'JOBNAME', 'SLURM_JOB_NAME', 'LSB_JOBNAME', 'JOB_NAME']:
      if key in os.environ:
        return os.environ.get(key)
    return None

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


def get_job_cores_per_node(defaultCores = get_hwd_cores_per_node()):
    """Query the job environment for the number of cores per node to use, if available"""
    # Only trust this environment variable from slurm, otherwise trust the hardware                                                                                                   
    ntasks_per_node = os.environ.get('SLURM_NTASKS_PER_NODE')
    if ntasks_per_node:
        return int(ntasks_per_node)
    # This SLURM variable defaults to all the hyperthreads if not overriden by the sbatch option --ntasks-per-node                                                                    
    ntasks_per_node = os.environ.get('SLURM_TASKS_PER_NODE') # SLURM_TASKS_PER_NODE=32(x4)                                                                                            
    if ntasks_per_node:
        if ntasks_per_node.find('(') > 0:
            ntasks_per_node = int(ntasks_per_node[:ntasks_per_node.find('(')])
        else:
            ntasks_per_node = int(ntasks_per_node)
        if ntasks_per_node <= defaultCores:
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


def get_job_nodes():
    """Query the job environment for the number of nodes"""
    if is_slurm_job():
        return get_slurm_job_nodes()
    if is_pbs_job():
        return get_pbs_job_nodes()
    if is_ge_job():
        return get_ge_job_nodes()
    print("Warning: could not determine the number of nodes in this unsupported scheduler - using 1")
    return 1


def which(file_name):
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
    print_red('\nFATAL ERROR:', *args)
    sys.stdout.flush()
    sys.stderr.flush()
    exit_all(1)


def check_exec(cmd, args, expected):
    test_exec = which(cmd)
    if not test_exec:
        die('Cannot find', cmd)
    try:
        result = subprocess.check_output([test_exec, args]).decode()
        if expected not in result:
            die(test_exec, 'failed to execute')
    except subprocess.CalledProcessError as err:
        die('Could not execute', test_exec +':', err)

def capture_err(err_msgs):
    global _proc
    global _stop_thread
    for line in iter(_proc.stderr.readline, b''):
        line = line.decode()
        if 'WARNING' in line:
            sys.stderr.write(line)
            sys.stderr.flush()
        err_msgs.append(line)
        if _stop_thread:
            return
    _proc.wait()

def print_err_msgs(err_msgs):
    global _output_dir
    err_msgs.append('==============================================')
    print_red("Check " + _output_dir + "err.log for details")
    # keep track of all msg copies so we don't print duplicates
    seen_msgs = {}
    with open(_output_dir + 'err.log', 'w') as f:
        for msg in err_msgs:
            clean_msg = msg.strip()
            #clean_msg = re.sub('\(proc .+\)', '(proc XX)', msg.strip())
            if clean_msg not in seen_msgs:
                f.write(clean_msg + '\n')
                f.flush()
                seen_msgs[clean_msg] = True

    
def main():
    global _orig_sighdlr
    global _proc
    global _output_dir
    global _err_thread
    
    _orig_sighdlr = signal.getsignal(signal.SIGINT)
    signal.signal(signal.SIGINT, handle_interrupt)

    argparser = argparse.ArgumentParser(add_help=False)
    argparser.add_argument("--auto-resume", action="store_true", help="Automatically resume after a failure")
    argparser.add_argument("--shared-heap", default="10%", help="Shared heap as a percentage of memory")
    
    options, unknown_options = argparser.parse_known_args()

    if options.auto_resume:
        print("--auto-resume is enabled: will try to restart if run fails")
    
    check_exec('upcxx-run', '-h', 'UPC++')
    # expect mhmxx to be in same directory as mhmxx.py
    mhmxx_binary_path = os.path.split(sys.argv[0])[0] + '/mhmxx'
    if not which(mhmxx_binary_path):
        die("Cannot find binary mhmxx")
    cores_per_node = get_hwd_cores_per_node()
    num_nodes = get_job_nodes()
    cmd = ['upcxx-run', '-n', str(cores_per_node * num_nodes), '-N', str(num_nodes), '-shared-heap', options.shared_heap, '--', 
           mhmxx_binary_path];
    cmd.extend(unknown_options)
    print('Executing:')
    print(' '.join(cmd))
    
    err_msgs = []
    while True:
      completed_round = False
      try:
          _proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
          # thread captures the error stream
          _err_thread = threading.Thread(target=capture_err, args=(err_msgs,))
          _err_thread.start()
          for line in iter(_proc.stdout.readline, b''):
              line = line.decode()
              sys.stdout.write(line)
              sys.stdout.flush()
              if line.startswith('  output = '):
                  _output_dir = line.split()[2]
                  if _output_dir[-1] != '/':
                      _output_dir += '/'
              if 'Completed ' in line and 'initialization' not in line:
                  completed_round = True
                  
          _err_thread.join()
          if _proc.returncode not in [0, -15] or not status:
              signame = ''
              if -_proc.returncode <= len(SIGNAMES) and _proc.returncode < 0:
                  signame = ' (' + SIGNAMES[-_proc.returncode - 1] + ')'
              print_red("\nERROR: subprocess terminated with return code ", -_proc.returncode, signame);
              err_msgs.append("ERROR: subprocess terminated with return code " + str(-_proc.returncode) + signame);
              print_err_msgs(err_msgs)
              if completed_round and options.auto_resume:
                  print_red('Trying to restart...')
                  cmd.append('--restart')
              else:
                  if options.auto_resume:
                      print_red("No additional completed round. Could not restart, exiting...")
                  return signal.SIGABRT
          else:
              for msg in err_msgs:
                  sys.stderr.write(msg)
                  sys.stderr.flush()
              break
      except:
          print_red("\nSubprocess failed to start")
          traceback.print_tb(sys.exc_info()[2], limit=100)
          print_err_msgs(err_msgs)
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
    # remove the .py from this script as the mhmxx wrapper needs to be excecuted for proper environment variable detection
    sys.argv[0] = os.path.splitext(sys.argv[0])[0] 
    status = 1
    try:
        status = main()
    except SystemExit:
        raise
    except:
        e = sys.exc_info()[0]
        print_red("\n", "\nCaught an exception %s in mhmxx.py!\n\n" % e); 
        traceback.print_exc(file=sys.stderr)
    finally:
        exit_all(status)

