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


orig_sighdlr = None
proc = None


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
    print('\n\nInterrupt received, signal', signum)
    signal.signal(signal.SIGINT, orig_sighdlr)
    exit_all(1)


def exit_all(status):
    if proc:
        #os.kill(proc.pid, signal.SIGINT)
        try:
            proc.terminate()
        except OSError:
            pass
            #print("Process ", proc, " is already terminated\n")
    sys.exit(status)


def die(*args):
    print('\nFATAL ERROR:', *args, file=sys.stderr)
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


def handle_failure_termination(proc, verbose=False):
    print("mhmxx detected a failure in process pid=%d, reading output and terminating at %s" % \
          (proc.pid, str(datetime.datetime.now())))
    # get any remaininginput, then send two interrupts and a term
    if proc:
        try:
            try:
                if proc.poll() is None:
                    time.sleep(1)
                    if proc.poll() is None:
                        print("Sending SIGINT at %s\n" % (str(datetime.datetime.now())))
                        proc.send_signal(2)
            except OSError:
                print("Process is already gone\n")
                pass
                  
            if proc.poll() is None:
                time.sleep(1)
                if proc.poll() is None:
                    print("Sending SIGINT again at %s\n" % (str(datetime.datetime.now())))
                    proc.send_signal(2)
                    time.sleep(1)
                    if proc.poll() is None:
                        time.sleep(1)
                        if proc.poll() is None:
                            print("Sending SIGTERM at %s\n" % (str(datetime.datetime.now())))
                            proc.terminate()
                            time.sleep(1)
                            if proc.poll() is None:
                                print("Sending SIGKILL at %s\n" % (str(datetime.datetime.now())))
                                proc.kill()
                                time.sleep(1)

            outs,errs = proc.communicate(None)
            if outs:
                outs = outs.decode()
                print("Got some stdout from the failed process at %s" % (str(datetime.datetime.now())))
                if verbose:
                    print(outs)
                    sys.stdout.flush()
                print(outs)
                for line in outs:
                    if 'srun: error' in line and 'REQUEST_FILE_BCAST' in line and os.environ.get('SLURM_BCAST') is not None: # Issue228
                        print("Detected a SLURM_BCAST error, disabling that optimization.")
                        del os.environ['SLURM_BCAST']
                        os.environ['NO_SLURM_BCAST'] = '1'
            if errs:
                errs = errs.decode()
                print("Got some stderr from the failed process at %s" % (str(datetime.datetime.now())))
                print(errs)
                sys.stderr.write(errs)
                sys.stderr.flush()
                sys.stdout.flush()
                for line in errs:
                    # Issue228
                    if 'srun: error' in line and 'REQUEST_FILE_BCAST' in line and os.environ.get('SLURM_BCAST') is not None: 
                        print("Detected a SLURM_BCAST error, disabling that optimization.")
                        del os.environ['SLURM_BCAST']
                        os.environ['NO_SLURM_BCAST'] = '1'
        except OSError:
            print("Process is gone\n")
            pass
        except:
            print("Unexpected error: ", sys.exc_info())
            traceback.print_tb(sys.exc_info()[2], limit=100)
            pass
      
      
def main():
    global orig_sighdlr
    global proc
    orig_sighdlr = signal.getsignal(signal.SIGINT)
    signal.signal(signal.SIGINT, handle_interrupt)

    argparser = argparse.ArgumentParser(add_help=False)
    argparser.add_argument("--auto-resume", action="store_true", help="Automatically resume after a failure")
    argparser.add_argument("--shared-heap", help="Shared heap as a percentage of memory")
    
    options, unknown_options = argparser.parse_known_args()

    if options.auto_resume:
        print("auto resume is enabled: will try to restart if run fails")
    
    check_exec('upcxx-run', '-h', 'UPC++')
    status = True
    # expect mhmxx to be in same directory as mhmxx.py
    mhmxx_binary_path = os.path.split(sys.argv[0])[0] + '/mhmxx'
    if not which(mhmxx_binary_path):
        die("Cannot find binary mhmxx")
    cores_per_node = get_hwd_cores_per_node()
    num_nodes = get_job_nodes()
    cmd = ['upcxx-run', '-n', str(cores_per_node * num_nodes), '-N', str(num_nodes), '-shared-heap', options.shared_heap, '--', 
           mhmxx_binary_path];
#    cmd.extend(sys.argv[1:])
    cmd.extend(unknown_options)
    print('Executing:')
    print(' '.join(cmd))
    
    while True:
      completed = False
      try:
          proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
          for line in iter(proc.stdout.readline, b''):
              line = line.decode()
              print(line, end='')
              sys.stdout.flush()
              # did we complete any new rounds?
              if "Completed " in line:
                  completed = True
          if status:
              proc.wait()
          else:
              handle_failure_termination(proc, run_log_file, options.verbose)
              proc.wait()

          if proc.returncode not in [0, -15] or not status:
              #print("ERROR: proc return code ", proc.returncode, "\n");
              # FIXME: should restart if it is not the same restart stage as before - need to parse output to 
              # find that value
              if completed and options.auto_resume:
                  cmd.append('--restart')
              else:
                  if options.auto_resume:
                      print("Could not restart, exiting...")
              return 1
          else:
              #print("SUCCESS: proc return code ", proc.returncode, "\n");
              break
      except:
          if proc:
              try:
                  print("Terminating process after exception: ", sys.exc_info(), "\n")
                  traceback.print_tb(sys.exc_info()[2], limit=100)
                  proc.terminate()
              except OSError:
                  pass
              except:
                  print("Unexpected error in final except: ", sys.exc_info())
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
        print("\nCaught an exception %s in mhmxx.py!\n\n" % e, file=sys.stderr); 
        traceback.print_exc(file=sys.stderr)
    finally:
        exit_all(status)

