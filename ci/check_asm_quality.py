#!/usr/bin/env python3

import subprocess
import sys
import os
import datetime
import argparse
import string


def get_qual_vals(fname):
    quals = {}
    with open(fname) as f:
        for line in f.readlines():
            if line.startswith('All') or line.startswith('Assembly') or len(line) < 2:
                continue
            check_multi = line.find(' + ')
            if check_multi > 0:
                line = line[:check_multi]
            fields = line.split()
            val = fields[-1]
            key = line[:(len(line.strip()) - len(val))]
            quals[key] = val
    return quals
    
def which(file_name):
    if os.path.exists(file_name) and os.access(file_name, os.X_OK):
        return file_name
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None

def main():
    argparser = argparse.ArgumentParser(add_help=False)
    argparser.add_argument("--asm-dir", dest='asm_dir', type=str, required=True, help='Directory containing assembly')
    argparser.add_argument("--expected-quals", dest='quals_fname', type=str, required=True, help='File containing expected metaquast qualities')
    argparser.add_argument("--refs", dest='refs', type=str, required=True, help='Directory containing references for metaquast')
    argparser.add_argument("--thres", dest='thres', type=float, default=0.01, help='Threshold of minimum difference, fraction')
    argparser.add_argument("--rna", action='store_true', default=False, help='Include rna predictions')

    options = argparser.parse_args()

    quals = get_qual_vals(options.quals_fname)
    # first, run metaquast on the assembly
    report_path = 'mq.out/combined_reference/report.txt'
    report_exists = os.path.exists(options.asm_dir + "/" + report_path)
    options.refs = os.path.realpath(options.refs)
    os.chdir(options.asm_dir)

    if not report_exists:
        pwd=os.getcwd()
        cmd = ['metaquast.py', '--fast', '-o', '%s/mq.out'%(pwd), '-r', options.refs, '%s/final_assembly.fasta'%(pwd)]
        if options.rna:
            cmd.append('--rna-finding')
        orig_mq_cmd = cmd
        test_exec_mq = which('metaquast.py')
        test_exec_shifter = which('shifter')
        test_exec_docker = which('docker')
        if test_exec_shifter:
            shifter = ['shifter', '--image=robegan21/quast:latest']
            shifter.extend(cmd)
            cmd = shifter
        elif test_exec_docker:
            user=os.getuid()
            refpath=os.path.dirname(options.refs)
            docker = ['docker', 'run', '-i', '--tty=false', '-a', 'STDIN', '-a', 'STDOUT', '-a', 'STDERR', '--user', '%s:%s' %(user,user), '--volume=%s:%s' % (refpath,refpath), '--volume=%s:%s' % (pwd,pwd), '--workdir=%s' % (pwd), 'robegan21/quast:latest']
            docker.extend(cmd)
            cmd = docker
        elif not test_exec_mq:
            sys.exit('ERROR: requires shifter, docker or metaquast.py in the path to check')
        while True:
            print('Running metaquast:', cmd)
            try:
                subprocess.check_output(cmd)#, stderr=subprocess.STDOUT)
                break
            except (subprocess.CalledProcessError):
                if cmd == orig_mq_cmd:
                    raise
                elif test_exec_docker and test_exec_mq:
                    print('Docker failed, trying fallback to metaquast.py directly')
                    cmd = orig_mq_cmd
                else:
                    raise
            
    new_quals = get_qual_vals(report_path)
    num_mismatches = 0
    for key, val in quals.items():
        if key not in new_quals:
            if 'rRNA' in key:
                print("WARN: '" + key + "' not found")
            else:
                print("ERROR: '" + key + "' not found")
        else:
            d = 0
            thres = options.thres
            if val != new_quals[key]:
                d_real = float(new_quals[key]) - float(val)
                d_abs = abs(d_real)
                d_max = max(float(val), float(new_quals[key]))
                d = d_abs / d_max
                if (key.startswith('# contigs') or key.startswith('# misassemblies') or key.startswith('# misassembled contigs') \
                    or key.startswith('# local misassemblies')) and d > thres and d < 4:
                    print("WARN: adjusted threshold: ", key, val, "!=", new_quals[key], 'd = %.3f' % d)
                    d = 0
                if key.startswith('Total length (>= 0 bp)') and d > thres:
                    print("WARN: adjusted threshold: ", key, val, "!=", new_quals[key], 'd = %.3f' % d)
                    thres = thres + 0.05 # + 5%
                    if d_real > 0 or d_abs < 3000: # < 3k diff on 0k total
                        d = 0                    
                if key.startswith('Total length (>= 1000 bp)') and d > thres:
                    print("WARN: adjusted threshold: ", key, val, "!=", new_quals[key], 'd = %.3f' % d)
                    thres = thres + 0.05 # + 5%
                    if d_real > 0 or d_abs < 3000: # < 3k diff on 1k total
                        d = 0
                if key.startswith('Total length (>= 5000 bp)') and d > thres:
                    print("WARN: adjusted threshold: ", key, val, "!=", new_quals[key], 'd = %.3f' % d)
                    thres = thres + 0.05 # + 5%
                    if d_real > 0 or d_abs < 10000: # < 10k diff on 5k total
                        d = 0    
                if key.startswith('Total length (>= 10000 bp)') and d > thres:
                    print("WARN: adjusted threshold: ", key, val, "!=", new_quals[key], 'd = %.3f' % d)
                    thres = thres + 0.125 # + 12.5%
                    if d_real > 0 or d_abs < 20000: # < 20k diff on 10k total
                        d = 0
                if key.startswith('Total length (>= 25000 bp)') and d > thres:
                    print("WARN: adjusted threshold: ", key, val, "!=", new_quals[key], 'd = %.3f' % d)
                    thres = thres + 0.075 # + 7.5%
                    if d_real > 0 or d_abs < 50000: # < 50k diff on 25k bp total
                        d = 0
                if key.startswith('Total length (>= 50000 bp)') and d > thres:
                    print("WARN: adjusted threshold: ", key, val, "!=", new_quals[key], 'd = %.3f' % d)
                    thres = thres + 0.075 # 7.5%
                    if d_real > 0 or d_abs < 100000: # < 100k diff on 50k bp total
                        d = 0
                if key.startswith("# N's per 100 kbp") and d > thres:
                    print("WARN: adjusted threshold: ", key, val, "!=", new_quals[key], 'd = %.3f' % d)
                    thres = thres + 0.05 # 5%
                    if d_real < 0 or d_abs < 3: # < 3 per 100k diff
                        d = 0
                if key.startswith("# indels per 100 kbp") and d > thres:
                    print("WARN: adjusted threshold: ", key, val, "!=", new_quals[key], 'd = %.3f' % d)
                    thres = thres + 0.05 # 5%
                    if d_real < 0 or d_abs < 3: # < 3 per 100k diff
                        d = 0
                if key.startswith("Misassembled contigs length") and d > thres:
                    print("WARN: adjusted threshold: ", key, val, "!=", new_quals[key], 'd = %.3f' % d)
                    thres = thres + 0.25 # + 25%
                    if d_real < 0:
                        d = 0
                if key.startswith("# predicted rRNA genes") and d > thres:
                    print("WARN: adjusted threshold: ", key, val, "!=", new_quals[key], 'd = %.3f' % d)
                    thres = thres + 0.05 # + 5%
                    if d_real > 0 or d_abs < 2: # < 2 diff
                        d = 0
                if key.startswith("# mismatches per 100 kbp") and d > thres:
                    print("WARN: adjusted threshold: ", key, val, "!=", new_quals[key], 'd = %.3f' % d)
                    thres = thres + 0.05 # + 5%
                    if d_real < 0 or d_abs < 2: # < 2 diff
                        d = 0
                if key.startswith("NA50") and d > thres:
                    print("WARN: adjusted threshold: ", key, val, "!=", new_quals[key], 'd = %.3f' % d)
                    thres = thres + 0.025 # + 2.5%
            if d > thres:
                print("MISMATCH:", key, val, "!=", new_quals[key], 'd = %.3f' % d, ' thres = %0.3f' % thres)
                num_mismatches += 1
            #print(key, val, d)
    print("Comparison yielded", num_mismatches, "mismatches")
    if num_mismatches > 0: 
        sys.exit("Detected %d mismatches!" % (num_mismatches)) # exit non-zero
    
if __name__ == "__main__":
    main()

