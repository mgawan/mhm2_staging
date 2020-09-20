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

    options = argparser.parse_args()

    quals = get_qual_vals(options.quals_fname)
    # first, run metaquast on the assembly
    os.chdir(options.asm_dir)
    cmd = ['metaquast.py', '--fast', '-o', 'mq.out', '-r', options.refs, 'final_assembly.fasta']
    test_exec_mq = which('metaquast.py')
    if not test_exec_mq:
        test_exec_shifter = which('shifter')
        # test_exec_docker = which('docker')
        if test_exec_shifter:
            shifter = ['shifter', '--image=robegan21/quast:latest']
            shifter.extend(cmd)
            cmd = shifter
    print('Running metaquast...', cmd)
    subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    new_quals = get_qual_vals('mq.out/combined_reference/report.txt')
    num_mismatches = 0
    for key, val in quals.items():
        if key not in new_quals:
            print("ERROR: '" + key + "' not found")
        else:
            d = 0
            thres = options.thres
            if val != new_quals[key]:
                d = abs(float(val) - float(new_quals[key])) / max(float(val), float(new_quals[key]))
                if (key.startswith('# contigs') or key.startswith('# misassemblies') or key.startswith('# misassembled contigs') \
                    or key.startswith('# local misassemblies')) and d < 4:
                    d = 0
                if key.startswith('Total length (>= 10000 bp)') or key.startswith('Total length (>= 25000 bp)') \
                   or key.startswith('Misassembled contigs length') or key.startswith('# indels per 100 kbp'):
                    thres = 0.05
                if key.startswith("# N's per 100 kbp") and d < 1.5:
                    d = 0
            if d > thres:
                print("MISMATCH:", key, val, "!=", new_quals[key], 'd = %.3f' % d)
                num_mismatches += 1
            #print(key, val, d)
    print("Comparison yielded", num_mismatches, "mismatches")
    
if __name__ == "__main__":
    main()

