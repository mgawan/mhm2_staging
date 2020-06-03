#!/usr/bin/env python3

import os
import sys
import shutil
import pandas
import subprocess
import glob
import time


def eval_clusters(fname_prefix):
    start_t = time.time()
    num_clusters = len(glob.glob(fname_prefix + '_*.fasta'))

    print('Found', num_clusters, 'clusters', file=sys.stderr)
    print('Running metaquast', file=sys.stderr)
    pids = []
    for i in range(num_clusters):
        mq_output = 'mq-cluster_' + str(i)
        shutil.rmtree(mq_output, ignore_errors=True)
        pids.append(subprocess.Popen(['metaquast.py', '--rna-finding', '--no-icarus', '--fragmented', '-t', '80', '-o', mq_output,
                                      '--fast', '-r', '/scratch2/shofmeyr/cami-data/mar_ref/arctic_samples/source_genomes/',
                                      fname_prefix + '_' + str(i) + '.fasta'], stdout=subprocess.DEVNULL))
    for pid in pids:
        pid.wait()


    bin_to_genomes = {}
    num_genomes = 0
    # extract genome fractions from mq results
    for i in range(num_clusters):
        fname = 'mq-' + fname_prefix +'_' + str(i) + '/summary/TXT/Genome_fraction.txt'
        data = pandas.read_csv(fname, delim_whitespace=True)
        genomes = list(data['Assemblies'])
        num_genomes = len(genomes)
        genfracs = list(data['cluster_' + str(i)])
        for j, genfrac in enumerate(genfracs):
            if genfrac != '-':
                if i not in bin_to_genomes:
                    bin_to_genomes[i] = []
                bin_to_genomes[i].append((genomes[j], genfrac))


    # The primary genome for a cluster is the one with the highest genfrac in that cluster.
    # A cluster is a primary cluster if the primary genome in it has the highest genfrac of any clusters
    # containing that genome.
    #
    # So we can calculate two metrics of accuracy:
    # - Average genome fraction that are primary in primary clusters (coverage)
    # - Average genome fraction not primaries or not in primary clusters (errors)

    primary_genomes = {}
    genome_to_bin_counts = {}
    num_errors = 0
    genfrac_error = 0
    num_pure = 0
    num_high_qual = 0
    num_medium_high_qual = 0
    num_medium_qual = 0
    all_results = []
    for bin, genomes in bin_to_genomes.items():
        sorted_genomes = sorted(genomes, key=lambda x: x[1])
        genome_to_genfracs = {}
        for genome, genfrac in sorted_genomes:
            if genome not in genome_to_genfracs:
                genome_to_genfracs[genome] = 0.0
            genome_to_genfracs[genome] += float(genfrac)
        sorted_genomes_genfracs = sorted(genome_to_genfracs.items(), key=lambda x: x[1], reverse=True)
        primary_genome, primary_genfrac = sorted_genomes_genfracs[0]
        if primary_genome in primary_genomes:
            if primary_genomes[primary_genome] < primary_genfrac:
                primary_genomes[primary_genome] = primary_genfrac
        else:
            primary_genomes[primary_genome] = primary_genfrac;
        this_genfrac_error = 0
        if len(sorted_genomes_genfracs) == 1:
            num_pure += 1
        else:
            this_genfrac_error = sum([x[1] for x in sorted_genomes_genfracs[1:]])
        genfrac_error += this_genfrac_error
        if primary_genfrac >= 90 and this_genfrac_error < 5:
            num_high_qual += 1
        elif primary_genfrac >= 70 and this_genfrac_error < 10:
            num_medium_high_qual += 1
        elif primary_genfrac >= 50 and this_genfrac_error < 10:
            num_medium_qual += 1        
        all_results.append((bin, primary_genfrac, this_genfrac_error, sorted_genomes_genfracs))

    coverage = 0.0
    for genome, genfrac in primary_genomes.items():
        coverage += genfrac
        
    for result in sorted(all_results, key=lambda x: x[1], reverse=True):
        print(result[0], result[1], '%.2f' % result[2], result[3])
        
    print('Number of pure clusters:', num_pure, '( %.2f %%)' % (float(num_pure) / len(bin_to_genomes)))
    print('Quality of clusters:')
    print('    high (frac > 90%, error < 5%)', num_high_qual)
    print('    medium high (frac > 70%, error < 10%)', num_medium_high_qual)
    print('    medium quality (frac > 50%, error < 10%)', num_medium_qual)
    print('Average genome fraction errors: %.3f %%' % (genfrac_error / num_genomes))
    print('Average genome fraction in primary clusters: %.3f %%' % (coverage / num_genomes))
    print('Evaluation took %.3f s' % (time.time() - start_t))

if __name__ == "__main__":
    fname_prefix = sys.argv[1]
    eval_clusters(fname_prefix)
