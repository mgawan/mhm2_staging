#!/usr/bin/env python3

import os
import sys
import shutil
import pandas
import subprocess
import glob
import time


def eval_clusters():
    start_t = time.time()
    num_clusters = len(glob.glob('cluster_*.fasta'))

    print('Found', num_clusters, 'clusters')

    print('Running metaquast')
    pids = []
    for i in range(num_clusters):
        mq_output = 'mq-cluster_' + str(i)
        shutil.rmtree(mq_output, ignore_errors=True)
        pids.append(subprocess.Popen(['metaquast.py', '--rna-finding', '--no-icarus', '--fragmented', '-t', '80', '-o', mq_output,
                                      '--fast', '-r', '/scratch2/shofmeyr/cami-data/mar_ref/arctic_samples/source_genomes/',
                                      'cluster_' + str(i) + '.fasta'], stdout=subprocess.DEVNULL))
    for pid in pids:
        pid.wait()


    bin_to_genomes = {}

    # extract genome fractions from mq results
    for i in range(num_clusters):
        fname = 'mq-cluster_' + str(i) + '/summary/TXT/Genome_fraction.txt'
        data = pandas.read_csv(fname, delim_whitespace=True)
        genomes = list(data['Assemblies'])
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
    for bin, genomes in bin_to_genomes.items():
        #print(bin, genomes)
        sorted_genomes = sorted(genomes, key=lambda x: x[1])
        genome_to_genfracs = {}
        for genome, genfrac in sorted_genomes:
            if genome not in genome_to_genfracs:
                genome_to_genfracs[genome] = 0.0
            genome_to_genfracs[genome] += float(genfrac)
        sorted_genomes_genfracs = sorted(genome_to_genfracs.items(), key=lambda x: x[1], reverse=True)
        print(bin, sorted_genomes_genfracs)
        primary_genome, primary_genfrac = sorted_genomes_genfracs[0]
        if primary_genome in primary_genomes:
            if primary_genomes[primary_genome] < primary_genfrac:
                primary_genomes[primary_genome] = primary_genfrac
        else:
            primary_genomes[primary_genome] = primary_genfrac;
        if len(sorted_genomes_genfracs) > 1:
            genfrac_error += sum([x[1] for x in sorted_genomes_genfracs[1:]])

    coverage = 0.0
    for genome, genfrac in primary_genomes.items():
        coverage += genfrac

    print('Average genome fraction errors: %.3f %%' % (genfrac_error / len(bin_to_genomes)))
    print('Average genome fraction in primary clusters: %.3f %%' % (coverage / len(bin_to_genomes)))
    print('Evaluation took %.3f s' % (time.time() - start_t))

if __name__ == "__main__":
    eval_clusters()
