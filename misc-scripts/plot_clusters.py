#!/usr/bin/env python3

import pandas
import matplotlib.pyplot as plt
import numpy
import sys
import argparse
import matplotlib

def normalize(xs):
    min_x = min(xs)
    max_x = max(xs)
    x_range = max_x - min_x
    return [(x - min_x) / x_range for x in xs]


def plot_clusters(cluster_stats_fname, mapped_refs_fname, refs_depths_fname, output_fname):
    plt.style.use('qpaper')

    cluster_data = pandas.read_csv(cluster_stats_fname, delim_whitespace=True)
    cluster_data.sort_values('bin')
    bins = list(cluster_data['bin'])
    #bin num_bins clen depth aln_depth entropyTNF gc_count A C G T
    aln_depths = list(cluster_data['aln_depth'])
    entropy_TNF = list(cluster_data['entropyTNF'])
    gc_counts = list(cluster_data['gc_count'])

    points_colors = []
    genomes_depths_strs = []
    if mapped_refs_fname != '' and refs_depths_fname != '':
        print('Using refs names and depths for colorbar')
        refs_depths_data = pandas.read_csv(refs_depths_fname, delim_whitespace=True)
        ref_names = list(refs_depths_data['genome'])
        ref_depths = list(refs_depths_data['depth'])
        genome_depths = {}
        for i, genome in enumerate(ref_names):
            if genome not in genome_depths:
                genome_depths[genome] = ref_depths[i]

        for genome in sorted(ref_names, reverse=True):
            genomes_depths_strs.append(genome + ' ' + str(genome_depths[genome]))
        genome_colors = {}
        for i, genome in enumerate(sorted(ref_names, reverse=True)):
            genome_colors[genome] = i

        #tab delimited!
        #bin     genfrac genfrac_err     refs_list
        refs_data = pandas.read_csv(mapped_refs_fname, sep='\t')
        refs_lists = list(refs_data['refs_list'])
        
        for refs_list in refs_lists:
            # choose the most abundant in a list - always the first one
            genome = refs_list.split()[0]
            points_colors.append(genome_colors[genome])
            
    else:
        points_colors = range(1, len(bins))    

    marker_sizes = [3] * len(aln_depths)
    plt.scatter(aln_depths, gc_counts, marker_sizes, lw=0.05, edgecolor='black', marker='.', c=points_colors, cmap='nipy_spectral')
    if mapped_refs_fname != '' and refs_depths_fname != '':
        plt.colorbar(ticks=range(len(genomes_depths_strs))).ax.set_yticklabels(genomes_depths_strs, fontsize=3)
    plt.xlabel('Depth')
    plt.ylabel('GC count')
    #plt.xlim([min(depths), 12.5])
    plt.tight_layout()
    print('Saving plot to', output_fname)
    plt.savefig(output_fname, dpi=300)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot clusters')
    parser.add_argument('-f', dest='cluster_stats_fname', required=True,
                        help='File containing cluster information - from bin_gfa.py (stats_bins.txt)')
    parser.add_argument('-r', dest='mapped_refs_fname', 
                        help='File containing references mapped to clusters - from eval_clusters.py (eg bin-eval.txt)')
    parser.add_argument('-d', dest='refs_depths_fname',
                        help='File containing references with original depths (eg arctic-samples-depths.txt')
    parser.add_argument('-o', dest='output_fname', required=True,
                        help='File with the output plot')
    opts = parser.parse_args()
    plot_clusters(opts.cluster_stats_fname, opts.mapped_refs_fname, opts.refs_depths_fname, opts.output_fname)
