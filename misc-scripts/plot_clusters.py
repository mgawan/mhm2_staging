#!/usr/bin/env python3

import pandas
import matplotlib.pyplot as plt
import numpy
import sys
import argparse
import matplotlib


def normalize(xs, min_x, max_x, min_val=0.0, max_val=1.0):
    x_range = max_x - min_x
    return [(float(x) - min_x) / x_range * (max_val - min_val) + min_val for x in xs]


def plot_clusters(opts):
    plt.style.use('qpaper')

    cluster_data = pandas.read_csv(opts.cluster_stats_fname, delim_whitespace=True)
    cluster_data.sort_values('bin')
    bins = list(cluster_data['bin'])
    #bin num_bins clen depth aln_depth entropyTNF gc_count A C G T
    aln_depths = list(cluster_data['aln_depth'])
    entropy_TNF = list(cluster_data['entropyTNF'])
    gc_counts = list(cluster_data['gc_count'])
        
    points_colors = []
    genomes_depths_strs = []
    if opts.mapped_refs_fname != '' and opts.refs_depths_fname != '':
        print('Using refs names and depths for colorbar')
        refs_depths_data = pandas.read_csv(opts.refs_depths_fname, delim_whitespace=True)
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
        refs_data = pandas.read_csv(opts.mapped_refs_fname, sep='\t')
        refs_lists = list(refs_data['refs_list'])
        
        for refs_list in refs_lists:
            # choose the most abundant in a list - always the first one
            genome = refs_list.split()[0]
            points_colors.append(genome_colors[genome])
    else:
        points_colors = range(1, len(bins))    

    min_clen = min(cluster_data['clen'])
    max_clen = max(cluster_data['clen'])
    max_point_size = 100
    other_cluster_data = None
    if opts.other_cluster_stats_fname:
        other_cluster_data =  pandas.read_csv(opts.other_cluster_stats_fname, delim_whitespace=True)
        min_clen = min(min_clen, min(other_cluster_data['clen']))
        max_clen = max(max_clen, max(other_cluster_data['clen']))
    clens = normalize(cluster_data['clen'], min_clen, max_clen, 2, max_point_size)
    zipped_lists = zip(clens, aln_depths, gc_counts, points_colors)
    sorted_tuples = sorted(zipped_lists, reverse=True)
    tuples = zip(*sorted_tuples)
    clens, aln_depths, gc_counts, points_colors = [list(tuple) for tuple in tuples]
    
    plt.scatter(aln_depths, gc_counts, clens, lw=0.05, edgecolor='black', marker='.', c=points_colors, cmap='nipy_spectral',
                alpha=0.7)
    if opts.mapped_refs_fname != '' and opts.refs_depths_fname != '':
        plt.colorbar(ticks=range(len(genomes_depths_strs))).ax.set_yticklabels(genomes_depths_strs, fontsize=3)

    if other_cluster_data is not None:
        other_aln_depths = other_cluster_data['aln_depth']
        other_gc_counts = other_cluster_data['gc_count']
        other_clens = normalize(other_cluster_data['clen'], min_clen, max_clen, 2, max_point_size)
        colors = ['white'] * len(other_clens)
        plt.scatter(other_aln_depths, other_gc_counts, other_clens, c=colors, lw=0.2, edgecolor='black', marker='.', alpha=0.5)
    plt.xlabel('Depth')
    plt.ylabel('GC count')
    plt.tight_layout()
    print('Saving plot to', opts.output_fname)
    plt.savefig(opts.output_fname, dpi=300)

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot clusters')
    parser.add_argument('-f', dest='cluster_stats_fname', required=True,
                        help='File containing cluster information - from bin_gfa.py (eg stats_bins.txt)')
    parser.add_argument('-r', dest='mapped_refs_fname', 
                        help='File containing references mapped to clusters - from eval_clusters.py (eg bin-eval.txt)')
    parser.add_argument('-d', dest='refs_depths_fname',
                        help='File containing references with original depths (eg arctic-samples-depths.txt')
    parser.add_argument('-o', dest='output_fname', required=True,
                        help='File with the output plot')
    parser.add_argument('-b', dest='other_cluster_stats_fname',
                        help='File containing cluster information for another method - also from bin_gfa.py')
    opts = parser.parse_args()
    plot_clusters(opts)
