#!/usr/bin/env python3

import sys
import time
import math
import argparse
import Bio.SeqIO as SeqIO
import networkx as nx
import networkx.algorithms.community as cm
import networkx.algorithms.components as cp
import scipy.stats


parser = argparse.ArgumentParser(description='Find bins using GFA output')
parser.add_argument('-f', dest='fasta_fname', help='file containing fasta assembly')
parser.add_argument('-g', dest='gfa_fname', help='file containing GFA2 format corresponding to assembly')
parser.add_argument('-t', dest='cum_len_thres', type=int, default=1000000, help='Cumulative contig length to accept a bin')    

opts = parser.parse_args()


def analyze_bins(g, bins):
    bin_counts = {}
    for i, bin in enumerate(bins):
        tot_clen = 0
        tot_depth = 0.0
        tot_aln_depth = 0.0
        for cid in bin:
            clen = g.nodes[cid]['clen']
            tot_clen += clen
            tot_depth +=  g.nodes[cid]['depth'] * clen
            tot_aln_depth +=  g.nodes[cid]['aln_depth'] * clen
        bin_counts[i] = (len(bin), tot_clen, tot_depth / tot_clen, tot_aln_depth / tot_clen)

    count_thres = [1000000, 500000, 200000, 100000, 50000, 20000, 10000, 0]
    counts = [0] * len(count_thres)
    for bin_i in sorted(bin_counts.keys(), key=lambda x: bin_counts[x][1], reverse=True):
        #print(bin_counts[bin_i])
        for i, thres in enumerate(count_thres):
            if bin_counts[bin_i][1] >= thres:
                counts[i] += 1
                break

    print('  found', len(bins), 'bins')
    for i, thres in enumerate(count_thres):
        print('    clen >= %-10d %-10d' % (thres, counts[i]))
        
    return bins, bin_counts
        

def longest_edge(g):
    max_edge = (0, 0)
    max_edge_cclen = 0
    for edge in g.edges.data():
        if edge[2]['cclen'] > max_edge_cclen:
            max_edge = (edge[0], edge[1])
    return max_edge


def get_kmer_distr_entropy(records, k):
    n = 0
    tnfs = {}
    tot_count = 0
    tot_gc = 0.0
    nuc_counts = {}
    for record in records:
        for i in range(len(record.seq) - k):
            if record.seq[i] in ['G', 'C']:
                tot_gc += 1
            if record.seq[i] not in nuc_counts:
                nuc_counts[record.seq[i]] = 0
            nuc_counts[record.seq[i]] += 1
            tn = str(record.seq[i:i+k])
            if tn not in tnfs:
                tnfs[tn] = 0
            tnfs[tn] += 1
            tot_count += 1
    # normalize the tnfs
    for key, val in tnfs.items():
        tnfs[key] = float(val) / tot_count
    nuc_freqs = []
    for key in sorted(nuc_counts.keys()):
        nuc_freqs.append(float(nuc_counts[key]) / tot_count)
    # now compute shannon entropy (base 2) of the distribution
    return scipy.stats.entropy(list(tnfs.values()), base=2), tot_gc / tot_count, nuc_freqs
    

def main():
    g = nx.Graph()

    num_edges = 0
    num_vertices = 0
    with open(opts.gfa_fname) as f:
        for line in f.readlines():
            fields = line.split('\t')
            if fields[0][0] in ['E', 'G']:
                g.add_edge(int(fields[2][:-1]), int(fields[3][:-1]))
            elif fields[0][0] == 'S':
                tags = fields[4].split()
                g.add_node(int(fields[1]), clen=int(fields[2]), depth=float(tags[1]), aln_depth=float(tags[3]))

    # edge attribute is the sum of the two contig lengths
    for edge in g.edges.data():
        g[edge[0]][edge[1]]['cclen'] = g.nodes[edge[0]]['clen'] + g.nodes[edge[1]]['clen']

    print('Graph has', g.number_of_nodes(), 'vertices and', g.number_of_edges(), 'edges')

#    for method in [cp.connected_components, cm.greedy_modularity_communities, cm.label_propagation_communities, cm.girvan_newman]:
    for method in [cp.connected_components]:
        cid_to_bin = {}
        t = time.perf_counter()
        print(method.__name__, flush=True)
        if method == cm.girvan_newman:
            cms = method(g, most_valuable_edge=longest_edge)
            bins = [list(s) for s in next(cms)]
        else:
            cms = method(g)
            bins = [list(s) for s in cms]
        print('  Took %.4f s' % (time.perf_counter() - t))
        bins, bin_counts = analyze_bins(g, bins)
        num_bins = 0
        for i, bin in enumerate(bins):
            for cid in bin:
                cid_to_bin[cid] = i
            if bin_counts[i][1] >= opts.cum_len_thres:
                num_bins += 1

        print('Writing bins')
        bin_records = [[] for i in range(len(bins))]
        for record in SeqIO.parse(opts.fasta_fname, "fasta"):
            cid = int(record.id[len("Contig"):])
            if cid not in cid_to_bin:
                continue
            bin_records[cid_to_bin[cid]].append(record)
        num_bins_written = 0
        print('Computing TNF entropy for', num_bins, 'bins')
        entropiesTNF = {}
        gc_content = {}
        nuc_freqs = {}
        ten_perc = int(num_bins / 10)
        t = time.perf_counter()
        for i, ith_bin_records in enumerate(bin_records):
            if bin_counts[i][1] >= opts.cum_len_thres:
                SeqIO.write(ith_bin_records, 'bin_' + str(num_bins_written) + '.fasta', 'fasta')
                num_bins_written += 1
                # compute the TNF
                entropiesTNF[i], gc_content[i], nuc_freqs[i] = get_kmer_distr_entropy(ith_bin_records, 4)
                if num_bins_written % ten_perc == 0:
                    print(10 * int(num_bins_written / ten_perc), end='% ', flush=True)
        print('')
        print('Wrote', num_bins_written, 'bins with cumulative clen >=', opts.cum_len_thres)
        print('Took %.4f s' % (time.perf_counter() - t))
        print('Writing stats_bins.txt')
        with open('stats_bins.txt', 'w') as f:
            print('bin num_bins clen depth aln_depth entropyTNF gc_count A C G T', file=f)
            for i, bin in enumerate(bins):
                if bin_counts[i][1] >= opts.cum_len_thres:
                    print(i, bin_counts[i][0], bin_counts[i][1], '%.3f' % bin_counts[i][2], '%.3f' % bin_counts[i][3],
                          '%.4f' % entropiesTNF[i], '%.4f' % gc_content[i], end=' ', file=f)
                    for nuc_freq in nuc_freqs[i]:
                        print('%.4f' % nuc_freq, end=' ', file=f)
                    print('', file=f)
                for cid in bin:
                    cid_to_bin[cid] = i
        
                
if __name__ == "__main__":
    main()
