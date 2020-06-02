#!/usr/bin/env python3

import pandas
import matplotlib.pyplot as plt
import numpy
import sys
import matplotlib

def normalize(xs):
    min_x = min(xs)
    max_x = max(xs)
    x_range = max_x - min_x
    return [(x - min_x) / x_range for x in xs]


plt.style.use('qpaper')

fname = sys.argv[1]
data = pandas.read_csv(fname, delim_whitespace=True)
    
#bin1 name refdepth genfrac bin num_bins clen depth aln_depth entropyTNF gc_count A C G T
aln_depths = list(data['aln_depth'])
entropy_TNF = list(data['entropyTNF'])
gc_counts = list(data['gc_count'])

genome_depths = {}
ref_depths = list(data['refdepth'])
for i, genome in enumerate(list(data['name'])):
    if genome not in genome_depths:
        genome_depths[genome] = ref_depths[i]
    
genomes = list(set(data['name']))
genomes_and_depths = []
for genome in genomes:
    genomes_and_depths.append(genome + ' ' + str(genome_depths[genome]))

genome_colors = {}
for i, genome in enumerate(genomes):
    genome_colors[genome] = i + 1

points_colors = []
for genome in list(data['name']):
    points_colors.append(genome_colors[genome])


#plt.scatter(aln_depths, entropy_TNF, lw=0.1, edgecolor='black', marker='.', c=points_colors, cmap='nipy_spectral')
plt.scatter(aln_depths, gc_counts, lw=0.1, edgecolor='black', marker='.', c=points_colors, cmap='nipy_spectral')
plt.colorbar(ticks=range(len(genomes))).ax.set_yticklabels(genomes_and_depths, fontsize=3)
plt.xlabel('Depth')
#plt.ylabel('TriNF entropy')
plt.ylabel('GC count')
#plt.xlim([min(depths), 12.5])
plt.tight_layout()
#plt.savefig('fig-depth-4k-' + fname + '.png', dpi=200)
plt.savefig('fig-depth-3k-' + fname + '.pdf')

