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
data = pandas.read_csv(fname, delim_whitespace=True, header=None)
    
#depths = normalize(data[7])
#tnfs = normalize(data[9])
depths = list(data[8])
tnfs = list(data[9])

genome_depths = {}
ref_depths = list(data[2])
for i, genome in enumerate(list(data[1])):
    if genome not in genome_depths:
        genome_depths[genome] = ref_depths[i]
    
genomes = list(set(data[1]))
genomes_and_depths = []
for genome in genomes:
    genomes_and_depths.append(genome + ' ' + str(genome_depths[genome]))

genome_colors = {}
for i, genome in enumerate(genomes):
    genome_colors[genome] = i + 1

points_colors = []
for genome in list(data[1]):
    points_colors.append(genome_colors[genome])


plt.scatter(depths, tnfs, lw=0, marker='.', c=points_colors, cmap='nipy_spectral')
plt.colorbar(ticks=range(len(genomes))).ax.set_yticklabels(genomes_and_depths, fontsize=3)
plt.xlabel('Depth')
plt.ylabel('4-mer entropy')
#plt.xlim([min(depths), 12.5])
plt.tight_layout()
plt.savefig('fig-depth-4k-' + fname + '.png', dpi=200)

