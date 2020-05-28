#!/usr/bin/env python3

import pandas
import matplotlib.pyplot as plt
import numpy
import sys
import sklearn.mixture
import sklearn.cluster as cluster
#import sklearn.preprocessing 


def normalize(xs):
    min_x = min(xs)
    max_x = max(xs)
    x_range = max_x - min_x
    return [(x - min_x) / x_range for x in xs]

plt.style.use('qpaper')
fname = sys.argv[1]

# contents of file:
# bini name ref_depth genfrac bini num_bins clen depth aln_depth entropy3k entropy2k gc_count

all_data = pandas.read_csv(fname, delim_whitespace=True, header=None)
kmer_depths = normalize(all_data[7])
kmer_depths = normalize(all_data[7])
tnfs = normalize(all_data[11])
#kmer_depths = normalize(all_data[7])
#depths = list(all_data[8])
#tnfs = list(all_data[9])

#data = list(map(list, zip(depths, tnfs, kmer_depths)))
data = list(map(list, zip(depths, tnfs)))

# from https://scikit-learn.org/stable/modules/clustering.html

# requires specifying the number of clusters. Not great either
#model = cluster.KMeans(n_clusters=25)

# good but too few clusters
#model = cluster.AffinityPropagation()

# not bad
#model = cluster.MeanShift(bandwidth=0.03)

# requires setting the number of clusters - seems really poor
#model = cluster.SpectralClustering(n_clusters=25)

# not bad
#model = cluster.AgglomerativeClustering(n_clusters=None, distance_threshold=0.15)

# hopeless
#model = cluster.DBSCAN(eps=0.001)

# not terrible, but not so great either
#model = cluster.OPTICS()

# quite decent, but requires specifying the number of clusters
model = sklearn.mixture.GaussianMixture(n_components=25)

# ok, but not great
#model = cluster.Birch(n_clusters=None, threshold=0.04)


bins = model.fit_predict(data)
clusters = numpy.unique(bins)
print('Found', len(clusters), 'clusters')
series_x = [[] for i in range(len(clusters))]
series_y = [[] for i in range(len(clusters))]

for i, d in enumerate(data):
    series_x[bins[i]].append(d[0])
    series_y[bins[i]].append(d[1])

#colors = ['red', 'green', 'blue', 'cyan', 'magenta', 'orange', 'yellow', 'black', 'brown', 'lightgrey']    
min_x = 10000
for i in range(len(series_x)):
    min_x = min(min_x, min(series_x[i]))
    plt.scatter(series_x[i], series_y[i], lw=0, marker='.', label='cluster ' + str(i), cmap='gist_ncar')

#plt.colorbar(ticks=range(len(clusters))).ax.set_yticklabels(clusters, fontsize=3)
plt.xlabel('Depth')
plt.ylabel('4-mer entropy')
#plt.xlim([min_x, 12.5])
plt.tight_layout()
#plt.savefig('fig-clustering-' + model.__name__ + '-' + fname + '.pdf')
plt.savefig('fig-clustering-' + fname + '.png', dpi=300)


