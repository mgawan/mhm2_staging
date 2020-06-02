#!/usr/bin/env python3

import pandas
import matplotlib
import matplotlib.pyplot as plt
import numpy
import sys
import sklearn.mixture
import sklearn.cluster as cluster
#import sklearn.preprocessing 
import Bio.SeqIO as SeqIO


def normalize(xs, max_val=1.0):
    min_x = min(xs)
    max_x = max(xs)
    x_range = max_x - min_x
    return [max_val * (x - min_x) / x_range for x in xs]


def test_model(model_data, model):
    model_name = type(model).__name__
    
    bins = model.fit_predict(model_data)
    clusters = numpy.unique(bins)

    print(model_name + ', dimensions', len(model_data[0]))
    print('  found', len(clusters), 'clusters')

    series_x = [[] for i in range(len(clusters))]
    series_y = [[] for i in range(len(clusters))]

    for i, d in enumerate(model_data):
        series_x[bins[i]].append(d[0])
        series_y[bins[i]].append(d[1])

    cmap = matplotlib.cm.get_cmap('gist_ncar')
    colors = [cmap(float(i+2)/len(clusters)) for i in range(len(clusters))]

    min_x = 10000
    for i in range(len(series_x)):
        min_x = min(min_x, min(series_x[i]))
        plt.scatter(series_x[i], series_y[i], lw=0.1, edgecolor='black', marker='.', label='cluster ' + str(i), color=colors[i])

    #plt.colorbar(ticks=range(len(clusters))).ax.set_yticklabels(clusters, fontsize=3)
    plt.xlabel('Depth')
    #plt.ylabel('TriNF entropy')
    plt.ylabel('GC count')
    #plt.xlim([min_x, 12.5])
    plt.tight_layout()
    plt.savefig('fig-clustering-' + model_name + '-' + fname + '.pdf')
    #plt.savefig('fig-clustering-' + fname + '.png', dpi=300)

    return bins


plt.style.use('qpaper')
fname = sys.argv[1]

# contents of file:
#bin num_bins clen depth aln_depth entropyTNF gc_count A C G T
data = pandas.read_csv(fname, delim_whitespace=True)
kmer_depths = normalize(data['depth'])
# weight depths a bit
aln_depths = normalize(data['aln_depth'], 5)
entropy_TNFs = normalize(data['entropyTNF'])
gc_counts = list(data['gc_count'])
a_counts = list(data['A'])
c_counts = list(data['C'])
g_counts = list(data['G'])
t_counts = list(data['T'])
bin_idxs = list(data['bin'])

model_data_nuc_freqs = list(map(list, zip(aln_depths, a_counts, c_counts, g_counts, t_counts, entropy_TNFs)))
model_data_3d = list(map(list, zip(aln_depths, gc_counts, entropy_TNFs)))
# well, it transpires that the entropy is really rather poor
#model_data_2d = list(map(list, zip(aln_depths, entropy_TNFs)))
model_data_2d = list(map(list, zip(aln_depths, gc_counts)))

# from https://scikit-learn.org/stable/modules/clustering.html

#test_model(model_data_3d, cluster.KMeans(n_clusters=25))
#test_model(model_data_2d, cluster.KMeans(n_clusters=25))

#test_model(model_data_3d, cluster.AffinityPropagation(damping=0.6))
#test_model(model_data_2d, cluster.AffinityPropagation(damping=0.6))

#test_model(model_data_3d, cluster.MeanShift(bandwidth=0.06))
#test_model(model_data_2d, cluster.MeanShift(bandwidth=0.04))

#test_model(model_data_3d, cluster.SpectralClustering(n_clusters=25))
#test_model(model_data_2d, cluster.SpectralClustering(n_clusters=25))

#test_model(model_data_3d, cluster.AgglomerativeClustering(n_clusters=25))
#test_model(model_data_2d, cluster.AgglomerativeClustering(n_clusters=25))

#test_model(model_data_3d, cluster.DBSCAN(eps=0.1))
#test_model(model_data_2d, cluster.DBSCAN(eps=0.05))

#test_model(model_data_3d, cluster.OPTICS())
#test_model(model_data_2d, cluster.OPTICS())

#test_model(model_data_nuc_freqs, sklearn.mixture.GaussianMixture(n_components=24, n_init=100))
#test_model(model_data_3d, sklearn.mixture.GaussianMixture(n_components=24, n_init=100))
#test_model(model_data_2d, sklearn.mixture.GaussianMixture(n_components=12, n_init=100))
#bins = test_model(model_data_3d, sklearn.mixture.GaussianMixture(n_components=25, n_init=100))
bins = test_model(model_data_2d, sklearn.mixture.GaussianMixture(n_components=25, n_init=100))
#test_model(model_data_2d, sklearn.mixture.GaussianMixture(n_components=100, n_init=100))

#test_model(model_data_3d, cluster.Birch(n_clusters=25, threshold=0.04))
#test_model(model_data_2d, cluster.Birch(n_clusters=25, threshold=0.03))

clusters = {}
print('Getting clusters', flush=True)
for i, bin in enumerate(bins):
#    print(i, bin_idxs[i], bin)
    if bin not in clusters:
        clusters[bin] = []
    for record in SeqIO.parse('bin_' + str(bin_idxs[i]) + '.fasta', "fasta"):
        clusters[bin].append(record)
print('Writing', len(clusters), 'clusters', flush=True)
for bin, cluster in clusters.items():
#    print(bin, len(cluster))
    SeqIO.write(cluster, 'cluster_' + str(bin) + '.fasta', 'fasta')

    
