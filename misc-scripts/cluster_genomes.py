#!/usr/bin/env python3

import pandas
import matplotlib
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


def test_model(model_data, model, genfracs):
    model_name = type(model).__name__
    
    bins = model.fit_predict(model_data)
    clusters = numpy.unique(bins)

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
        plt.scatter(series_x[i], series_y[i], lw=0.1, edgecolor='black',marker='.', label='cluster ' + str(i), color=colors[i])

    #plt.colorbar(ticks=range(len(clusters))).ax.set_yticklabels(clusters, fontsize=3)
    plt.xlabel('Depth')
    #plt.ylabel('TriNF entropy')
    plt.ylabel('GC count')
    #plt.xlim([min_x, 12.5])
    plt.tight_layout()
    plt.savefig('fig-clustering-' + model_name + '-' + fname + '.pdf')
    #plt.savefig('fig-clustering-' + fname + '.png', dpi=300)

    # now compute the accuracy: number of misclassifications, number of failed classifications
    bin_to_genome = {}
    for i, bin in enumerate(bins):
        if not bin in bin_to_genome:
            bin_to_genome[bin] = []
        bin_to_genome[bin].append([data['name'][i], data['genfrac'][i]])

    genome_to_bin_counts = {}
    num_errors = 0
    genfrac_error = 0

    bins_for_output = []
    for k, v in bin_to_genome.items():
        genomes_in_bin = set([x[0] for x in v])
        genfrac_in_bin = [x[1] for x in v]
        num_errors += (len(genomes_in_bin) - 1)
        #genfrac_error += genfracs[k]
        bins_for_output.append([k, sum(genfrac_in_bin), ', '.join(genomes_in_bin)])
        for genome in genomes_in_bin:
             if genome not in genome_to_bin_counts:
                 genome_to_bin_counts[genome] = 0
             genome_to_bin_counts[genome] += 1

    num_ungrouped = 0
    for val in genome_to_bin_counts.values():
        num_ungrouped += val - 1

    print(model_name + ', dimensions', len(model_data[0]))
    for bin in sorted(bins_for_output):
        print('  %4d %.2f' % (bin[0], bin[1]), bin[2])
    print('  found', len(clusters), 'clusters')
    print('  misclassifications:', num_errors, '(%.3f %%)' % (100.0 * num_errors / len(bins)))
    print('  ungrouped:', num_ungrouped, '(%.3f %%)' % (100.0 * num_ungrouped / len(bins)), flush=True)
    print('  genfrac error:', genfrac_error)

    
    

plt.style.use('qpaper')
fname = sys.argv[1]

# contents of file:

#bin name refdepth genfrac bin num_bins clen depth aln_depth entropy3k entropy2k gc_count
data = pandas.read_csv(fname, delim_whitespace=True)
kmer_depths = normalize(data['depth'])
aln_depths = normalize(data['aln_depth'])
entropy_3nfs = normalize(data['entropy3k'])
gc_counts = normalize(data['gc_count'])

model_data_3d = list(map(list, zip(aln_depths, gc_counts, entropy_3nfs)))
# well, it transpires that the entropy is really rather poor
#model_data_2d = list(map(list, zip(aln_depths, entropy_3nfs)))
model_data_2d = list(map(list, zip(aln_depths, gc_counts)))

print('There are', len(set(data['name'])), 'reference genomes')

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

#test_model(model_data_3d, sklearn.mixture.GaussianMixture(n_components=24, n_init=100))
test_model(model_data_2d, sklearn.mixture.GaussianMixture(n_components=24, n_init=10), genfracs=data['genfrac'])

#test_model(model_data_3d, cluster.Birch(n_clusters=25, threshold=0.04))
#test_model(model_data_2d, cluster.Birch(n_clusters=25, threshold=0.03))
