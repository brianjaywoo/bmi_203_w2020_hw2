import sys
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning, cluster_hierarchically, return_silhouette_score

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file>")
    sys.exit(0)

active_sites = read_active_sites(sys.argv[2])
clustering_type = sys.argv[1][0:2]

# Choose clustering algorithm
if clustering_type == '-P':
    print("Clustering using Partitioning method")
    clustering = cluster_by_partitioning(active_sites)
    silhouette_score = return_silhouette_score(clustering, clustering_type)
    write_clustering(sys.argv[3], clustering, silhouette_score)

if clustering_type == '-H':
    print("Clustering using hierarchical method")
    clusterings = cluster_hierarchically(active_sites)
    silhouette_score = return_silhouette_score(clustering, clustering_type)
    write_mult_clusterings(sys.argv[3], clusterings)
