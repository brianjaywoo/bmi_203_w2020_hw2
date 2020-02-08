from hw2skeleton import cluster
from hw2skeleton import io
import pandas as pd
import os

def test_similarity():
    filename_a = os.path.join("data", "276.pdb")
    filename_b = os.path.join("data", "4629.pdb")

    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)

    # update this assertion

    #First test: reflexivity of this distance metric
    assert cluster.compute_similarity(activesite_a, activesite_a) == 0.0

    #Second test: symmetric
    trans_1 = cluster.compute_similarity(activesite_a, activesite_b)
    trans_2 = cluster.compute_similarity(activesite_b, activesite_a)
    assert trans_1 == trans_2

    #Third test: non-negativity
    assert trans_1, trans_2 >= 0



def test_partition_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion

    #Slice out clusters that have an assignment
    clusters_out = cluster.cluster_by_partitioning(active_sites)
    non_empty_clusters = [i for i in clusters_out if len(i) > 1]
    distances_from_centroids = [i[1][1] for i in non_empty_clusters]
    print('non-empty clusters for k-means is', non_empty_clusters)
    print('non-empty cluster distance from centroids are', distances_from_centroids)

    #Check that the length of non-empty cluster list is <= 3 (in case where e.g. two or three residues get assigned to same cluster)
    assert len(non_empty_clusters) <= 3

    #For my particular variant of k-means, check that the distance between the cluster residue
    #and the centroid is >= 0 (nonnegativity)
    for i in distances_from_centroids:
        assert i >= 0

def test_hierarchical_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    clusters_out = cluster.cluster_hierarchically(active_sites)
    print('the clusters out look like', clusters_out)

    #Check length and structure
    hierarchical_structure_check(clusters_out)

    #Let's also check what this hierarchical test looks like for much longer, even-length list:
    pdb_ids = [1806, 3458, 3733, 4629, 6040, 7674, 7780, 8208, 8304, 9776, 10701, 10814]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    clusters_out = cluster.cluster_hierarchically(active_sites)
    print('the clusters out look like', clusters_out)

    #Check length and structure
    hierarchical_structure_check(clusters_out)

def hierarchical_structure_check (clusters):
    #First check that the final list only has 1 cluster (equivalent to length 1):
    assert len(clusters) == 1
    #Check the structure of my list:
    for i in clusters:
        cluster.compute_cut(clusters)
        assert i[0] == 'cluster'
        assert isinstance(i[1], list)
        assert isinstance(i[2], pd.DataFrame)
