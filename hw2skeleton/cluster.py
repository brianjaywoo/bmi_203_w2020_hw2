from .utils import Atom, Residue, ActiveSite
import random
import statistics

def return_shorter (seq_a, seq_b):
    if len(seq_a) <= len(seq_b):
        return seq_a
    return seq_b

def return_longer (seq_a, seq_b):
    if len(seq_a) >= len(seq_b):
        return seq_a
    return seq_b

def calculate_distance (seq):
    seq = [int(i) for i in seq]
    new = []
    for i in range(0, len(seq) - 1):
        distance = seq[i + 1] - seq[i]
        new.append(distance)
    return new

#This function works on two sequences of exactly the same length, with the structure
#defined in the compute_similarity function below.
def return_custom_similarity (seq_a, seq_b):
    #Get residue names and also inter-residue distance
    seq_a_residues, seq_b_residues = [i[0] for i in seq_a], [i[0] for i in seq_b]
    seq_a_prim, seq_b_prim = [i[1] for i in seq_a], [i[1] for i in seq_b]
    distance_a, distance_b = calculate_distance(seq_a_prim), calculate_distance(seq_b_prim)

    #Linearly add up the number of matches, first by residue, then by inter-residue distance
    counter = 0
    residues, distances = zip (seq_a_residues, seq_b_residues), zip(distance_a, distance_b)
    for i in residues:
        if i[0] == i[1]:
            counter += 1
    for i in distances:
        if i[0] == i[1]:
            counter += 1
    return counter

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """

    similarity = 0.0
    residues_a = [str(i) for i in site_a.residues]
    residues_b = [str(i) for i in site_b.residues]

    #Split into strings
    residues_a_names = [i.split() for i in residues_a]
    residues_b_names = [i.split() for i in residues_b]

    #residues_a_prim_seq = [i.split()[1] for i in residues_a]
    #residues_b_prim_seq = [i.split()[1] for i in residues_b]

    #Finding a seed
    shorter_seq = return_shorter(residues_a_names, residues_b_names)
    longer_seq = return_longer(residues_a_names, residues_b_names)
    seed = shorter_seq[0][0]
    score_chart = []
    for i in range(0, len(longer_seq) - len(shorter_seq) + 1):
        if seed == longer_seq[i][0]: #means the amino acid matches
            part_of_longer_seq = longer_seq[i:len(shorter_seq)]
            score = return_custom_similarity(shorter_seq, part_of_longer_seq)
            score_chart.append(score)

    # Fill in your code here!
    # Return similarity score that is maximal
    if len(score_chart) > 0:
        similarity = max(score_chart)
        #print(similarity)
        return similarity
    return 0

def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    # Fill in your code here!
    #Initialize centroids and have a
    centroids = []
    blank_cluster = []
    for i in range(0,10): #k-means with k of 10
        n = random.choice(active_sites)
        #n = random.randint(0, max(distances)) #max number possible is max possible distance
        centroids.append(['cluster ' + str(i), n]) #so 'n' is the actual centroid active site
        blank_cluster.append(['cluster ' + str(i)])

    #Call to compute_similarity between active sites and the cluster centroids
    for _ in range (0, 50):
        print(_)
        cluster_call = blank_cluster.copy()
    #if cluster_call == :
        for i in active_sites:
            distance_list = []
            for j in centroids: #Calculate distance between residue in active_sites and all centroids
                centroid_cluster = j[0]
                centroid_residue = j[1]
                name_i, name_j = i.name, centroid_residue.name
                if name_i != name_j: #Don't calculate distance between something & itself
                    distance_i_j = compute_similarity (i, centroid_residue)
                    distance_list.append([centroid_cluster, distance_i_j])
            highest_similarity = max([i[1] for i in distance_list])
            #Get the cluster the active site 'belongs' to; just the cluster name of the first matching cluster
            centroid_cluster = [i for i in distance_list if i[1] == highest_similarity][0]
            #Append it to the list defined above
            for cluster in cluster_call:
                if centroid_cluster[0] == cluster[0]:
                    cluster.append([i, centroid_cluster[1]]) #this was the max distance that got it assigned
        #Update step: update the centroids. Choose the residue with median distance from centroids
        centroids = []
        for i in cluster_call:
            cluster_call_name = i[0]
            if len(i) > 1:
                residues = i[1:]
                median_distance = int(statistics.median([i[1] for i in residues]))
                print(median_distance)
                new_centroid = [i[0] for i in residues if i[1] == median_distance][0]
                print(new_centroid.residues) #just the first element is fine
                centroids.append([cluster_call_name, new_centroid])
            else:
                n = random.choice(active_sites)
                centroids.append([cluster_call_name, n])

    #print(cluster_call)



    return []


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []
