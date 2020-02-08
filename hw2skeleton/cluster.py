from .utils import Atom, Residue, ActiveSite
import random
import statistics
import pandas as pd
import numpy as np
import math
import copy
import sklearn
from sklearn import metrics
from sklearn.metrics import pairwise_distances

#This function works on two sequences of exactly the same length, with the structure
#defined in the compute_similarity function below.

#---------------------SECTION START---------------------#


#-----HELPER FUNCTIONS FOR FORMATTING PURPOSES----------#

#This function returns a dataframe from a given active site instance. This DataFrame
#is of the form: rownames being the name of the residues present in the active site;
#first column being the number of that specific amino acid residue present in the
#active site.

def return_dataframe (active_site):

    #String format so I can split
    site = [str(i) for i in active_site.residues]

    #Split into strings, get first element to get only residues
    active_site_names = [i.split()[0] for i in site]
    set_site = set(active_site_names)
    residues_unique = list(set_site)

    #Want to encapsulate each residue in list so I can extend the list with counts
    residues_unique = [[i] for i in residues_unique]

    #Count the occurrences of residues
    for i in residues_unique:
        #Pass unique entries as keys into .count method
        residues_counts = active_site_names.count(str(i[0]))
        #Extend with the counts for the residue
        i.extend([residues_counts])

    #Create dataframe and set index (rownames) to be first column of residues.
    frame = pd.DataFrame(residues_unique)
    frame.columns = ['residue', 'number']

    #Set the index to be the row names and return the frame
    frame = frame.set_index('residue')

    return frame

#This function returns a joined dataframe (either inner or outer join) and optionally
#aggregates the results of the join into one column (takes the mean of the two resulting
#columns). It calls the following procedure to agglomerate the two dataframes.

def return_joined_df (frame_a, frame_b, how = 'outer', agg = False):

    if how == 'outer':
        #Use merge instead of join, because of key problem...
        frame = frame_a.merge(frame_b, on='residue', how='outer')

    elif how == 'inner':
        frame = frame_a.merge(frame_b, on='residue', how='inner')

    if agg:
        frame = agglomerate_df(frame)

    return frame

#This function aggregates two dataframes' columns together by taking means of the two
#columns and replacing the first column with that mean (keeping the structure of the
#input joined dataframe). Used for hierarchical clustering.

def agglomerate_df (joined_frame):

    #Annotate the columns
    joined_frame.columns = ['number', 'to_drop']

    #Calculate the means of the two columns, along the columns (axis = 1)
    joined_mean = joined_frame.mean(axis = 1)

    #Replace the first column with mean value
    joined_frame['number'] = joined_mean

    #Then drop the second column of observations
    joined_frame = joined_frame.drop(columns = ['to_drop'])

    return joined_frame

#This function operates on the level of two dataframes (ie from above procedure),
#and then returns their similarity by my chosen distance metric, described in more
#detail in the homework assignment itself.

def return_custom_similarity (frame_a, frame_b, dissimilarity = 0):

    #Call the above procedure
    joined = return_joined_df (frame_a, frame_b)

    #Iterate through the rows and count instances of mismatches

    for index, row in joined.iterrows():
    #add 1 to this score if either the join leads to a situation in where the residue was only present in one active site and not the other (the isnan check)
        if math.isnan(row[0]) or math.isnan(row[1]):
            dissimilarity += 1
    #or add the difference between the two columns (meaning, the difference in amounts of the amino acid residue present)
        else:
            dissimilarity += abs(row[0] - row[1])

    return dissimilarity

#---------------------SECTION END---------------------#

#---------------------SECTION START---------------------#

#-----FUNCTIONS TO: COMPUTE SIMILARITY, CLUSTER (K-means + Agglomerative 'average' linkage)----------#

#Function to compute similarity between two protein active sites. Returns a single scalar;
#dissimilarity is actually what I'm returning though.

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    # Fill in your code here!
    similarity = 0.0

    #Need to convert to my custom dataframe format for both sites if they aren't formatted already
    if not isinstance(site_a, pd.DataFrame):
        site_a = return_dataframe(site_a)
    if not isinstance(site_b, pd.DataFrame):
        site_b = return_dataframe(site_b)

    #Call above custom similarity calculator
    similarity = return_custom_similarity(site_a, site_b)

    return similarity

#Perform k-means clustering on the (full) set of active sites that I feed into this.
#Works by: initializing centroids; -> iteratively, calculate centroid-active_site distance and
#assign clusters by choosing most similar centroid; -> reestimate centroids as median ->
#terminate if either all cycles are up or the centroids don't change

def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    # Fill in your code here!

    #Parameters for populating centroids.
    protein_abbr = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                    'THR', 'TRP', 'TYR', 'VAL']

    residue_count_max = 10
    cluster_number = 10

    #Initialize centroids: randomly make some dataframes that take the residues present in that active site,
    #as well as random numbers (0-residue_count_max) for the number column.

    centroids = []
    cluster_call = []

    #Populate centroids
    for _ in range(0, cluster_number):
        #Sample some residues to populate
        centroid_residues = random.sample(protein_abbr, residue_count_max)
        #Place into list format for list extension
        centroid_residues = [[i] for i in centroid_residues]

        for i in centroid_residues:
            #Add a random number of 'counts' per residue
            i.extend([random.randint(0, residue_count_max)])

        #Append the centroid in dataframe format to list of centroids
        centroids.append(['cluster' + str(_), pd.DataFrame(centroid_residues)])

        #Also populate an 'empty list' that is just annotated with cluster names
        cluster_call.append(['cluster' + str(_)])

    #Format the centroid dataframe to match the helper functions above
    for i in centroids:
        i[1].columns = ['residue', 'number']
        i[1] = i[1].set_index('residue')

    #Copy my cluster_call object in order to repopulate it on each iteration
    fresh_clusters = copy.deepcopy(cluster_call)

    #Wraps around ($cycles) times at which it will return

    cycles = 20
    for iterations in range (0, cycles):
        print ('Iteration is number', iterations)
        for i in active_sites:
            #Empty distance list in inner loop to keep track of distances between centroid-residues
            distance_list = []
            for j in centroids:
            #Call to compute_similarity to calculate distance between all active sites and all centroid residues
                centroid_cluster = j[0]
                centroid_df = j[1]
                distance_i_j = compute_similarity (i, centroid_df)
                distance_list.append([centroid_cluster, distance_i_j])

            #Take the minimum here: I've calculated dissimilarity above
            highest_similarity = min([i[1] for i in distance_list])

            #Get the cluster the active site 'belongs' to; just the cluster name of the first matching cluster
            centroid_cluster_min = [i for i in distance_list if i[1] == highest_similarity][0]

            #Append it to the list defined above
            for cluster in cluster_call:

                #Check if the name matches
                if centroid_cluster_min[0] == cluster[0]:

                    #Append in list format, the residue (and all its info, as well as the centroid distance)
                    cluster.append([i, centroid_cluster_min[1]])

        print('Updated clusters are', cluster_call)
        #Update step: update the centroids. Choose the residue with median distance from centroids

        #Call the update_centroids function defined below
        new_centroids = update_centroids(cluster_call, [])

        #Parse out old and new centroid information to comapre
        old_names, new_names = [i[0] for i in centroids], [i[0] for i in new_centroids]
        old_dfs, new_dfs = [i[1] for i in centroids], [i[1] for i in new_centroids]

        #tests that the two centroid dataframes are equivalent
        test = zip(old_dfs, new_dfs)

        #A list comparing if the dataframes from the old centroids is the same as the new centroids
        exit = [i[0].equals(i[1]) for i in test]
        if old_names == new_names:
            print('Cluster names unchanged!')
            #Only exit if all clusters are unchanged
            if all(exit):
                print('Stable clusters!')
                return cluster_call

        #Guaranteed termination condition if the cycle number is met
        if iterations == cycles - 1:
            print('Termination!')
            return cluster_call

        #If neither of the above conditions are met: do two things.
        #1. Assign centroids to be new centroids
        #2. Get rid of all of the current assignments
        centroids = new_centroids
        cluster_call = copy.deepcopy(fresh_clusters)

#Function that updates centroids from an old input: this function takes in the current
#cluster_call and then takes the median of the current cluster assignments, assigning one
#of these median points as the new centroid of the cluster.

def update_centroids (cluster_call, centroids):
    #Update step: update the centroids. Choose the residue with median distance from centroids
    for i in cluster_call:
        #Extract relevant labels and residues
        cluster_name = i[0]
        residues = i[1:]

        #Checks for the case in where the cluster is empty, in which case I leave it alone
        if len(residues) > 0:
            #Extract the distances which is in the second entry of each entry in residues
            distance_of_cluster = [distance[1] for distance in residues]

            #Use the median_low function in order to have the datapoint returned be an actual datapoint
            median_distance = int(statistics.median_low(distance_of_cluster))

            #Only consider the first of these 'candidate centroids'
            candidate_centroids = [i[0] for i in residues if i[1] == median_distance]
            new_centroid = return_dataframe(candidate_centroids[0])

            #I pass in the empty list as centroids originally; if the cluster is empty, drop it
            centroids.append([cluster_name, new_centroid])

    return centroids


#This function uses a 'median'-linkage, hierarchical/agglomerative approach to cluster
#all active sites fed in. While all active sites do not belong to a single cluster:
#distance between clusters and all other clusters are computed. The minimal distance-cluster
#is then chosen to lump the current cluster to. Finally, the result of outer-joining the two
#clusters' DataFrame representations is then used as the new representation of the cluster
#going forward.

def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!
    #First: initialize a 'data structure' for how I will represent this.
    active_sites = [['cluster', i, return_dataframe(i)] for i in active_sites]
    counter = 0

    #While loop loops until all points belong to single cluster
    while True:
        counter += 1
        print('Iteration is', counter)

        #Choose a randomly selected element to assign to another cluster
        index = random.randrange(len(active_sites))
        site_of_interest = active_sites[index]
        print('Index of the active site chosen is', index)

        #Computing cluster-cluster distance
        distance_list = []
        for other_site in active_sites:
            df_site = site_of_interest[2] #exists in second entry
            df_other_site = other_site[2]

            #Checks whether the dataframes are equal. Assumes that all starting clusters' dataframes
            #are unique, which could be an inaccurate assumption to make.
            if not df_site.equals(df_other_site):
                distance = compute_similarity(df_site, df_other_site)
                #Append a list of the residue and its distance
                distance_list.append([other_site, distance])

        #Now parse list for the closest cluster
        distances_only = [i[1] for i in distance_list]
        closest_dist = min(distances_only)

        #Get a list of closest neighbors
        closest_neighbors = [i[0] for i in distance_list if i[1] == closest_dist]
        assigned = closest_neighbors[0]

        #Get the index of the to-be assigned active site so we can pop it off later
        index_assigned = active_sites.index(assigned)

        #print('old site looks like', site_of_interest)
        #Notes: 'site_of_interest' is the current active site we're considering. 'assigned' is the closest neighbor we got
        site_of_interest = ['cluster', [site_of_interest, assigned], return_joined_df(df_site, df_other_site, agg = True)]
        #print('new site looks like', site_of_interest)

        #Assign this other active site to the same cluster as site
        active_sites[index] = site_of_interest

        #Pop off the other site from the outer for loop so we don't consider it anymore
        active_sites.pop(index_assigned)
        print('Length of active sites now is', len(active_sites)) #should decrement by 1

        #At the very end, should have one list containing all the active sites.
        if len(active_sites) == 1:
            #Terminate if all points belong to a single cluster.
            print('Final active sites look like', active_sites)
            return active_sites

#---------------------SECTION END---------------------#

#---------------------SECTION START---------------------#

#-----FUNCTIONS TO COMPUTE QUALITY OF CLUSTERING AND COMPARE CLUSTERING METHODS----------#

def return_silhouette_score (clusters, type_of_clustering):
    #Call below procedures to handle this logic to feed into sklearn.metrics.silhouette_score.

    #Extract labels in relevant format for calculating score.
    print('type of clustering is', type_of_clustering)
    labels = extract_labels (clusters, type_of_clustering)

    #Compute distance between all active sites and all other active sites.
    print('The labels are', labels, 'and the length of labels is', len(labels))
    distance_df = compute_distance(clusters, type_of_clustering)

    #Compute a single scalar of silhouette_score.
    silhouette_score = metrics.silhouette_score(distance_df, labels, metric = 'precomputed')
    print('Silhouette score for k_means is', silhouette_score)

    return silhouette_score

#This function extracts relevant labels for feeding into sklearn.metric's silhouette_score calculator.

def extract_labels (clusters, type_of_clustering):
    #Two different ways to extract labels depending on whether k-means or agglomerative clustering chosen
    if type_of_clustering== '-H':
        clusters = bisect_cut(clusters)

    labels = []
    #Clusters here follow [['cluster0', [residues].....]]
    for i in clusters:
        name = i[0]
        #Get length because we just want to label all active sites with their corresponding cluster name
        length_cluster = len(i[1:])
        for i in range(0, length_cluster):
            #Would be implemented as 'rep()' in R
            labels.append(name)

    return labels

#This function is specific to return exactly two clusters from a hierarchical_cluster that arises from
#output from my hierarchical_cluster function. It works by extending two lists with
#the relevant amino acid residues, using nonlocal as a way to bypass messy data structure problems
#that could arise from having to traverse lists (a nested version of which stores the hierarchical_cluster output).

def bisect_cut (hierarchical_cluster):

    #Naively: let's just compare the distance between the left and rightmost clusters.

    left_cluster = []
    right_cluster = []

    #Define a 'tree traversal' helper function within the outer one for use with nonlocal
    def retrieve_clusters (structure, right = False):
        nonlocal left_cluster
        nonlocal right_cluster
        if isinstance(structure[1], ActiveSite):
            if right:
                right_cluster.extend([structure[1]])
            else:
                left_cluster.extend([structure[1]])
            pass
        else:
            left = structure[1][0]
            right = structure[1][1]
            return [retrieve_clusters(left), retrieve_clusters(right)]

    #The top level is a list of length 1
    top_level = hierarchical_cluster
    tree = retrieve_clusters(top_level)

    #The correct structure to feed into compute_distance below
    cluster_list = [['cluster0', left_cluster], ['cluster1', right_cluster]]
    print('the clusters look like', cluster_list)

    return cluster_list

#This function computes distance between all clusters fed in and returns an
#n x n data matrix with distances (so 0's along the diagonals).

def compute_distance (clusters, type_of_clustering):
    residues_only = []
    if type_of_clustering== '-H':
        clusters = bisect_cut(clusters)
        for i in clusters:
        #Just get residues for purposes of calculating distances
            residues = i[1]
            residues_only.extend(residues)
    else:
        for i in clusters:
        #Just get residues for purposes of calculating distances
            print('i is like', i)
            residues = i[1:]
            residues = [i[0] for i in residues] #only get the AA residues
            residues_only.extend(residues)

    #Calculate distance and convert the list of distance lists, to a pandas DataFrame
    convert_to_df = []
    for i in residues_only:
        add_to_df = []
        for j in residues_only:
            distance = compute_similarity(i, j)
            add_to_df.append(distance)
        #Second append statement to append list of distances for a particular residue, to the 'outer list'
        convert_to_df.append(add_to_df)

    #Convert the list of distance lists, to a pd.DataFrame
    distance_df = pd.DataFrame(convert_to_df)

    #0's on the diagonal
    print('The final dataframe looks like', distance_df)
    return distance_df

def compare_clusterings ():
    #This is not implemented due to inability on my end. However, I describe what would have been done here in the
    #actual assignment pdf.
    pass

#---------------------SECTION END---------------------#
