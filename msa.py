# this runs the algorithm
# FASTA file >> alfatclust

from lindvall import get_alignment_string
from run_lindvall import run_lindvall
from alfatclust import get_clusters_and_centers

"""
Input cluster is tuple
"""
def get_sequence_strings(cluster, additional_sequences = list()):
    # if not cluster:
    #     return list()
    sequence_string_set = additional_sequences
    
    for item in cluster:
        sequence_string_set = sequence_string_set + [str(item.seq)]
    return sequence_string_set


if __name__ == '__main__':
    # turn this into an input file
    seq_file_path = "/Users/melod/Desktop/msa/test_data/test_file.fa"

    # settings
    min_cluster_processing_size = 2
    gaps = 0

    # run the clustering algorithm to get a dictionary of cluster and their centers
    id_to_center_and_cluster_map = get_clusters_and_centers(seq_file_path)  
    
    print("---------------------------------------------")
    print()
    print("Getting sequence strings")

    # list of sequences in singletons or small-sized clusters
    additional_sequences = list()

    center_sequence = None
    max_local_cluster_size = -1
    old_center = None
    count = 0

    # empty final output 
    aligned_final = [[]]

    print("Getting alignments")
    print()
    for id in id_to_center_and_cluster_map:
        count = count + 1
        # initialize the set of strings
        print("Sequence string set initialized")
        sequence_string_set = list()

        # if this is the first sequence, set the center sequence
        if not center_sequence:
            current_center_sequence = id_to_center_and_cluster_map[id][0]
            center_sequence = id_to_center_and_cluster_map[id][0]
        else:
            current_center_sequence = id_to_center_and_cluster_map[id][0]

        # get the size of the current cluster (list)
        cluster = id_to_center_and_cluster_map[id][1]
        current_cluster_size = len(cluster)
        print("cluster type " + str(type(cluster)))

        if current_cluster_size >= min_cluster_processing_size:
            # get a list of the sequences
            sequence_string_set = get_sequence_strings(cluster)
            # set the center sequence
            center_sequence = current_center_sequence
        elif current_cluster_size + len(additional_sequences) >= min_cluster_processing_size or count == len(id_to_center_and_cluster_map):
            sequence_string_set = get_sequence_strings(cluster, additional_sequences)
            if current_cluster_size > max_local_cluster_size:
                center_sequence = current_center_sequence
            # reset setting for next iteration
            additional_sequences.clear()
            max_local_cluster_size = -1
        else:
            # find the largest cluster size. The center will be that of the largest cluster. 
            # SHORT COMING: favors clusters that were prioritized earlier, does not account for how closely related clusters are
            if current_cluster_size > max_local_cluster_size:
                center_sequence = current_center_sequence
                max_local_cluster_size = current_cluster_size
            # if threshold is not met, append additional sequences and continue onto next cluster
            additional_sequences = get_sequence_strings(cluster, additional_sequences)
            continue
        # print(sequence_string_set)
        print("Running Lindvall algorithm...")
        results = run_lindvall(sequence_string_set, old_center)
        print("Results retrieved")
        # gets the lowest energy solution, converts to dataframe
        positions = results.lowest().to_pandas_dataframe()
        # NOTE: ASSUMES IS NOT SIMULATION, removes last three columns
        positions = positions[:(len(positions)-3)]

        # convert process to strings
        align_strings = get_alignment_string(sequence_string_set, gaps, positions)
        


        # return the alignment
        # if this is not the first cluster
        # if old_center:
            
        #     print()

        # somehow process so that it can be added to some existing database
        
        # for future alignment purposes
        old_center = center_sequence
        center_sequence = None

        # cluster tracking
        if count % 100 == 0 and count != 0:
            print(str(count) + " clusters aligned...")

        # and then run it again

