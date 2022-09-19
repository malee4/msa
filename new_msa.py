from new_Constants import *
from new_cluster import *
from run_lindvall import run_lindvall
from lindvall import get_alignment_string

#from alfatclust import get_clusters_and_centers
def get_sequence_strings(cluster, center_sequence, additional_sequences = list()):
    # if not cluster:
    #     return list()
    sequence_string_set = additional_sequences

    # if cluster is a singleton
    if not cluster:
        return [center_sequence]
    for item in cluster:
        sequence_string_set = sequence_string_set + [str(item.seq)]
    return sequence_string_set + [center_sequence]

def find_spaces(sequence):
    shifts = []

    for i in range(len(sequence)):
        if sequence[i] == "-":
            shifts = shifts + [True]
        else:
            shifts = shifts + [False]  
    return shifts

def shift_spaces(seq_set, shifts): # seq_set is a list of strings
    shifted_set = []

    # manually shift each sequence
    for seq in seq_set:
        shifted_seq = []
        for i in range(len(shifts)):
            # if a "-" exists at the position
            if shifts[i]:
                seq = seq[:i] + ["-"] + seq[i:]
            else:
                continue
        # add the sequence to shifted_set
        shifted_set = shifted_set + [seq]

    return shifted_set

"""
Default format for old_aligned: [aligned_body] + [old_center]
Default format for new_aligned: [old_center] + cluster + [current_center]
"""
def merge_seq_sets(old_aligned, new_aligned):
    # shift old_aligned to match spaces in new_aligned
    old_center_new_aligned = new_aligned[0]
    shifts_old_to_new = find_spaces(old_center_new_aligned)
    old_aligned_shifted_set = shift_spaces(old_aligned, shifts_old_to_new)

    # shift new_aligned to match spaces in old_aligned
    old_center_old_aligned = old_aligned[len(old_aligned) - 1]
    shifts_new_to_old = find_spaces(old_center_old_aligned)
    new_aligned_shifted_set = shift_spaces(new_aligned, shifts_new_to_old)

    return old_aligned_shifted_set + new_aligned_shifted_set[1:] # remove repeat of old_center

if __name__ == '__main__':
    seq_file_path = "/Users/melod/Desktop/msa/test_data/test_file.fa"

    # get clusters
    print("Getting Alignments")
    id_to_center_and_cluster_map = get_clusters_and_centers(seq_file_path)  

    # list of sequences in singletons or small-sized clusters
    additional_sequences = list()

    # declare variables for reformatting
    center_sequence = None
    max_local_cluster_size = -1
    old_center = None
    count = 0

    # empty final output 
    aligned_final = []

    #####################################
    # Default format: [old_center] + cluster + [current_center]
    ####################################
    for id in id_to_center_and_cluster_map:
        count = count + 1
        # initialize the set of strings
        sequence_string_set = list()

        # if this is the first sequence, set the center sequence
        if not center_sequence:
            current_center_sequence = id_to_center_and_cluster_map[id][0]
            center_sequence = id_to_center_and_cluster_map[id][0]
        else:
            current_center_sequence = id_to_center_and_cluster_map[id][0]

        # get the size of the current cluster (list)
        do_not_skip = True # if additional checks are to be skipped
        cluster = id_to_center_and_cluster_map[id][1] 

        if cluster:
            current_cluster_size = len(cluster) + 1 # recall that the center is not a part of the cluster
        else:
            current_cluster_size = 1 
            if count == len(id_to_center_and_cluster_map):
                sequence_string_set = additional_sequences + [str(current_center_sequence.seq)]
                do_not_skip = False
            else:
                # if cluster is singleton, pair it with another cluster
                additional_sequences = additional_sequences + [str(center_sequence.seq)]
                continue
        
        # additional checks
        if do_not_skip:
            if current_cluster_size >= MIN_CLUSTER_PROCESSING_SIZE:
                # get a list of the sequences
                sequence_string_set = get_sequence_strings(cluster, center_sequence)
                # set the center sequence
                center_sequence = current_center_sequence
            elif current_cluster_size + len(additional_sequences) >= MIN_CLUSTER_PROCESSING_SIZE or count == len(id_to_center_and_cluster_map): # not -1 because it starts w/ 1 with the first sequence
                sequence_string_set = get_sequence_strings(cluster, center_sequence, additional_sequences) 
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
                additional_sequences = get_sequence_strings(cluster, center_sequence, additional_sequences)
                continue
        # print("Length of string set: " + str(len(sequence_string_set)))
        print(sequence_string_set)
        results = run_lindvall(sequence_string_set, old_center, simulation=IS_SIMULATION)

        # gets the lowest energy solution, converts to dataframe
        positions = results.lowest().samples()[0]

        # get the alignments
        aligned_strings = get_alignment_string(sequence_string_set, GAPS, positions)
        # print(aligned_strings)
        # if this is the first cluster
        if not old_center:
            aligned_final = aligned_final + aligned_strings
        else:
            # merge aligned strings with past aligned sequences
            aligned_final = merge_seq_sets(aligned_final, aligned_strings)
        
        # for future alignment purposes
        old_center = center_sequence
        center_sequence = None

        # cluster tracking (user output)
        if count % 5 == 0 and count != 0:
            print(str(count) + " clusters aligned...")
    
    for sequence in aligned_final:
        print(sequence)
