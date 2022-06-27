# this runs the algorithm
# FASTA file >> alfatclust

import lindvall as lindvall
from alfatclust import get_clusters_and_centers

"""
Input cluster is tuple
"""
def get_sequence_strings(cluster):
    sequence_string_set = list()
    cluster_center = cluster[0]
    for id in cluster[1]:
        break
    return sequence_string_set, cluster_center

seq_file_path = "/Users/melod/Desktop/msa/test_data/test_file.fa"

id_to_center_and_cluster_map = get_clusters_and_centers(seq_file_path)

# print(id_to_center_and_cluster_map)

for item in id_to_center_and_cluster_map:
    print(type(id_to_center_and_cluster_map[item][1][0]))