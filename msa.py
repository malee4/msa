# this runs the algorithm
# FASTA file >> alfatclust

import modules.lindvall as lindvall
from modules.alfatclust import get_clusters_and_centers


seq_file_path = "/Users/melod/Desktop/msa/test_data/argdit_aa_06feb2020_full.fa"

id_to_center_and_cluster_map = get_clusters_and_centers(seq_file_path)

