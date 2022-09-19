from new_Constants import *
from new_cluster import *
#from alfatclust import get_clusters_and_centers

if __name__ == '__main__':
    seq_file_path = "/Users/melod/Desktop/msa/test_data/test_file.fa"

    # get clusters
    id_to_center_and_cluster_map = get_clusters_and_centers(seq_file_path)  



