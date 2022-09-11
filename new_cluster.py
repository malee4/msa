from shutil import rmtree
# from alfatclust import internal_parse_to_user_params
from collections import namedtuple
from new_Constants import *
from SeqCluster import *
from SeqSimilarity import *
import os
import sys

def internal_parse_to_user_params(seq_file_path, config):
    param_error_log = list() 

    UserParams = namedtuple('UserParams', ['res_param_start', 'res_param_end', 'res_param_step_size', 'precision',
                                           'precluster_thres', 'min_shared_hash_ratio', 'kmer_size',
                                           'default_dna_kmer_size', 'default_protein_kmer_size', 'sketch_size',
                                           'default_dna_sketch_size', 'default_protein_sketch_size',
                                           'noise_filter_thres', 'num_of_threads', 'seed'])
    
def get_clusters_and_centers(seq_file_path, is_precluster_mode = False):
    cluster_ids_to_centers_and_cluster_seqs = dict()
    # file locations
    main_dir_path = os.path.dirname(os.path.realpath(__file__))

    # SET UP CONFIGURATIONS
    config = getConfig()

    try:
        # TODO: run function to check configurations
        # if config not valid us sys.exit()
        user_params, param_error_log = internal_parse_to_user_params(seq_file_path, config)

        # TODO: print out user parameters

        # initialize SeqClust and SeqSimilarity
        SeqSimilarity.init(user_params)
        SeqCluster.init(user_params)

        print()

    except KeyboardInterrupt:
        print()
        print('Process aborted due to keyboard interrupt')
    except SystemExit as sys_exit:
        if sys_exit.code != 0:
            print(sys_exit.code)
    except:
        print()
        print('Process aborted due to error occurred: {}'.format(sys.exc_info()[1]))
    # finally:
    #     Precluster.clear_temp_data()