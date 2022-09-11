from shutil import rmtree
# from alfatclust import internal_parse_to_user_params
from collections import namedtuple
from new_Constants import *
from new_SeqCluster import *
from new_SeqSimilarity import *
from new_Utils import get_max_precision
import os
import sys

# MOFIFIED FROM ALFATCLUST FILE
def internal_parse_to_user_params(seq_file_path, config):
    param_error_log = list() 

    UserParams = namedtuple('UserParams', ['res_param_start', 'res_param_end', 'res_param_step_size', 'precision',
                                           'precluster_thres', 'min_shared_hash_ratio', 'kmer_size',
                                           'default_dna_kmer_size', 'default_protein_kmer_size', 'sketch_size',
                                           'default_dna_sketch_size', 'default_protein_sketch_size',
                                           'noise_filter_thres', 'num_of_threads', 'seed'])
    
    # make sure any configurations are valid
    if not os.path.isfile(seq_file_path):
        param_error_log.append('Sequence file \'{}\' does not exist'.format(seq_file_path))

    if config['estimated_similarity']['high'] > 1 or config['estimated_similarity']['high'] <= 0:
        param_error_log.append('Upper bound for estimated similarity range must be > 0 and <= 1')
        param_error_log.append('Please check the configuration file')

    if config['threshold_precluster'] <= 0:
        param_error_log.append('Precluster threshold must be positive integer')
        param_error_log.append('Please check the configuration file')

    # if any errors have been identified, return None
    if len(param_error_log) > 0:
        return None, param_error_log

    # replace any args with config defaults, assume filter margin =  0.2, args.kmer = 17, sketch = 2000, seed = 0
    noise_filter_thres = round(max(0, config['estimated_similarity']['low'] - config['noise_filter_margin']), get_max_precision(config['estimated_similarity']['low'], config['noise_filter_margin']))


    return UserParams(config['estimated_similarity']['high'], config['estimated_similarity']['low'], -1 * config['estimated_similarity']['stepsize'], get_max_precision(config.res_param_step_size),
                     config['threshold_precluster'],  None, None, config['dna_mash']['kmer'],
                      config['protein_mash']['kmer'], None, config['dna_mash']['sketch'],
                      config['protein_mash']['sketch'], noise_filter_thres, os.cpu_count(), None), None 

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
        # STOPPED HERE
        new_SeqSimilarity.init(user_params)
        SeqSimilarity.init(user_params)
        # new_SeqCluster.init(user_params)


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