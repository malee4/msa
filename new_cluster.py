from shutil import rmtree
# from alfatclust import internal_parse_to_user_params
from collections import namedtuple
from new_Constants import *
from SeqCluster import *
from SeqSimilarity import *
from Precluster import *
from new_Utils import get_max_precision, read_seq_file, convert_to_seq_clusters
import os
import sys
from ClusterEval import *

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
    

    return UserParams(config['estimated_similarity']['high'], config['estimated_similarity']['low'], -1 * config['estimated_similarity']['stepsize'], get_max_precision(config['estimated_similarity']['stepsize']),
                     config['threshold_precluster'],  None, None, config['dna_mash']['kmer'],
                      config['protein_mash']['kmer'], None, config['dna_mash']['sketch'],
                      config['protein_mash']['sketch'], noise_filter_thres, os.cpu_count(), None), None 

def cluster_seqs_in_precluster(precluster_seq_records):
    if len(precluster_seq_records) == 1:
        return [['{}{}'.format(precluster_seq_records[0].description, os.linesep)]], None, list()

    temp_seq_file_path = Precluster.write_precluster_seq_records(precluster_seq_records)
    seq_file_info = read_seq_file(temp_seq_file_path)

    if len(seq_file_info.error_log) > 0:
        os.remove(temp_seq_file_path)

        return list(), None, seq_file_info.error_log

    global_edge_weight_mtrx = SeqSimilarity.get_pairwise_similarity(seq_file_info)
    sparse_edge_weight_mtrx = coo_matrix(global_edge_weight_mtrx, shape=global_edge_weight_mtrx.shape)

    seq_cluster_ptrs = SeqCluster.cluster_seqs(global_edge_weight_mtrx)
    cluster_eval_output_df = \
        ClusterEval.eval_clusters_single_thread(seq_cluster_ptrs, sparse_edge_weight_mtrx.toarray(), seq_file_info)

    os.remove(temp_seq_file_path)

    return convert_to_seq_clusters(seq_cluster_ptrs, seq_file_info.seq_id_to_seq_name_map), cluster_eval_output_df, \
        list()

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
        
        # validating input sequence file
        seq_file_info = read_seq_file(seq_file_path, user_params)
        
        # if no sequences have been read in
        if seq_file_info.seq_count == 0:
            sys.exit('No sequence found in \'{}\''.format(seq_file_path))

        # if an error is encountered
        if len(seq_file_info.error_log) > 0:
            sys.exit(os.linesep.join(seq_file_info.error_log))
        
        cluster_eval_output_df = None
        output_seq_clusters = list()
        overall_error_log = list() # to keep track of errors encountered
        is_precluster_mode = seq_file_info.seq_count > user_params.precluster_thres 
        
        # TODO: precluster mode functionality
        if is_precluster_mode:
            print('Pre-clustering sequences into subsets...')
            precluster_to_seq_recs_map, max_precluster_size = \
                Precluster.precluster_seq_file(user_params, seq_file_path, seq_file_info.max_seq_len)
            num_of_preclusters = len(precluster_to_seq_recs_map)
            print('{} individual subsets to be clustered'.format(num_of_preclusters))

            if max_precluster_size < user_params.precluster_thres:
                SeqSimilarity.set_to_run_in_single_thread()
                num_of_threads_for_main_loop = user_params.num_of_threads
            else:
                num_of_threads_for_main_loop = 1

            Precluster.create_temp_dir()
            SeqCluster.disable_verbose()

            chunk_size = min(ceil(num_of_preclusters / num_of_threads_for_main_loop), 500)
            last_max_cluster_id = 0
            process_count = 0

            with Pool(processes=num_of_threads_for_main_loop, maxtasksperchild=40) as pool:
                for seq_clusters, block_cluster_eval_output_df, error_log in \
                    pool.imap_unordered(cluster_seqs_in_precluster, precluster_to_seq_recs_map.values(), chunk_size):
                    output_seq_clusters += seq_clusters
                    overall_error_log += error_log

                    if block_cluster_eval_output_df is not None:
                        block_cluster_eval_output_df.index += last_max_cluster_id

                        if cluster_eval_output_df is None:
                            cluster_eval_output_df = block_cluster_eval_output_df
                        else:
                            cluster_eval_output_df = pd.concat([cluster_eval_output_df, block_cluster_eval_output_df])

                    last_max_cluster_id = len(output_seq_clusters)
                    process_count += 1
                    print('{} / {} subsets processed'.format(process_count, num_of_preclusters), end='\r')

        else:
            
            print("Estimating pairwise sequence distances") # for progress tracking purposes
            global_edge_weight_mtrx = SeqSimilarity.get_pairwise_similarity(seq_file_info)
            
            print("Finding raw sequence cluster")
            seq_cluster_ptrs = SeqCluster.cluster_seqs(global_edge_weight_mtrx)

            # KEY DIFFERENCE: find the centers first

            # get the centers
            print("Getting cluster centers...")
            cluster_ids_to_centers_and_cluster_seqs, count = ClusterEval.get_centers(seq_cluster_ptrs, global_edge_weight_mtrx, seq_file_info.seq_file_path)
            
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
    
    return cluster_ids_to_centers_and_cluster_seqs