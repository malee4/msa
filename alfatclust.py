# AUTHOR:Jimmy Ka Ho Chiu & Rick Twee-Hee Ong
# CODE VERSION: 2022, last updated locally June 21, 2022
# LANGUAGE: Python
# SOURCE: "Clustering biological sequences with dynamic sequence similarity threshold"
# URL: https://doi.org/10.1186/s12859-022-04643-9

# The following file is meant to be a modified form of the ALFATClust tool.
# The tool has been altered to allow for direct piping into a test file
# as a preprocessing algorithm for Lindvall's annealing MSA algorithm.

# import packages
from ClusterEval import *
from Config import *
from Constants import *
from Precluster import *
from SeqCluster import *
from SeqSimilarity import *
from Utils import read_seq_file, convert_to_seq_clusters, get_max_precision
from collections import namedtuple
from math import ceil
from multiprocessing import Pool
from scipy.sparse import coo_matrix
import argparse
import os
import pandas as pd
import sys

# file locations
config_file_path = "settings.cfg"

def set_and_parse_args(config):
    num_of_threads = len(os.cpu_count())

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='seq_file_path', required=True,
                        help='DNA/Protein sequence FASTA file path')
    parser.add_argument('-o', '--output', action='store', dest='seq_cluster_file_path', required=True,
                        help='Output cluster file path')
    parser.add_argument('-e', '--evaluate', action='store', dest='cluster_eval_csv_file_path',
                        help='Sequence cluster evaluation output CSV file path')
    help_msg = 'Lower bound for estimated similarity range [{}]'
    parser.add_argument('-l', '--low', type=float, help=help_msg.format(config.res_param_end),
                        default=config.res_param_end)
    help_msg = 'Estimated similarity step size [{}]'
    parser.add_argument('-d', '--step', type=float, help=help_msg.format(config.res_param_step_size),
                        default=config.res_param_step_size)
    parser.add_argument('-f', '--filter', type=float, help='Pairwise shared hash ratio threshold for filtering')
    parser.add_argument('-p', '--precluster', action='store_true', dest='is_precluster',
                        help='Always pre-cluster sequences into individual subsets for separate clustering')
    parser.add_argument('-k', '--kmer', type=int, help='K-mer size for Mash')
    parser.add_argument('-s', '--sketch', type=int, help='Sketch size for Mash')
    help_msg = 'Discard estimated similarity values below (lower bound of estimated similarity - margin [{}])'
    parser.add_argument('-m', '--margin', type=float, help=help_msg.format(config.noise_filter_margin),
                        default=config.noise_filter_margin)
    help_msg = 'No. of threads to be used [{}]'
    parser.add_argument('-t', '--thread', type=int, help=help_msg.format(num_of_threads), default=num_of_threads)
    parser.add_argument('-S', '--seed', type=int, help='Seed value')

    return parser.parse_args()

# takes in the parameters from settings.cfg, makes sure params are valid values, returns user params as a named tuple
def parse_to_user_params(args, config):
    param_error_log = list() # this will be returned alongside any other outputs

    UserParams = namedtuple('UserParams', ['res_param_start', 'res_param_end', 'res_param_step_size', 'precision',
                                           'precluster_thres', 'min_shared_hash_ratio', 'kmer_size',
                                           'default_dna_kmer_size', 'default_protein_kmer_size', 'sketch_size',
                                           'default_dna_sketch_size', 'default_protein_sketch_size',
                                           'noise_filter_thres', 'num_of_threads', 'seed'])

    if not os.path.isfile(args.seq_file_path):
        param_error_log.append('Sequence file \'{}\' does not exist'.format(args.seq_file_path))

    # makes sure the configurations are valid
    if config.res_param_start > 1 or config.res_param_start <= 0:
        param_error_log.append('Upper bound for estimated similarity range must be > 0 and <= 1')
        param_error_log.append('Please check the configuration file')

    if config.precluster_thres <= 0:
        param_error_log.append('Precluster threshold must be positive integer')
        param_error_log.append('Please check the configuration file')

    if args.low > 1 or args.low <= 0:
        param_error_log.append('Lower bound for estimated similarity range must be > 0 and <= 1')

    if args.low >= config.res_param_start:
        param_error_log.append('Lower bound for estimated similarity range must be smaller than upper bound')

    if args.step > 1 or args.step <= 0:
        param_error_log.append('Estimated similarity step size must be > 0 and <= 1')

    if args.filter is not None and (args.filter >= 1 or args.filter <= 0):
        param_error_log.append('Pairwise shared hash ratio threshold must be > 0 and < 1')

    if args.kmer is None:
        if config.default_dna_kmer_size <= 0:
            param_error_log.append('Default k-mer size for DNA sequences must be a positive integer')
            param_error_log.append('Please check the configuration file')

        if config.default_protein_kmer_size <= 0:
            param_error_log.append('Default k-mer size for protein sequences must be a positive integer')
            param_error_log.append('Please check the configuration file')
    elif args.kmer <= 0:
        param_error_log.append('K-mer size must be a positive integer')

    if args.sketch is None:
        if config.default_dna_sketch_size <= 0:
            param_error_log.append('Default sketch size for DNA sequences must be a positive integer')
            param_error_log.append('Please check the configuration file')

        if config.default_protein_sketch_size <= 0:
            param_error_log.append('Default sketch size for protein sequences must be a positive integer')
            param_error_log.append('Please check the configuration file')
    elif args.sketch <= 0:
        param_error_log.append('Sketch size must be a positive integer')

    if args.margin >= 1 or args.margin < 0:
        param_error_log.append('Estimated similarity filtering margin must be >= 0 and < 1')

    if args.thread <= 0:
        param_error_log.append('No. of threads must be a positive integer')

    if args.seed is not None and args.seed < 0:
        param_error_log.append('Seed value must be a non-negative integer')

    if len(param_error_log) > 0:
        return None, param_error_log

    noise_filter_thres = round(max(0, args.low - args.margin), get_max_precision(args.low, args.margin))

    return UserParams(config.res_param_start, args.low, -1 * args.step, get_max_precision(args.step),
                      config.precluster_thres, args.filter, args.kmer, config.default_dna_kmer_size,
                      config.default_protein_kmer_size, args.sketch, config.default_dna_sketch_size,
                      config.default_protein_sketch_size, noise_filter_thres, args.thread, args.seed), None

def display_user_params(user_params):
    print('---------------------------------------------')
    print('Estimated similarity range = [{}, {}]'.format(user_params.res_param_start, user_params.res_param_end))
    print('Estimated similarity step size = {}'.format(user_params.res_param_step_size * -1))

    if user_params.min_shared_hash_ratio is not None:
        print('Pairwise min. shared hash ratio = {}'.format(user_params.min_shared_hash_ratio))

    if user_params.kmer_size is None:
        print('Default DNA k-mer size = {}'.format(user_params.default_dna_kmer_size))
        print('Default protein k-mer size = {}'.format(user_params.default_protein_kmer_size))
    else:
        print('K-mer size = {}'.format(user_params.kmer_size))

    if user_params.sketch_size is None:
        print('Default DNA sketch size = {}'.format(user_params.default_dna_sketch_size))
        print('Default protein sketch size = {}'.format(user_params.default_protein_sketch_size))
    else:
        print('Sketch size = {}'.format(user_params.sketch_size))

    print('Min. estimated similarity considered = {}'.format(user_params.noise_filter_thres))
    print('No. of threads = {}'.format(user_params.num_of_threads))
    print('---------------------------------------------')
    print()

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

## NOTE TO REMOVE?
if __name__ == '__main__':
    main_dir_path = os.path.dirname(os.path.realpath(__file__))
    config_file_path = os.path.join(main_dir_path, 'settings.cfg')

    try:
        # this creates an object that will store the settings based on the input file
        config = Config(config_file_path)
    except:
        # if the file path is incorrect or the settings.cfg
        sys.exit('Configuration file \'{}\' is not set properly'.format(config_file_path))

    temp_seq_file_path = None

    try:
        # based on the input 
        args = set_and_parse_args(config)
        # check if user params are valid; if yes, return as named tuple
        user_params, param_error_log = parse_to_user_params(args, config)
        # if an error is encountered, exit
        if param_error_log is not None:
            sys.exit(os.linesep.join(param_error_log))

        # else show user the input parameters
        display_user_params(user_params)

        SeqSimilarity.init(user_params)
        SeqCluster.init(user_params)

        print('Validating input sequence file \'{}\'...'.format(args.seq_file_path))
        seq_file_info = read_seq_file(args.seq_file_path, user_params)
        if seq_file_info.seq_count == 0:
            sys.exit('No sequence found in \'{}\''.format(args.seq_file_path))

        if len(seq_file_info.error_log) > 0:
            sys.exit(os.linesep.join(seq_file_info.error_log))

        cluster_eval_output_df = None
        output_seq_clusters = list()
        overall_error_log = list()
        is_precluster_mode = args.is_precluster or seq_file_info.seq_count > user_params.precluster_thres

        # from the --precluster help function: "Always pre-cluster sequences into individual subsets for separate clustering"
        if is_precluster_mode:
            print('Pre-clustering sequences into subsets...')
            # 
            precluster_to_seq_recs_map, max_precluster_size = \
                Precluster.precluster_seq_file(user_params, args.seq_file_path, seq_file_info.max_seq_len)
            num_of_preclusters = len(precluster_to_seq_recs_map)
            print('{} individual subsets to be clustered'.format(num_of_preclusters))

            if max_precluster_size < user_params.precluster_thres:
                SeqSimilarity.set_to_run_in_single_thread()
                num_of_threads_for_main_loop = user_params.num_of_threads
            else:
                num_of_threads_for_main_loop = 1

            Precluster.create_temp_dir()
            SeqCluster.disable_verbose()

            if args.cluster_eval_csv_file_path is not None:
                ClusterEval.init(config, 1)

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
            print('Estimating pairwise sequence distances...')
            global_edge_weight_mtrx = SeqSimilarity.get_pairwise_similarity(seq_file_info)

            if args.cluster_eval_csv_file_path is not None:
                sparse_edge_weight_mtrx = coo_matrix(global_edge_weight_mtrx, shape=global_edge_weight_mtrx.shape)

            seq_cluster_ptrs = SeqCluster.cluster_seqs(global_edge_weight_mtrx)
            output_seq_clusters = convert_to_seq_clusters(seq_cluster_ptrs, seq_file_info.seq_id_to_seq_name_map)

            if args.cluster_eval_csv_file_path is not None:
                print('Evaluating sequence clusters...')
                ClusterEval.init(config, user_params.num_of_threads)
                cluster_eval_output_df = \
                    ClusterEval.eval_clusters(seq_cluster_ptrs, sparse_edge_weight_mtrx.toarray(), seq_file_info)

        with open(args.seq_cluster_file_path, 'w') as f_out:
            cluster_count = 0

            for seq_cluster in output_seq_clusters:
                cluster_count += 1
                f_out.write('#Cluster {}{}'.format(cluster_count, os.linesep))
                f_out.writelines(seq_cluster)

        if cluster_eval_output_df is not None:
            cluster_eval_output_df.to_csv(args.cluster_eval_csv_file_path)

        print()
        print('Process completed. No. of sequence clusters = {}'.format(cluster_count))
    except KeyboardInterrupt:
        print()
        print('Process aborted due to keyboard interrupt')
    except SystemExit as sys_exit:
        if sys_exit.code != 0:
            print(sys_exit.code)
    except:
        print()
        print('Process aborted due to error occurred: {}'.format(sys.exc_info()[1]))
    finally:
        Precluster.clear_temp_data()