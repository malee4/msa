# AUTHOR:Jimmy Ka Ho Chiu & Rick Twee-Hee Ong
# CODE VERSION: 2022, last updated locally June 21, 2022
# LANGUAGE: Python
# SOURCE: "Clustering biological sequences with dynamic sequence similarity threshold"
# URL: https://doi.org/10.1186/s12859-022-04643-9

# The following file is meant to be a modified form of the ALFATClust tool.
# The tool has been altered to allow for direct piping into a test file
# as a preprocessing algorithm for Lindvall's annealing MSA algorithm.

# import the packages
from ClusterEval import *
from Config import *
# from Constants import *
# from Precluster import *
from SeqCluster import *
from SeqSimilarity import *
from Utils import read_seq_file, convert_to_seq_clusters, get_max_precision
from collections import namedtuple

# from math import ceil
# from multiprocessing import Pool
from scipy.sparse import coo_matrix
# import argparse
# import os
import pandas as pd
import sys
import numpy as np



# read in file
seq_file_path = "/Users/melod/Desktop/msa/test_data/argdit_aa_06feb2020_full.fa.txt" # input file path
cluster_eval_csv_file_path = "" # output file path

# read in user parameters

main_dir_path = os.path.dirname(os.path.realpath(__file__))
config_file_path = os.path.join(main_dir_path, 'settings.cfg')

try:
    # this creates an object that will store the settings based on the input file
    config = Config(config_file_path)
except:
    # if the file path is incorrect or the settings.cfg
    sys.exit('Configuration file \'{}\' is not set properly'.format(config_file_path))


# read in file and get info using Utils.read_seq_file
# seq_file_info = read_seq_file(seq_file_path, user_params)
seq_file_info = read_seq_file(seq_file_path)

# get the pairwise sequence similarity
global_edge_weight_mtrx = SeqSimilarity.get_pairwise_similarity(seq_file_info)
print(global_edge_weight_mtrx)
# run it through the clustering algorithm
# seq_cluster_ptrs = SeqCluster.cluster_seqs(global_edge_weight_mtrx)
# sparse_edge_weight_mtrx = coo_matrix(global_edge_weight_mtrx, shape=global_edge_weight_mtrx.shape)

