# AUTHOR:Jimmy Ka Ho Chiu & Rick Twee-Hee Ong
# CODE VERSION: 2022, last updated locally June 21, 2022
# LANGUAGE: Python
# SOURCE: "Clustering biological sequences with dynamic sequence similarity threshold"
# URL: https://doi.org/10.1186/s12859-022-04643-9

# The following file is meant to be a modified form of the ALFATClust tool.
# The tool has been altered to allow for direct piping into a test file
# as a preprocessing algorithm for Lindvall's annealing MSA algorithm.

# import packages
# from ClusterEval import *
# from Config import *
# from Constants import *
# from Precluster import *
# from SeqCluster import *
# from SeqSimilarity import *
# from Utils import read_seq_file, convert_to_seq_clusters, get_max_precision
# from collections import namedtuple
# from math import ceil
# from multiprocessing import Pool
# from scipy.sparse import coo_matrix
# import argparse
# import os
# import pandas as pd
# import sys


# configure settings


# read in file
sequence_file_path = "/Users/melod/Desktop/msa/test_data/argdit_aa_06feb2020_full.fa.txt" # input file path
cluster_file_path = "" # output file path

# stores sequences in a dictionary
sequences = {}

with open(sequence_file_path) as f:
    print("hi")
    segments = f.read().split(">")
    for item in segments:
        if len(item) == 0:
            continue
        key = item[:item.index("\n")]
        value = item[item.index("\n"):].replace("\n", "")

