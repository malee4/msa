# AUTHOR:Jimmy Ka Ho Chiu & Rick Twee-Hee Ong
# CODE VERSION: 2022, last updated locally June 21, 2022
# LANGUAGE: Python
# SOURCE: "Clustering biological sequences with dynamic sequence similarity threshold"
# URL: https://doi.org/10.1186/s12859-022-04643-9

from Constants import AA
import numpy as np
import os
import re
import shlex
import subprocess
import Mash

class SeqSimilarity:
    _dna_kmer_size = None
    _protein_kmer_size = None
    _dna_sketch_size = None
    _protein_sketch_size = None
    _seed = None
    _min_shared_hash_ratio = None
    _noise_filter_thres = None
    _max_dist = None
    _num_of_threads = None
    _p_value = None
    _is_init = False

    @classmethod
    def init(cls, user_params, p_value=0.0001):
        # if certain parameters are not give, use the default values
        if user_params.kmer_size is None:
            cls._dna_kmer_size = user_params.default_dna_kmer_size
            cls._protein_kmer_size = user_params.default_protein_kmer_size
        else:
            cls._dna_kmer_size = user_params.kmer_size
            cls._protein_kmer_size = user_params.kmer_size

        if user_params.sketch_size is None:
            cls._dna_sketch_size = user_params.default_dna_sketch_size
            cls._protein_sketch_size = user_params.default_protein_sketch_size
        else:
            cls._dna_sketch_size = user_params.sketch_size
            cls._protein_sketch_size = user_params.sketch_size

        cls._seed = user_params.seed
        cls._min_shared_hash_ratio = user_params.min_shared_hash_ratio
        cls._max_dist = 1 - user_params.noise_filter_thres
        cls._num_of_threads = user_params.num_of_threads
        cls._p_value = p_value
        cls._is_init = True

    @classmethod
    def set_to_run_in_single_thread(cls):
        cls._num_of_threads = 1

    # 
    @classmethod
    def _parse_mash_output(cls, fid, mash_seq_name_to_seq_id_map, seq_count):
        max_seq_id = seq_count - 1
        global_edge_weight_mtrx = np.zeros((seq_count, seq_count), dtype=np.float32)
        # print(fid)
        # readFile = os.fdopen(fid)
        # print(readFile.readline())
        # return global_edge_weight_mtrx
        # i = 0
        with os.fdopen(fid) as f:
            while True:
                line = f.readline()
                line_fields = line.rstrip().split('\t')
                seq_name1 = line_fields[0]
                # print(r"line {} reached", i)
                # i = i + 1
                seq_name2 = line_fields[1]
                # print(r"line {} reached", i)
                # i = i + 1
                seq_id1 = mash_seq_name_to_seq_id_map[seq_name1]
                seq_id2 = mash_seq_name_to_seq_id_map[seq_name2]
                # print(r"line {} reached", i)
                # i = i + 1
                # when we have arrived at the end
                if seq_id1 == max_seq_id and seq_id2 == max_seq_id:
                    break
                # print(r"line {} reached", i)
                # i = i + 1
                # skip comparing sequences that refer to the same object
                if seq_id1 == seq_id2:
                    continue
                # print(r"line {} reached", i)
                # i = i + 1
                if cls._min_shared_hash_ratio is not None:
                    m = re.match(r'(\d+)/(\d+)', line_fields[4])
                    if m:
                        if int(m.group(1)) / int(m.group(2)) < cls._min_shared_hash_ratio:
                            continue
                    else:
                        continue
                # print(r"line {} reached", i)
                # i = i + 1
                global_edge_weight_mtrx[seq_id1, seq_id2] = 1 - float(line_fields[2]) # this turns it into a maximization problem
                # print(r"line {} reached", i)
                # i = i + 1
        return global_edge_weight_mtrx

    @classmethod
    def get_pairwise_similarity(cls, seq_file_info):
    
        if not cls._is_init:
            return None

        # format the input
        mash_command = 'mash dist -i -v {} -d {} -p {} {} {}'
        mash_command = mash_command.format(cls._p_value, cls._max_dist, cls._num_of_threads,
                                           seq_file_info.seq_file_path, seq_file_info.seq_file_path)

        if seq_file_info.seq_type == AA:
            mash_command = '{} -a -k {} -s {}'.format(mash_command, cls._protein_kmer_size, cls._protein_sketch_size)
        else:
            mash_command = '{} -k {} -s {}'.format(mash_command, cls._dna_kmer_size, cls._dna_sketch_size)

        if cls._seed is not None:
            mash_command = '{} -S {}'.format(mash_command, cls._seed)
        
        print("Mash commands prepared")
        # print(mash_command)
        # returns read, write definitions
        fr, fw = os.pipe()

        # catches any less common cases (documentation: https://docs.python.org/3/library/subprocess.html#subprocess.Popen)
        args=shlex.split(mash_command)
        with subprocess.Popen(args, stdout=fw, stderr=subprocess.DEVNULL) as p:
            global_edge_weight_mtrx = cls._parse_mash_output(fr, seq_file_info.mash_seq_name_to_seq_id_map,
                                                        seq_file_info.seq_count)
         
        
        os.close(fw)
        # print("Hello")
        # return 0
        return global_edge_weight_mtrx
