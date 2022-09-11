# AUTHOR:Jimmy Ka Ho Chiu & Rick Twee-Hee Ong
# CODE VERSION: 2022, last updated locally June 21, 2022
# LANGUAGE: Python
# SOURCE: "Clustering biological sequences with dynamic sequence similarity threshold"
# URL: https://doi.org/10.1186/s12859-022-04643-9

from new_Constants import AA

class SeqSimilarity:
    # default parameters prior to initialization
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

    # @classmethod 
    # def init(cls, config, p_value=0.0001):
    #     if config['kmer_size'] is None:
    #         cls._dna_kmer_size = config['default_dna_kmer_size']
    #         cls._protein_kmer_size = config['default_protein_kmer_size']
    #     else:
    #         cls._dna_kmer_size = config['kmer_size'] 
    #         cls._protein_kmer_size = config['kmer_size'] 

    #     # MAJOR CHANGE: no check default/general sketch size, directly assign dna and protein sketch sizes
    #     cls._dna_sketch_size = config['default_dna_sketch_size']
    #     cls._protein_sketch_size = config['default_protein_sketch_size']

        
    #     cls._seed = user_params.seed
    #     cls._min_shared_hash_ratio = user_params.min_shared_hash_ratio
    #     cls._max_dist = 1 - user_params.noise_filter_thres
    #     cls._num_of_threads = user_params.num_of_threads
    #     cls._p_value = p_value
    #     cls._is_init = True
