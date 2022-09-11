from Bio import SeqIO
from new_Constants import *
from collections import namedtuple
import re

def get_max_precision(*vals):
    max_precision = 2
    
    for val in vals:
        m = re.match(r'\d+\.(\d*[1-9])0*', str(val))
        if m:
            precision = len(m.group(1))
            if precision > max_precision:
                max_precision = precision

    return max_precision

def check_seq_type(seq_str):
    if re.match(IUPAC_DNA_STR_PATTERN, seq_str):
        unambig_dna_seq_str = re.sub(IUPAC_AMBIG_DNA_BASES, '', seq_str)
        if len(unambig_dna_seq_str) / len(seq_str) >= 0.95:
            return DNA
        else:
            return AA # amino acids
    elif re.match(IUPAC_AA_STR_PATTERN, seq_str):
        unambig_protein_seq_str = re.sub(IUPAC_AMBIG_AA_BASES, '', seq_str)
        if len(unambig_protein_seq_str) / len(seq_str) >= 0.95:
            return AA
    return None

def read_seq_file(seq_file_path, user_params):
    SeqFileInfo = namedtuple('SeqFileInfo', ['mash_seq_name_to_seq_id_map', 'seq_id_to_seq_name_map',
                                             'error_log', 'seq_file_path', 'seq_type', 'mash_last_seq_name',
                                             'seq_count', 'max_seq_len'])
    mash_seq_name_to_seq_id_map = dict()
    seq_id_to_seq_name_map = dict()
    seq_error_log = list()
    file_seq_type = None
    last_seq_name = None
    mash_last_seq_name = None
    seq_count = 0
    max_seq_len = 0

    check_seq_len = user_params is not None # if configurations exist

    if check_seq_len:
        kmer_size = user_params.kmer_size
        default_dna_kmer_size = user_params.default_dna_kmer_size
        default_protein_kmer_size = user_params.default_protein_kmer_size


    # read in the file
    with open(seq_file_path, 'r') as f:

        for seq_record in SeqIO.parse(f, 'fasta'): # parse using fasta file format
            # error messages
            unknown_seq_type_msg = "Cannot determine whether it is DNA or protein sequence"
            short_seq_len_msg = "Sequence length shorter than the required k-mer size"

            last_seq_name = seq_record.description
            seq_type = check_seq_type(str(seq_record.seq))

            # if the sequence type can not be determined
            if seq_type is None:                 
                seq_error_log.append(unknown_seq_type_msg.format(last_seq_name))

            # ensure the file sequence type remains constant
            if file_seq_type is None:
                file_seq_type = seq_type
            elif file_seq_type != seq_type:
                seq_error_log.append('The input sequences consist of both DNA and protein sequences')
                break
            
            seq_len = len(seq_record.seq) # get and store the length of the current sequence

            # ensure the sequence length is valid
            if check_seq_len: # if configuration exists
                # MAJOR CHANGE: moved the logic if statements for _check_kmer_size into the original function
                if kmer_size is not None:
                    is_invalid_seq_len = seq_len < kmer_size
                    true_kmer_size = kmer_size
                elif seq_type == DNA:
                    is_invalid_seq_len = seq_len < default_dna_kmer_size
                    true_kmer_size = default_dna_kmer_size
                else:
                    is_invalid_seq_len = seq_len < default_protein_kmer_size
                    true_kmer_size = default_protein_kmer_size

                if is_invalid_seq_len:
                    seq_error_log.append(short_seq_len_msg.format(seq_record.description, true_kmer_size))

            # update new maximum sequence
            if seq_len > max_seq_len:
                max_seq_len = seq_len


            # MAJOR CHANGE: get and store name of the sequence within read_seq_file function
            m = re.match(FASTA_SEQ_NAME_WITH_COMMENT_PATTERN, last_seq_name)
            if m:
                mash_last_seq_name = '{}{}{}'.format(m.group(1), MASH_COMMENT_FIELD_SEP, m.group(2))
            else:
                mash_last_seq_name = last_seq_name

            mash_seq_name_to_seq_id_map[mash_last_seq_name] = seq_count
            seq_id_to_seq_name_map[str(seq_count)] = last_seq_name
            seq_count += 1

    if file_seq_type is None:
        seq_error_log.append('No valid sequences found in \'{}\''.format(seq_file_path))

    return SeqFileInfo(mash_seq_name_to_seq_id_map, seq_id_to_seq_name_map, seq_error_log, seq_file_path, \
                        file_seq_type, mash_last_seq_name, seq_count, max_seq_len)
