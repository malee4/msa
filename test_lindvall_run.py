import test_lindvall_columns
#from data_formats import SeqQuery 
import numpy as np
import dimod
import neal
import pickle

from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

samples = 1000 # number of samples to run 
save_file = "filename"

# sequences = ["AT", "T"]
# sequences = ["AT", "T", "T", "A"]
sequences = ["MVRMMN", "MNVRLMMN", "NVRLVMMN", "PVRMMN", "MVRMMNLMM", "VRMMLMM"] # write your sequences here
simulation = True # flag determining to run simulated or quantum annealer

# rewards (negative) and costs (positive)
match_cost = -1
mismatch_cost = 1


sizes = [len(sequences[i]) for i in range(len(sequences))]
# calculate weights for matching
matchings = np.zeros((len(sequences), max(sizes), len(sequences), max(sizes))) 
for s1 in range(len(sequences)):
    for s2 in range(len(sequences)):
        for n1 in range(sizes[s1]): 
            for n2 in range(sizes[s2]):
                if sequences[s1][n1] == sequences[s2][n2]: 
                    matchings[s1,n1,s2,n2] = match_cost
                else:
                    matchings[s1,n1,s2,n2] = mismatch_cost

inserts = 0
gap_penalty = 0
params = {"gap_pen": gap_penalty, "extra_inserts": inserts}
mat, shift, rev_inds = test_lindvall_columns.get_MSA_qubitmat(sizes, matchings,\
gap_pen=gap_penalty, extra_inserts=inserts)

def mat_to_dimod_format(matrix): 
    """interprets the QUBO matrix of interactions into separate format for linear and interaction terms """
    n = matrix.shape[0]
    linear = {}
    interaction = {} 
    for i in range(n):
        linear[i] = matrix[i,i] 
        for j in range(i+1,n):
            interaction[(i,j)] = matrix[i,j] 
    return linear, interaction

h, J = mat_to_dimod_format(mat) 
if not simulation:
    cont = input("Continue? y/n ") 
    if cont != "y":
        exit()

bqm = dimod.BinaryQuadraticModel(h,J, shift, dimod.BINARY) 
if simulation:
    solver = neal.SimulatedAnnealingSampler() 
else:
    # WARNING! requires setup of D-wave certificate file beforehand
    solver = EmbeddingComposite(DWaveSampler())
response = solver.sample(bqm, num_reads=samples)
positions = response.lowest().samples()[0]

# here's the brute force assignment of elements to positions part

# get the number of items per alignment
total_len = 0
all_seq = ""


def get_alignment_string(sequence_string_set, gaps, positions):
    # group positions based on sequence
    string_size = max([len(s) for s in sequence_string_set]) + gaps

    # get the positions, based on sequence
    organized_positions = get_positions(string_size, sequence_string_set, positions)

    # create an empty matrix
    align_strings = [["-" for i in range(string_size)] for i in range(len(sequence_string_set))]

    # fill in the matrix
    for key in organized_positions.keys():
        for result_id in range(len(organized_positions[key])):
            if organized_positions[key][result_id]:
                align_strings[key[0]][result_id] = sequence_string_set[key[0]][key[1]]
    return align_strings

def get_positions(string_size, sequence_string_set, positions):
    count = 0
    # split into positions
    organized_positions = dict()
    
    for seq_number in range(len(sequence_string_set)):
        for i in range(len(sequence_string_set[seq_number])):
            # create list of items
            temp = list()
            start = count
            for j in range(count, count + string_size):
                count = count + 1
                temp = temp + [positions[j]]
            organized_positions[(seq_number, i)] = temp
    return organized_positions

total_len = 0 
max_len = -1
for s in sequences:
    total_len += len(s)
    if len(s) > max_len:
        max_len = len(s)

max_seq_len = len(positions) // total_len
gaps = max_seq_len - max_len
print(max_seq_len)

output = get_alignment_string(sequences, gaps, positions)
for item in output:
    print(item)