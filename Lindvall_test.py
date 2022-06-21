# AUTHOR: Oscar Bulancea Lindvall, KTH Royal Institute of Technology
# CODE VERSION: 2019, last updated locally June 16, 2022
# LANGUAGE: Python, Leap IDE
# SOURCE: "Quantum Methods for Sequence Alignment and Metagenomics"
# URL: https://www.diva-portal.org/smash/get/diva2:1345195/FULLTEXT02

import lindvall
# from data_formats import SeqQuery
import numpy as np
import dimod
import neal
import pickle

from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

samples = 1000
save_file = "msa_output"

sequences = ["AT", "T", "A", "A"] # put your input here
simulation = False # flag determining to run simulated or quantum annealer

# rewards and costs (recall this is a minimization function, reward = positive)
match_cost = -1
mismatch_cost = 1

# quantum spin column version
sizes = [len(sequences[i]) for i in range(len(sequences))]
# calc weights for matching
matchings = np.zeros((len(sequences), max(sizes), len(sequences), max(sizes)))
for s1 in range(len(sequences)):
    for s2 in range(len(sequences)):
        for n1 in range(sizes[1]):
            for n2 in range(sizes[s2]):
                if sequences[s1][n1] == sequences[s2][n2]:
                    matchings[s1, n1, s2, n2] = match_cost
                else:
                    matchings[s1, n1, s2, n2] = mismatch_cost

inserts = 0
gap_penalty = 0
params = {"gap_pen": gap_penalty, "extra_inserts": inserts}
mat, shift, rev_inds = Lindvall.get_MSA_qubitmat(sizes, matchings,\
gap_pen=gap_penalty, extra_inserts=inserts)

def mat_to_dimod_format(matrix): 
    """
    interprets the QUBO matrix of interactions into separate format for linear and interaction terms 
    """
    n = matrix.shape[0]
    linear = {}
    interaction = {}
    for i in range(n):
        linear[i] = matrix[i, i]
        for j in range(i+1, n):
            interaction[(i, j)] = matrix[i, j]
    return linear, interaction

h, J = mat_to_dimod_format(mat)

if not simulation:
    cont = input("continue: y/n\n")
    if cont != "y":
        exit()

bqm = dimod.BinaryQuadraticModel(h, J, shift, dimod.BINARY)
if simulation:
    solver = neal.SimulatedAnnealingSampler()
else:
    solver = EmbeddingComposite(DWaveSampler())

response = solver.sample(bqm, num_reads = samples)