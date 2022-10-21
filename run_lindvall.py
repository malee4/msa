# AUTHOR: Oscar Bulancea Lindvall, KTH Royal Institute of Technology
# CODE VERSION: 2019, last updated locally June 16, 2022
# LANGUAGE: Python, Leap IDE
# SOURCE: "Quantum Methods for Sequence Alignment and Metagenomics"
# URL: https://www.diva-portal.org/smash/get/diva2:1345195/FULLTEXT02

import Lindvall
# from data_formats import SeqQuery
import numpy as np
import dimod
import neal
import pickle

from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite


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

"""
Samples: Number of samples taken on annealer
match_cost and mismatch_cost: Weight for matching or mismatching base. Rewards and costs (recall this is a minimization function, reward = positive)
Simulation: flag determining to run simulated or quantum annealer
"""
def run_lindvall(sequence_string_set, old_center = None, samples = 1000, match_cost = -1, mismatch_cost = 1, simulation = False): 
    # if an old_center is provided, append to the sequence_string_set
    if old_center:
        sequence_string_set = [old_center.seq] + sequence_string_set
    # quantum spin column version
    sizes = [len(sequence_string_set[i]) for i in range(len(sequence_string_set))]

    # calc weights for matching
    matchings = np.zeros((len(sequence_string_set), max(sizes), len(sequence_string_set), max(sizes)))
    for s1 in range(len(sequence_string_set)):
        for s2 in range(len(sequence_string_set)):
            for n1 in range(sizes[s1]):
                for n2 in range(sizes[s2]):

                    if sequence_string_set[s1][n1] == sequence_string_set[s2][n2]:
                        matchings[s1, n1, s2, n2] = match_cost
                    else:
                        matchings[s1, n1, s2, n2] = mismatch_cost
    inserts = 0
    gap_penalty = 0
    params = {"gap_pen": gap_penalty, "extra_inserts": inserts}
    mat, shift, rev_inds = Lindvall.get_MSA_qubitmat(sizes, matchings,\
    gap_pen=gap_penalty, extra_inserts=inserts)

    h, J = mat_to_dimod_format(mat)

    if not simulation:
        cont = input("Continue with processing on annealer? y/n\n")
        if cont != "y":
            exit()
    
    bqm = dimod.BinaryQuadraticModel(h, J, shift, dimod.BINARY)
    print("Model formed")
    if simulation:
        print("Solver set-up for simulation")
        solver = neal.SimulatedAnnealingSampler()
    else:
        print("Solver set-up for QC")
        solver = EmbeddingComposite(DWaveSampler())
    
    response = solver.sample(bqm, num_reads = samples)
    return response

