# import packages
import numpy as np # for pre- and post-processing
from dimod import BinaryQuadraticModel
from dwave.system import DWaveSampler, EmbeddingComposite

# insert inputs
strands = [[]] # iterable that stores DNA strands as lists
# number of columns 
columns = 20 # currently some arbitrary value

# create pairwise weights

# final weights matrix; strand 1 x original index of base x strand 2 x original index of base [weight]
weights = [
    [
        [
            [
                # insert weight
            ]
        ]
    ]
]

# final locations matrix; strand x original index of base x final column
table = [
    [
        [
            # 
        ]
    ]
]

# create empty model
bqm = BinaryQuadraticModel('BINARY') # {0, 1}

# incorporate objective function
## rewarding term
for s_1 in range(len(strands)):
    for s_2 in range(len(strands) - s_1):
        for n_1 in range(len(strands[s_1])):
            for n_2 in range(len(strands[s_2])):
                for i in range(columns):
                    print()
        
## penalty score

