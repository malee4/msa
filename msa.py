# import packages
import numpy as np # for pre- and post-processing
from dimod import BinaryQuadraticModel
from dwave.system import DWaveSampler, EmbeddingComposite

# insert inputs
strands = [[]] # iterable that stores DNA strands as lists
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
for s in range(len(strands)):
    for i in range(len(strands) - s):

## penalty score

