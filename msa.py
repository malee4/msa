# import packages
import numpy as np # for pre- and post-processing
from dimod import BinaryQuadraticModel
from dwave.system import DWaveSampler, EmbeddingComposite

# insert inputs
strands = [[]] # iterable that stores DNA strands as lists
# create pairwise weights

# final weights matrix
weights = [[[]]]

# final locations matrix
table = [[[]]]

# create empty model
bqm = BinaryQuadraticModel('BINARY') # {0, 1}

# incorporate objective function
## rewarding term

## penalty score

# incorporate constraints
## 

# sum of x_s,i = |s|

# sum of all x in one column >= 1