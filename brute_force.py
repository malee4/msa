# import packages
import numpy as np # for pre- and post-processing
from dimod import BinaryQuadraticModel
from dwave.system import DWaveSampler, EmbeddingComposite

# inputs: 2D np array with pairwise scores

# objective: select one pairing for each