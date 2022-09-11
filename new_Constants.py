MIN_CLUSTER_PROCESSING_SIZE = 2
GAPS = 0
# is_simulation = True

# SEQUENCE CONSTANTS
DNA = 1
AA = 2

IUPAC_DNA_STR_PATTERN = r'^[ACGTRYSWKMBDHVN]+$'
IUPAC_AA_STR_PATTERN = r'^[ABCDEFGHIKLMNPQRSTVWXYZ]+$'
IUPAC_AMBIG_DNA_BASES = r'[RYSWKMBDHVN]'
IUPAC_AMBIG_AA_BASES = r'[X]'

FASTA_SEQ_NAME_WITH_COMMENT_PATTERN = r'^(\S*)\s?(.*)$'
MASH_COMMENT_FIELD_SEP = ':'

# for intermediate processing between clustering and MSA
MAX_LOCAL_CLUSTER_SIZE = -1

# from settings
# EstimatedSimilarity
ESTIMATED_SIMILARITY_HIGH = 0.95
ESTIMATED_SIMILARITY_LOW = 0.75
ESTIMATED_SIMILARITY_STEPSIZE = 0.025

# Threshold
THRESHOLD_PRECLUSTER = 20000

# DNAMash
DNAMash_Kmer = 7
DNAMash_Sketch = 2000

# ProteinMash
ProteinMash_KMER = 5
ProteinMash_Sketch = 2000

# NoiseFilter
NOISE_FILTER_MARGIN = 0.2

# DNAEvaluation
DNA_MatchScore = 2
DNA_MismatchPenalty = 3
DNA_GapOpeningPenalty = 5
DNA_GapExtensionPenalty = 2

# ProteinEvaluation
Protein_ScoreMatrix = "blosum62"
Protein_GapOpeningPenalty = 11
Protein_GapExtensionPenalty = 1

MIN_CLUSTER_PROCESSING_SIZE = 2
GAPS = 0
SIMULATION_SETTING = True

def getConfig():
    # TODO: Implement configurations validation

    config = dict()

    # settings for estimated similarity metrics
    estSim = dict()
    estSim['high'] = ESTIMATED_SIMILARITY_HIGH
    estSim['low'] = ESTIMATED_SIMILARITY_LOW
    estSim['stepsize'] = ESTIMATED_SIMILARITY_STEPSIZE
    config['estimated_similarity'] = "estSim"

    config['threshold_precluster'] = THRESHOLD_PRECLUSTER

    config['kmer_size'] = None

    dna_mash = dict()
    dna_mash['kmer'] = DNAMash_Kmer
    dna_mash['sketch'] = DNAMash_Sketch
    config['dna_mash'] = dna_mash
    
    protein_mash = dict()
    protein_mash['kmer'] = ProteinMash_KMER
    protein_mash['sketch'] = ProteinMash_Sketch
    config['protein_mash'] = protein_mash

    config['noise_filter_margin'] = NOISE_FILTER_MARGIN

    dna_eval = dict()
    dna_eval['match_score'] = DNA_MatchScore
    dna_eval['mismatch_penalty'] = DNA_MismatchPenalty
    dna_eval['gap_opening_penalty'] = DNA_GapOpeningPenalty
    dna_eval['gap_extension_penalty'] = DNA_GapExtensionPenalty
    config['dna_evaluation'] = dna_eval

    protein_eval = dict()
    protein_eval['score_matrix'] = Protein_ScoreMatrix
    protein_eval['gap_opening_penalty'] = Protein_GapOpeningPenalty
    protein_eval['gap_extension_penalty'] = Protein_GapExtensionPenalty
    config['protein_evaluation'] = protein_eval

    return config


