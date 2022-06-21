# AUTHOR:Jimmy Ka Ho Chiu & Rick Twee-Hee Ong
# CODE VERSION: 2022, last updated locally June 21, 2022
# LANGUAGE: Python
# SOURCE: "Clustering biological sequences with dynamic sequence similarity threshold"
# URL: https://doi.org/10.1186/s12859-022-04643-9

DNA = 1
AA = 2

IUPAC_DNA_STR_PATTERN = r'^[ACGTRYSWKMBDHVN]+$'
IUPAC_AA_STR_PATTERN = r'^[ABCDEFGHIKLMNPQRSTVWXYZ]+$'
IUPAC_AMBIG_DNA_BASES = r'[RYSWKMBDHVN]'
IUPAC_AMBIG_AA_BASES = r'[X]'

FASTA_SEQ_NAME_WITH_COMMENT_PATTERN = r'^(\S*)\s?(.*)$'
MASH_COMMENT_FIELD_SEP = ':'