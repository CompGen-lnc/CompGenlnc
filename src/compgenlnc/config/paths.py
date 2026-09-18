from pathlib import Path
from typing import Final

SEQ_SS_PREFIX = "seq+ss_join_"
LOOPS_PREFIX = "loops_"

BASE_PATH: Final[Path] = Path(__file__).resolve().parent.parent

# Main data folders
DATA_FOLDER = BASE_PATH / "data"
EMBEDDING_FOLDER = DATA_FOLDER / "embedding"
INTERACTIONS_FOLDER = DATA_FOLDER / "interactions"
LNCRNA_FOLDER = DATA_FOLDER / "lncrna"
MIRNA_FOLDER = DATA_FOLDER / "mirna"

# Embedding
MODELS_FOLDER = EMBEDDING_FOLDER / "models"
LNCRNA_DOC2VEC_MODEL = MODELS_FOLDER / "doc2vec.model"
PRE_MIRNA_DOC2VEC_MODEL = MODELS_FOLDER / "doc2vec.model"
MAT_MIRNA_DOC2VEC_MODEL = MODELS_FOLDER / "doc2vec.model"

# Interactions
PAIRS_FILE = INTERACTIONS_FOLDER / "interactions.csv"
PAIRS_FILE_POSITIVE = INTERACTIONS_FOLDER / "validated_pairs_negative.csv"
PAIRS_FILE_NEGATIVE = INTERACTIONS_FOLDER / "validated_pairs_positive.csv"
PAIRS_FILES = [PAIRS_FILE_POSITIVE, PAIRS_FILE_NEGATIVE]

# lncRNA
LNCRNA_FASTA = LNCRNA_FOLDER / "fasta" / "lncRNA.fa"
LNCRNA_SS_FASTA = LNCRNA_FOLDER / "fasta" / "lncrna_2d.fa"
LNCRNA_SEQ_SS_FOLDER = LNCRNA_FOLDER / "characterization" / "seq+ss_join"
LNCRNA_LOOPS_FOLDER = LNCRNA_FOLDER / "characterization" / "ss_loops"
LNCRNA_DICTS_FOLDER = LNCRNA_FOLDER / "dicts"

# miRNA
PRE_MIRNA_FASTA = MIRNA_FOLDER / "fasta" / "miRNA.fa"
MAT_MIRNA_FASTA = MIRNA_FOLDER / "fasta" / "mature_miRNA.fa"
MIRNA_SS_FASTA = MIRNA_FOLDER / "fasta" / "miRNA_2d.fa"
MIRNA_SEQ_SS_FOLDER = MIRNA_FOLDER / "characterization" / "seq+ss_join"
MIRNA_LOOPS_FOLDER = MIRNA_FOLDER / "characterization" / "ss_loops"
PRE_MIRNA_DICTS_FOLDER = MIRNA_FOLDER / "dicts" / "precursor"
MAT_MIRNA_DICTS_FOLDER = MIRNA_FOLDER / "dicts" / "mature"

# Dictionaries
KMER_DICT = "kmer.dict"
CTD_DICT = "ctd.dict"
DOC2VEC_DICT = "doc2vec.dict"
ROLE2VEC_DICT = "role2vec.dict"
NORMAL_SEQ_DICT = "normal_seq.dict"
