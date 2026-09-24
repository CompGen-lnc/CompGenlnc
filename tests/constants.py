from pathlib import Path

BASE_FOLDER = Path(__file__).resolve().parent

TESTDATA_FOLDER = BASE_FOLDER / "testdata"
EMBEDDING_FOLDER = TESTDATA_FOLDER / "embedding"
INTERACTIONS_FOLDER = TESTDATA_FOLDER / "interactions"
LNCRNA_FOLDER = TESTDATA_FOLDER / "lncrna"
MIRNA_FOLDER = TESTDATA_FOLDER / "mirna"

EXAMPLE_FOLDER = MIRNA_FOLDER

# Examples
SEQ_FASTA_EXAMPLE = EXAMPLE_FOLDER / "mirna_pre.fa"
SS_FASTA_EXAMPLE = EXAMPLE_FOLDER / "mirna_2d.fa"
SEQ_SS_FOLDER = EXAMPLE_FOLDER / "characterization" / "seq+ss_join"
SS_LOOPS_FOLDER = EXAMPLE_FOLDER / "characterization" / "ss_loops"

# Interactions
INTERACTIONS_CSV = INTERACTIONS_FOLDER / "interactions.csv"
POS_INTERACTIONS_CSV = INTERACTIONS_FOLDER / "positive_interactions.csv"
NEG_INTERACTIONS_CSV = INTERACTIONS_FOLDER / "negative_interactions.csv"

# lncRNA
LNCRNA_FASTA = LNCRNA_FOLDER / "lncrna.fa"
LNCRNA_SS_FASTA = LNCRNA_FOLDER / "lncrna_2d.fa"

# miRNA
PRE_MIRNA_FASTA = MIRNA_FOLDER / "mirna_pre.fa"
MAT_MIRNA_FASTA = MIRNA_FOLDER / "mirna_mat.fa"
MIRNA_SS_FASTA = MIRNA_FOLDER / "mirna_2d.fa"
