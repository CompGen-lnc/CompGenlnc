import sys
from multiprocessing import cpu_count

# molecule type
LNCRNA = "lnc"
MIRNA = "mir"
PRECURSOR = "pre"
MATURE = "mat"
INTERACTION = "inter"

# record mode
SEQ_MODE = "seq"
SS_MODE = "ss"

# workers
WORKERS = cpu_count() - 2 if "pytest" not in sys.modules else 0

# url
LFS_URL = (
    "https://media.githubusercontent.com/media/CompGen-lnc/CompGenlnc/git-lfs-basefiles"
)
GIT_LNCRNA_FASTA = f"{LFS_URL}/basefiles/lncRNA.fa"
GIT_MAT_MIRNA_FASTA = f"{LFS_URL}/basefiles/mature_miRNA.fa"
GIT_PRE_MIRNA_FASTA = f"{LFS_URL}/basefiles/miRNA.fa"
GIT_DATA_MIRNA = f"{LFS_URL}/basefiles/data_miRNA.csv"
GIT_POSITIVE_PAIRS = f"{LFS_URL}/basefiles/validated_pairs_positive.csv"
