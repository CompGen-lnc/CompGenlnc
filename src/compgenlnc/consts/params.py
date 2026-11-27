import sys
from multiprocessing import cpu_count

# RNA type
LNCRNA = "lnc"
"""Indicate to process lncRNAs."""
MIRNA = "mir"
"""Indicate to process miRNAs."""
PRECURSOR = "pre"
"""Indicate to process precursor miRNAs."""
MATURE = "mat"
"""Indicate to process mature miRNAs."""
INTERACTION = "inter"
"""Indicate to process interactions."""

# record mode
SEQ_MODE = "seq"
"""Indicate to process sequences."""
SS_MODE = "ss"
"""Indicate to process 2D structures."""

# workers
WORKERS = cpu_count() - 2 if "pytest" not in sys.modules else 0
"""Number of processes to run on multiprocessing tasks."""

# url
LFS_URL = (
    "https://media.githubusercontent.com/media/CompGen-lnc/CompGenlnc/main"
)
"""URL for package repository git LFS."""
GIT_LNCRNA_FASTA = f"{LFS_URL}/basefiles/lncRNA.fa"
"""URL lncRNA sequence FASTA in package repository."""
GIT_MAT_MIRNA_FASTA = f"{LFS_URL}/basefiles/mature_miRNA.fa"
"""URL mature miRNA sequence FASTA in package repository."""
GIT_PRE_MIRNA_FASTA = f"{LFS_URL}/basefiles/miRNA.fa"
"""URL precursor miRNA sequence FASTA in package repository."""
GIT_DATA_MIRNA = f"{LFS_URL}/basefiles/data_miRNA.csv"
"""URL precursor-mature miRNA association CSV in package repository."""
GIT_POSITIVE_PAIRS = f"{LFS_URL}/basefiles/positive_pairs.csv"
"""URL positive interactions CSV in package repository."""
