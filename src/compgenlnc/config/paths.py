from pathlib import Path
from typing import Final


BASE_PATH: Final[Path] = Path(__file__).resolve().parent.parent


DATA_FOLDER = BASE_PATH / 'data'
EMBEDDING_FOLDER = DATA_FOLDER / 'embedding'
INTERACTIONS_FOLDER = DATA_FOLDER / 'interactions'
LNCRNA_FOLDER = DATA_FOLDER / 'lncrna'
MIRNA_FOLDER = DATA_FOLDER / 'mirna'

SEQ_SS_PREFIX = 'seq+ss_join_'

LNCRNA_FASTA = LNCRNA_FOLDER / 'fasta' / 'lncRNA.fa'
LNCRNA_SS_FASTA = LNCRNA_FOLDER / 'fasta' / 'lncrna_2d.fa'
LNCRNA_SEQ_SS_FOLDER = LNCRNA_FOLDER / 'characterization' / 'seq+ss_join'

PRE_MIRNA_FASTA = MIRNA_FOLDER / 'fasta' / 'miRNA.fa'
MAT_MIRNA_FASTA = MIRNA_FOLDER / 'fasta' / 'mature_miRNA.fa'
MIRNA_SS_FASTA = MIRNA_FOLDER / 'fasta' / 'miRNA_2d.fa'
MIRNA_SEQ_SS_FOLDER = MIRNA_FOLDER / 'characterization' / 'seq+ss_join'

PAIRS_FILE_POSITIVE = INTERACTIONS_FOLDER / 'validated_pairs_negative.csv'
PAIRS_FILE_NEGATIVE = INTERACTIONS_FOLDER / 'validated_pairs_positive.csv'
PAIRS_FILES = [PAIRS_FILE_POSITIVE, PAIRS_FILE_NEGATIVE]
