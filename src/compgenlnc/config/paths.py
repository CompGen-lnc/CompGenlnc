from pathlib import Path
from typing import Final


BASE_PATH: Final[Path] = Path(__file__).resolve().parent.parent


DATA_FOLDER = BASE_PATH / 'data'
EMBEDDING_FOLDER = DATA_FOLDER / 'embedding'
INTERACTIONS_FOLDER = DATA_FOLDER / 'interactions'
LNCRNA_FOLDER = DATA_FOLDER / 'lncrna'
MIRNA_FOLDER = DATA_FOLDER / 'mirna'

SEQ_SS_PREFIX = 'seq+ss_join'

PRE_MIRNA_FASTA = MIRNA_FOLDER / 'fasta' / 'mirna.fa'
MAT_MIRNA_FASTA = MIRNA_FOLDER / 'fasta' / 'mature_mirna.fa'
MIRNA_SEQ_SS = MIRNA_FOLDER / 'characterization' / 'seq+ss_join'

PAIRS_FILE_POSITIVE = INTERACTIONS_FOLDER / 'validated_pairs_negative.csv'
PAIRS_FILE_NEGATIVE = INTERACTIONS_FOLDER / 'validated_pairs_positive.csv'
PAIRS_FILES = [PAIRS_FILE_POSITIVE, PAIRS_FILE_NEGATIVE]
