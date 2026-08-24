from pathlib import Path

BASE_FOLDER = Path(__file__).resolve().parent

FIXTURE_FOLDER = BASE_FOLDER / 'fixtures'
SEQ_FASTA_EXAMPLE = FIXTURE_FOLDER / 'fasta_ejemplo.fa'
SS_FASTA_EXAMPLE = FIXTURE_FOLDER / 'fasta_2d_ejemplo.fa'

SEQ_SS_FOLDER = FIXTURE_FOLDER / 'characterization' / 'seq+ss_join'