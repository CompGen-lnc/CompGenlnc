from pathlib import Path

BASE_FOLDER = Path(__file__).resolve().parent

TESTDATA_FOLDER = BASE_FOLDER / "testdata"
SEQ_FASTA_EXAMPLE = TESTDATA_FOLDER / "fasta_ejemplo.fa"
SS_FASTA_EXAMPLE = TESTDATA_FOLDER / "fasta_2d_ejemplo.fa"

SEQ_SS_FOLDER = TESTDATA_FOLDER / "characterization" / "seq+ss_join"
SS_LOOPS_FOLDER = TESTDATA_FOLDER / "characterization" / "ss_loops"
