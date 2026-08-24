import os
from pathlib import Path
from typing import Iterator

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord

from compgenlnc.config.paths import (
    SEQ_SS_PREFIX,
)


def iter_seq_fasta(
        filename: str | os.PathLike,
        filter: list[str] | None = None
) -> Iterator[SeqRecord]:
    """Iterate all the RNA sequences from a FASTA file that fit the 
    filter.

    Args:
        filename: Path of the FASTA file.
        filter: A list with all substrings contained in the name of the 
            desired sequences. Iterate all if None.

    Yields:
        SeqRecord: A record with the sequence and its name.
    """
    filename = Path(filename).resolve()

    with open(filename, 'r') as file:
        for title, seq in SimpleFastaParser(file):
            id_ = title.split()[0]
            if (filter and
                not any(prefix in id_ for prefix in filter)):
                continue
            yield SeqRecord(Seq(seq), id_, description='')


def read_sequence(
        filename: str | os.PathLike,
        mature_only: bool = False
) -> str:
    """Read the first line of a seq+ss file and normalize to ACGT
    alphabet.

    Args:
        filename: Path of the seq+ss file.
        mature_only: Only saves the mature section of the sequence if
            True.
    """
    filename = Path(filename).resolve()

    with open(filename) as f:
        seq = f.readline().strip()
    if mature_only:
        seq = ''.join(c for c in seq if c.isupper())
    return seq.upper().replace('U', 'T')


def iter_seq_ss(
        folder: str | os.PathLike,
        keep: list[str] | None = None,
        mature_only: bool = False
) -> Iterator[SeqRecord]:
    """Iterate all the seq+ss files from the same folder that fit the
    filter.

    Args:
        folder: path of the folder containing seq+ss files.
        keep: list of the RNA names to iterate.
        mature_only: Only saves the mature section of the sequence if
            True.
    """
    folder = Path(folder).resolve()

    for filename in sorted(os.listdir(folder)):
        if not filename.startswith(SEQ_SS_PREFIX):
            continue
        id_ = os.path.splitext(filename[len(SEQ_SS_PREFIX):])[0]
        if keep is not None and id_ not in keep:
            continue
        seq = read_sequence(folder / filename, mature_only)
        yield SeqRecord(Seq(seq), id_, description='')