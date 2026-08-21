import os
from pathlib import Path
from typing import Iterator

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord

from compgenlnc.config.paths import (
    LNCRNA_FASTA,
    MAT_MIRNA_FASTA,
    PAIRS_FILES,
    PRE_MIRNA_FASTA,
    SEQ_SS_PREFIX,
)


def iter_seq_fasta(
        filename: str | Path,
        filter: list[str] | None = None
) -> Iterator[SeqRecord]:
    """Return all the RNA sequences from a FASTA file that fit the 
    filter.

    Args:
        filename: path of the FASTA file.
        filter: a list with all substrings contained in the name of the 
            desired sequences.

    Yields:
        SeqRecord: a record with the sequence and its name.
    """
    if not isinstance(filename, Path):
        filename = Path(filename).resolve()

    with open(filename, 'r') as file:
        for title, seq in SimpleFastaParser(file):
            id_ = title.split()[0]
            if (filter and
                not any(prefix in id_ for prefix in filter)):
                continue
            yield SeqRecord(Seq(seq), id_, description='')


def filter_fasta(
        filename: str | Path,
        filter: list[str],
        new_fasta: str | Path
) -> None:
    """Filter the RNA sequences from a FASTA file by the name.

    Args:
        filename: path of the FASTA file.
        filter: a list with all substrings contained in the name of the 
            desired sequences.
        new_fasta: path of the filtered FASTA file.
    """
    if not isinstance(new_fasta, Path):
        new_fasta = Path(new_fasta).resolve()
        
    fasta = [record for record in iter_seq_fasta(filename, filter)]
    fasta.sort(key=lambda record: record.id)

    new_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(fasta, new_fasta, 'fasta')


def filter_pre_mirna(
        filename: str | Path,
        filter: list[str]
) -> None:
    """Filter the precursor miRNA sequences from a FASTA file by the 
    name.

    Args:
        filename: path of the FASTA file.
        filter: a list with all substrings contained in the name of the 
            desired sequences.
    """
    filter_fasta(filename, filter, new_fasta=PRE_MIRNA_FASTA)


def filter_mat_mirna(
        filename: str | Path,
        filter: list[str]
) -> None:
    """Filter the mature miRNA sequences from a FASTA file by the name.
    
    Args:
        filename: path of the FASTA file.
        filter: a list with all substrings contained in the name of the 
            desired sequences.
    """
    filter_fasta(filename, filter, new_fasta=MAT_MIRNA_FASTA)


def filter_lncrna(
        filename: str | Path,
        filter: list[str]
) -> None:
    """Filter the lncRNA sequences from a FASTA file by the name.
    
    Args:
        filename: path of the FASTA file.
        filter: a list with all substrings contained in the name of the 
            desired sequences.
    """
    filter_fasta(filename, filter, new_fasta=LNCRNA_FASTA)



def read_sequence(filename: str | Path, mature_only=False) -> str:
    """Read the first line of the file and normalized to ACGT alphabet.

    Args:
        filename:
    """
    with open(filename) as f:
        seq = f.readline().strip()
    if mature_only:
        seq = ''.join(c for c in seq if c.isupper())
    return seq.upper().replace('U', 'T')


def list_seq_ss(
        folder: str | Path,
        new_fasta: str | Path,
        /,
        keep: list[str] | None = None,
        mature_only=False
) -> Iterator[SeqRecord]:
    with open(new_fasta, 'w') as out_file:
        for filename in sorted(os.listdir(folder)):
            if not filename.startswith(SEQ_SS_PREFIX):
                continue
            id_ = os.path.splitext(filename[len(SEQ_SS_PREFIX):])[0]
            if keep is not None and id_ not in keep:
                continue
            seq = read_sequence(folder / filename, mature_only)
            yield SeqRecord(Seq(seq), id_, description='')


def seq_ss_to_fasta(
        folder: str | Path,
        new_fasta: str | Path
) -> None:
    if not isinstance(folder, Path):
        folder = Path(str(folder)).resolve()
    if not isinstance(new_fasta, Path):
        new_fasta = Path(str(new_fasta)).resolve()

    needed_seq = []
    for pf in PAIRS_FILES:
        if pf.exists():
            for line in open(pf):
                parts = line.strip().split(',')
                if len(parts) >= 2:
                    needed_seq.add(parts[1])

    fasta = [
        record
        for record in list_seq_ss(folder, new_fasta, keep=needed_seq)
    ]
    fasta.sort(key=lambda record: record.id)

    new_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(fasta, new_fasta, 'fasta')
