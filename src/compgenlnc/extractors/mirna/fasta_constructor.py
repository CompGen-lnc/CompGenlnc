import os
from typing import Iterator

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord

from compgenlnc.config.paths import *


def iter_filtered_mirna(
        filename: str,
        filter: list[str]
) -> Iterator[SeqRecord]:
    """Return all the RNA sequences from a FASTA file that fit the 
    filter.

    Args:
        filename: the name of the FASTA file.
        filter: a list with all substrings contained in the name of the 
            desired sequences.

    Yields:
        SeqRecord: a record with the sequence and its name.
    """
    with open(filename, 'r') as out_file:
        for title, seq in SimpleFastaParser(out_file):
            id = title.split()[0]
            if not any(prefix in id for prefix in filter):
                continue
            yield SeqRecord(Seq(seq), id, description='')


def filter_pre_mirna(
        filename: str,
        filter: list[str],
        /,
        new_fasta: str | Path = PRE_MIRNA_FASTA
) -> None:
    """Filter the precursor miRNA sequences from a FASTA file by the 
    name.

    Args:
        filename: the name of the FASTA file.
        filter: a list with all substrings contained in the name of the 
            desired sequences.    
    """
    if not isinstance(new_fasta, Path):
        new_fasta = Path(str(new_fasta)).resolve()
        
    fasta = [record for record in iter_filtered_mirna(filename, filter)]
    fasta.sort(key=lambda record: record.id)

    new_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(fasta, new_fasta, 'fasta')


def filter_mat_mirna(
        filename: str,
        filter: list[str],
        /,
        new_fasta: str | Path = MAT_MIRNA_FASTA
) -> None:
    """Filter the mature miRNA sequences from a FASTA file by the 
    name.
    
    Args:
        filename: the name of the FASTA file.
        filter: a list with all substrings contained in the name of the 
            desired sequences.    
    """
    if not isinstance(new_fasta, Path):
        new_fasta = Path(str(new_fasta)).resolve()

    fasta = [record for record in iter_filtered_mirna(filename, filter)]
    fasta.sort(key=lambda record: record.id)

    new_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(fasta, new_fasta, 'fasta')


