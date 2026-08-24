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
)
from compgenlnc.utils.fasta_manager import iter_seq_fasta, list_seq_ss


def filter_fasta(
        filename: str | os.PathLike,
        filter: list[str],
        new_fasta: str | os.PathLike
) -> None:
    """Filter the RNA sequences from a FASTA file by the name.

    Args:
        filename: path of the FASTA file.
        filter: a list with all substrings contained in the name of the 
            desired sequences.
        new_fasta: path of the filtered FASTA file.
    """
    new_fasta = Path(new_fasta).resolve()
        
    fasta = [record for record in iter_seq_fasta(filename, filter)]
    fasta.sort(key=lambda record: record.id)

    new_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(fasta, new_fasta, 'fasta')


def filter_pre_mirna(
        filename: str | os.PathLike,
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
        filename: str | os.PathLike,
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
        filename: str | os.PathLike,
        filter: list[str]
) -> None:
    """Filter the lncRNA sequences from a FASTA file by the name.
    
    Args:
        filename: path of the FASTA file.
        filter: a list with all substrings contained in the name of the 
            desired sequences.
    """
    filter_fasta(filename, filter, new_fasta=LNCRNA_FASTA)



def seq_ss_to_fasta(
        folder: str | os.PathLike,
        new_fasta: str | os.PathLike
) -> None:
    folder = Path(str(folder)).resolve()
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
