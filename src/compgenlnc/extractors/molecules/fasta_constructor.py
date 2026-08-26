import os
from pathlib import Path
from typing import Iterable

from Bio import SeqIO

from compgenlnc.config.paths import (
    LNCRNA_FASTA,
    MAT_MIRNA_FASTA,
    PAIRS_FILES,
    PRE_MIRNA_FASTA,
)
from compgenlnc.utils.fasta_manager import (
    iter_seq_fasta,
    iter_seq_ss,
    save_fasta,
)


def filter_fasta(
        filename: str | os.PathLike,
        new_fasta: str | os.PathLike,
        filter: Iterable[str]
) -> None:
    """Filter the RNA sequences from a FASTA file by the name.

    Args:
        filename: Path of the FASTA file.
        filter: A list with all substrings contained in the name of the 
            desired sequences.
        new_fasta: Path of the filtered FASTA file.
    """
    filename = Path(filename).resolve()
    new_fasta = Path(new_fasta).resolve()
        
    fasta = [record for record in iter_seq_fasta(filename, filter)]
    fasta.sort(key=lambda record: record.id)

    new_fasta.parent.mkdir(parents=True, exist_ok=True)
    save_fasta(new_fasta, fasta, 'seq')


def filter_pre_mirna(
        filename: str | os.PathLike,
        filter: Iterable[str]
) -> None:
    """Filter the precursor miRNA sequences from a FASTA file by the 
    name.

    Args:
        filename: Path of the FASTA file.
        filter: A list with all substrings contained in the name of the 
            desired sequences.
    """
    filename = Path(filename).resolve()
    filter_fasta(filename, filter, new_fasta=PRE_MIRNA_FASTA)


def filter_mat_mirna(
        filename: str | os.PathLike,
        filter: Iterable[str]
) -> None:
    """Filter the mature miRNA sequences from a FASTA file by the name.
    
    Args:
        filename: Path of the FASTA file.
        filter: A list with all substrings contained in the name of the 
            desired sequences.
    """
    filename = Path(filename).resolve()
    filter_fasta(filename, filter, new_fasta=MAT_MIRNA_FASTA)


def filter_lncrna(
        filename: str | os.PathLike,
        filter: Iterable[str]
) -> None:
    """Filter the lncRNA sequences from a FASTA file by the name.
    
    Args:
        filename: Path of the FASTA file.
        filter: A list with all substrings contained in the name of the 
            desired sequences.
    """
    filename = Path(filename).resolve()
    filter_fasta(filename, filter, new_fasta=LNCRNA_FASTA)



def seq_ss_to_fasta(
        folder: str | os.PathLike,
        new_seq_fasta: str | os.PathLike,
        new_ss_fasta: str | os.PathLike
) -> None:
    """Join all the seq+ss files from the same folder into a FASTA file.

    Args:
        new_seq_fasta: Path for the new FASTA file with the sequences.
        folder: Path of the folder containing all the seq+ss files.
    """
    folder = Path(folder).resolve()
    new_seq_fasta = Path(new_seq_fasta).resolve()

    needed_seq = []
    for pf in PAIRS_FILES:
        if pf.exists():
            for line in open(pf):
                parts = line.strip().split(',')
                if len(parts) >= 2:
                    needed_seq.add(parts[1])

    fasta = [
        record
        for record in iter_seq_ss(folder, keep=needed_seq)
    ]
    fasta.sort(key=lambda record: record.id)

    new_seq_fasta.parent.mkdir(parents=True, exist_ok=True)
    save_fasta(new_seq_fasta, fasta, 'seq')
