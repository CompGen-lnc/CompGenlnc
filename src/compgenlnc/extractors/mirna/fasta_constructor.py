import os
from typing import Iterator

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord

from compgenlnc.config.paths import *


def iter_mirna(
        filename: str,
        filter: list[str] | None = None
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
            mir = title.split()[0]
            if (filter and
                not any(prefix in mir for prefix in filter)):
                continue
            yield SeqRecord(Seq(seq), mir, description='')


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
        
    fasta = [record for record in iter_mirna(filename, filter)]
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

    fasta = [record for record in iter_mirna(filename, filter)]
    fasta.sort(key=lambda record: record.id)

    new_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(fasta, new_fasta, 'fasta')



def read_sequence(path, mature_only=False):
    """Lee la linea 1 del archivo seq+ss_join_*.txt y la normaliza a
    alfabeto ACGT.
    """
    with open(path) as f:
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
            mir = os.path.splitext(filename[len(SEQ_SS_PREFIX):])[0]
            if keep is not None and mir not in keep:
                continue
            seq = read_sequence(folder / filename, mature_only)
            yield SeqRecord(Seq(seq), mir, description='')


def seq_ss_to_fasta(
        folder: str | Path = MIRNA_SEQ_SS_FOLDER,
        new_fasta: str | Path = PRE_MIRNA_FASTA
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
