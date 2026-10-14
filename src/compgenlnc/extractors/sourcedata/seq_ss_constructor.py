import os
from collections.abc import Iterable
from pathlib import Path

from compgenlnc.config.paths import (
    LNCRNA_FASTA,
    LNCRNA_SEQ_SS_FOLDER,
    MIRNA_SEQ_SS_FOLDER,
    PRE_MIRNA_FASTA,
    SEQ_SS_PREFIX,
)
from compgenlnc.fileman.fasta_manager import load_fasta


def join_seq_ss(
    seq: str, ss: str, id_: str, folder: str | os.PathLike
) -> None:
    """Save the sequence and structure of a RNA molecule in the same
    file.

    Args:
        seq: The RNA sequence on ACGT alphabet.
        ss: The RNA 2d structure on dot-bracket notation.
        id_: Name of the RNA molecule.
        folder: Path for the seq+ss file of the RNA molecule.
    """
    folder = Path(folder).resolve()

    path = folder / f"{SEQ_SS_PREFIX}{id_}.dat"
    with open(path, "w") as out_file:
        out_file.write(seq + "\n")
        out_file.write(ss)


def fasta_to_seq_ss(
    seq_file: str | os.PathLike,
    ss_file: str | os.PathLike,
    folder: str | os.PathLike,
    /,
    *,
    keep: Iterable[str] | None = None,
) -> int:
    """Convert a FASTA file into seq+ss files for each molecule.

    Args:
        seq_file: Path of the file containing all the RNA sequences.
        ss_file: Path of the file containing all the RNA structures.
        folder: Path for the seq+ss files.
    """
    seq_file = Path(seq_file).resolve()
    ss_file = Path(ss_file).resolve()
    folder = Path(folder).resolve()

    seq_dict = {
        record.id: str(record.seq).upper().replace("U", "T")
        for record in load_fasta(seq_file, mode="seq")
    }
    ss_dict = {
        record.id: str(record.seq)
        for record in load_fasta(ss_file, mode="seq")
    }
    id_list = seq_dict.keys() & ss_dict.keys()
    if keep is not None:
        id_list &= keep
    folder.mkdir(parents=True, exist_ok=True)
    total = 0

    for id_ in id_list:
        try:
            seq = seq_dict[id_]
            ss = ss_dict[id_]
            join_seq_ss(seq, ss, id_, folder)
            total += 1
        except:
            pass

    return total


def mirna_to_seq_ss() -> int:
    """Convert the pre-miRNA FASTA file into seq+ss files for each
    molecule.
    """
    return fasta_to_seq_ss(PRE_MIRNA_FASTA, MIRNA_SEQ_SS_FOLDER)


def lncrna_to_seq_ss() -> int:
    """Convert the lncRNA FASTA file into seq+ss files for each
    molecule.
    """
    return fasta_to_seq_ss(LNCRNA_FASTA, LNCRNA_SEQ_SS_FOLDER)
