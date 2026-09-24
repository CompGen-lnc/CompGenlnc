import multiprocessing as mp
import os
from collections.abc import Iterable
from pathlib import Path

from .structure_predictor import (
    extract_2d_structure_identified,
)
from compgenlnc.config.paths import (
    LNCRNA_FASTA,
    LNCRNA_SS_FASTA,
    MAT_MIRNA_FASTA,
    MIRNA_SS_FASTA,
    PRE_MIRNA_FASTA,
)
from compgenlnc.utils.fasta_manager import (
    iter_seq_fasta,
    iter_seq_ss,
    load_fasta,
    save_fasta,
)


def filter_fasta(
    filename: str | os.PathLike,
    new_fasta: str | os.PathLike,
    filter: Iterable[str],
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
    save_fasta(new_fasta, fasta, mode="seq")


def filter_pre_mirna(
    filename: str | os.PathLike, filter: Iterable[str]
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
    filename: str | os.PathLike, filter: Iterable[str]
) -> None:
    """Filter the mature miRNA sequences from a FASTA file by the name.

    Args:
        filename: Path of the FASTA file.
        filter: A list with all substrings contained in the name of the
            desired sequences.
    """
    filename = Path(filename).resolve()
    filter_fasta(filename, filter, new_fasta=MAT_MIRNA_FASTA)


def filter_lncrna(filename: str | os.PathLike, filter: Iterable[str]) -> None:
    """Filter the lncRNA sequences from a FASTA file by the name.

    Args:
        filename: Path of the FASTA file.
        filter: A list with all substrings contained in the name of the
            desired sequences.
    """
    filename = Path(filename).resolve()
    filter_fasta(filename, filter, new_fasta=LNCRNA_FASTA)


def gen_fasta_2d(
    seq_fasta: str | os.PathLike,
    ss_fasta: str | os.PathLike,
    /,
    *,
    keep: Iterable[str] | None = None,
) -> None:
    num_processes = max(1, mp.cpu_count() - 2)
    seq_fasta = Path(seq_fasta).resolve()
    ss_fasta = Path(ss_fasta).resolve()

    ss_fasta.parent.mkdir(parents=True, exist_ok=True)

    with mp.Pool(processes=num_processes) as pool:
        ss_list = pool.imap(
            extract_2d_structure_identified,
            load_fasta(seq_fasta, mode="seq", keep=keep),
            chunksize=10,
        )
        save_fasta(ss_fasta, ss_list, mode="ss")


def seq_ss_to_fasta(
    folder: str | os.PathLike,
    new_seq_fasta: str | os.PathLike,
    new_ss_fasta: str | os.PathLike,
    /,
    *,
    keep: Iterable[str] | None = None,
) -> None:
    """Join all the seq+ss files from the same folder into a FASTA file.

    Args:
        new_seq_fasta: Path for the new FASTA file with the sequences.
        folder: Path of the folder containing all the seq+ss files.
    """
    folder = Path(folder).resolve()
    new_seq_fasta = Path(new_seq_fasta).resolve()
    new_ss_fasta = Path(new_ss_fasta).resolve()

    new_seq_fasta.parent.mkdir(parents=True, exist_ok=True)
    new_ss_fasta.parent.mkdir(parents=True, exist_ok=True)
    save_fasta(new_seq_fasta, iter_seq_ss(folder, keep=keep), mode="seq")
    save_fasta(new_ss_fasta, iter_seq_ss(folder, keep=keep), mode="ss")


if __name__ == "__main__":
    gen_fasta_2d(LNCRNA_FASTA, LNCRNA_SS_FASTA)
    gen_fasta_2d(PRE_MIRNA_FASTA, MIRNA_SS_FASTA)
