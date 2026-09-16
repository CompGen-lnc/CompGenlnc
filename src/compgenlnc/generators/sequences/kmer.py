import os
from collections import Counter
from collections.abc import Iterable
from itertools import product
from pathlib import Path

import numpy as np

from compgenlnc.config.paths import SEQ_SS_PREFIX
from compgenlnc.structs import SeqRecord, SeqSSRecord
from compgenlnc.typing import SeqLike
from compgenlnc.utils.dict_manager import save_dict
from compgenlnc.utils.fasta_manager import (
    iter_seq_ss,
    load_fasta,
    read_sequence,
)


def get_kmer(seq: SeqLike, k: int) -> np.typing.NDArray[np.float32]:
    """Get the k-mers of a RNA sequences from 1 to k from the given
    sequence.

    Args:
        seq: The RNA sequence string on ATCG alphabet.
        k: Maximum length for k-mers. It calculates k-mers with length
            from 1 to k, included.

    Returns:
        NDArray: A concatenated numpy array with the k-mers of the
            sequence.
    """
    # Get the kmers of length n for the sequence
    if isinstance(seq, SeqSSRecord):
        seq = str(seq.seq)
    if isinstance(seq, SeqRecord):
        seq = str(seq)

    seq = seq.upper().replace("U", "T")

    def get_nmer(seq: str, n: int) -> np.typing.NDArray[np.float32]:
        if not seq or len(seq) < n:
            return np.array([0] * 4**n, np.float32)
        kmers = Counter((seq[i : i + n] for i in range(len(seq) - n + 1)))
        combinations = ("".join(p) for p in product("ATCG", repeat=n))
        return np.array(
            [kmers[sub_s] / kmers.total() for sub_s in combinations],
            np.float32,
        )

    return np.concatenate([get_nmer(seq, n) for n in range(1, k + 1)])


def get_kmer_by_name(
    id_: str, folder: str | os.PathLike, k: int, /, mature_only: bool = False
) -> np.typing.NDArray[np.float32]:
    """Get the k-mers of a RNA sequences from 1 to k from the seq+ss
    file of the given molecule.

    Args:
        id_: Name of the RNA sequence.
        folder: Path of the folder where the seq+ss file is.
        k: Maximum length for k-mers. It calculates k-mers with length
            from 1 to k, included.
        mature_only: Only saves the mature section of the sequence if
            True.

    Returns:
        NDArray: A concatenated numpy array with the k-mers of the
            sequence.
    """
    folder = Path(folder).resolve()

    filename = folder / f"{SEQ_SS_PREFIX}{id_}.dat"
    seq = read_sequence(filename, mature_only=mature_only)
    return get_kmer(seq, k)


def gen_kmer_dict(
    record_list: Iterable[SeqSSRecord],
    dict_file: str | os.PathLike,
    k: int,
    /,
    mature_only: bool = False,
) -> None:
    """Get the k-mers of RNA sequences from 1 to k for the given molecule.

    Args:
        record_list: Iterable that contains all the sequence records to
            save into the dictionary.
        dict_file: Path of the file where the dictionary will be saved.
        k: Maximum length for k-mers. It calculates k-mers with length
            from 1 to k, included.
        keep: List of the RNA names to save. Save all if None.
        mature_only: Only saves the mature section of the sequence if
            True.
    """
    dict_file = Path(dict_file).resolve()

    dict_folder = dict_file.parent
    dict_folder.mkdir(parents=True, exist_ok=True)

    kmer_dict = {record.id: get_kmer(record.seq, k) for record in record_list}
    save_dict(dict_file, kmer_dict)


def gen_kmer_dict_from_fasta(
    filename: str | os.PathLike,
    dict_file: str | os.PathLike,
    k: int,
    /,
    keep: Iterable[str] | None = None,
    mature_only: bool = False,
) -> None:
    """Get the k-mers of RNA sequences from 1 to k from a FASTA file with RNA
    sequences.

    Args:
        filename: Path of the FASTA file.
        dict_file: Path of the file where the dictionary will be saved.
        k: Maximum length for k-mers. It calculates k-mers with length
            from 1 to k, included.
        keep: List of the RNA names to save. Save all if None.
        mature_only: Only saves the mature section of the sequence if
            True.
    """
    filename = Path(filename).resolve()
    dict_file = Path(dict_file).resolve()

    dict_folder = dict_file.parent
    dict_folder.mkdir(parents=True, exist_ok=True)

    kmer_dict = {
        record.id: get_kmer(record.seq, k)
        for record in load_fasta(filename, mode="seq")
        if not keep or record.id in keep
    }
    save_dict(dict_file, kmer_dict)


def gen_kmer_dict_from_folder(
    folder: str | os.PathLike,
    dict_file: str | os.PathLike,
    k: int,
    /,
    keep: Iterable[str] | None = None,
    mature_only: bool = False,
) -> None:
    """Get the k-mers of RNA sequences from 1 to k from the seq+ss
    files.

    Args:
        folder: Path of the folder where the seq+ss file is.
        dict_file: Path of the file where the dictionary will be saved.
        k: Maximum length for k-mers. It calculates k-mers with length
            from 1 to k, included.
        keep: List of the RNA names to save. Save all if None.
        mature_only: Only saves the mature section of the sequence if
            True.
    """
    folder = Path(folder).resolve()
    dict_file = Path(dict_file).resolve()

    dict_folder = dict_file.parent
    dict_folder.mkdir(parents=True, exist_ok=True)

    kmer_dict = {
        record.id: get_kmer(record.seq, k)
        for record in iter_seq_ss(folder, keep=keep, mature_only=mature_only)
    }
    save_dict(dict_file, kmer_dict)
