from __future__ import annotations

from collections import Counter
from collections.abc import Collection
from itertools import product
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from compgenlnc.consts.params import SEQ_MODE
from compgenlnc.fileman.dict_manager import load_dict, save_dict
from compgenlnc.fileman.fasta_manager import (
    load_fasta,
    load_seq_ss,
    read_sequence,
)
from compgenlnc.fileman.log_manager import send_msg, track_progress
from compgenlnc.structs.seq_record import SeqRecord
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.molecules import SeqLike

if TYPE_CHECKING:
    from _typeshed import StrPath


def get_kmer(seq: SeqLike, /, k: int) -> NDArray[np.float32]:
    """Get the k-mers of a RNA sequence.

    Parameters
    ----------
    seq : SeqLike
        The RNA sequence string on ATCG alphabet.
    k : int
        Maximum length for k-mers. It calculates k-mers with length from
        1 to `k`, included.

    Returns
    -------
    ndarray of float32
        A 1D array with the k-mers of the sequence.
    """
    if isinstance(seq, SeqSSRecord):
        seq = str(seq.seq)
    if isinstance(seq, SeqRecord):
        seq = str(seq)

    seq = seq.upper().replace("T", "U")

    # Get the kmers of length n for the sequence
    def get_nmer(seq: str, n: int) -> NDArray[np.float32]:
        if not seq or len(seq) < n:
            return np.array([0] * 4**n, np.float32)
        kmers = Counter((seq[i : i + n] for i in range(len(seq) - n + 1)))
        combinations = ("".join(p) for p in product("AUCG", repeat=n))
        return np.array(
            [kmers[sub_s] / kmers.total() for sub_s in combinations],
            np.float32,
        )

    return np.concatenate([get_nmer(seq, n) for n in range(1, k + 1)])


def get_kmer_by_name(dir: StrPath, id_: str, k: int) -> NDArray[np.float32]:
    """Get the k-mers of a RNA sequence from the seq+ss file.

    Parameters
    ----------
    dir : str or os.PathLike
        Path of the directory where the seq+ss file is.
    id_ : str
        ID of the RNA sequence.
    k : int
        Maximum length for k-mers. It calculates k-mers with length from
        1 to `k`, included.

    Returns
    -------
    ndarray of float32
        A 1D array with the k-mers of the sequence.
    """
    dir = Path(dir).resolve()

    seq = read_sequence(dir, id_)
    return get_kmer(seq, k)


def gen_kmer_dict(
    dict_file: StrPath,
    records: Collection[SeqSSRecord],
    k: int,
    *,
    overwrite: bool = False,
) -> None:
    """Generate a dictionary of k-mers from a set of RNA sequences.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    records : collection of SeqSSRecord
        Set that contains all the sequence records to save into the
        dictionary.
    k : int
        Maximum length for k-mers. It calculates k-mers with length from
        1 to `k`, included.
    overwrite : bool, optional
        Append the items to the existing dictionary if `True`, ignoring
        those that have already been processed. Otherwise, replace the
        existing dictionary. Default: `False`.
    """
    send_msg("Generating k-mer dictionary\n")
    records = list(records)
    total = len(records)

    dict_file = Path(dict_file).resolve()
    if dict_file.exists():
        send_msg("  Analyzing existing dictionary...")

        existing_dict = load_dict(dict_file)
        records = [
            record for record in records if record.id not in existing_dict
        ]

        send_msg("\r  Existing dictionary analyzed    \n")

    initial = total - len(records)
    desc = "  Calculating k-mers"
    records = track_progress(records, desc, total, "seq", initial)

    kmer_dict = {record.id: get_kmer(record.seq, k) for record in records}
    save_dict(dict_file, kmer_dict, overwrite=overwrite)

    send_msg("  k-mer dictionary generated\n")


def gen_kmer_dict_from_fasta(
    dict_file: StrPath,
    fasta_file: StrPath,
    k: int,
    *,
    keep: Collection[str] | None = None,
    exclude: Collection[str] | None = None,
    overwrite: bool = False,
) -> None:
    """Generate a dictionary of k-mers from a FASTA file.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    fasta_file : str or os.PathLike
        Path of the FASTA file with the sequences.
    k : int
        Maximum length for k-mers. It calculates k-mers with length from
        1 to `k`, included.
    keep : collection of str, optional
        Set of the RNA IDs to save. Do not apply if `None` or empty.
        Default: `None`.
    exclude : collection of str, optional
        Set of the RNA IDs to ignore. Do not apply if `None` or empty.
        Default: `None`.
    overwrite : bool, optional
        Append the items to the existing dictionary if `True`, ignoring
        those that have already been processed. Otherwise, replace the
        existing dictionary. Default: `False`.
    """
    records = load_fasta(
        fasta_file,
        mode=SEQ_MODE,
        keep=keep,
        exclude=exclude,
    )
    gen_kmer_dict(dict_file, records, k, overwrite=overwrite)


def gen_kmer_dict_from_dir(
    dict_file: StrPath,
    seq_ss_dir: StrPath,
    k: int,
    *,
    keep: Collection[str] | None = None,
    exclude: Collection[str] | None = None,
    overwrite: bool = False,
) -> None:
    """Generate a dictionary of k-mers from seq+ss files.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    seq_ss_dir : str or os.PathLike
        Path of the directory where the seq+ss files are.
    k : int
        Maximum length for k-mers. It calculates k-mers with length from
        1 to `k`, included.
    keep : collection of str, optional
        Set of the RNA IDs to save. Do not apply if `None` or empty.
        Default: `None`.
    exclude : collection of str, optional
        Set of the RNA IDs to ignore. Do not apply if `None` or empty.
        Default: `None`.
    overwrite : bool, optional
        Append the items to the existing dictionary if `True`, ignoring
        those that have already been processed. Otherwise, replace the
        existing dictionary. Default: `False`.
    """
    records = load_seq_ss(seq_ss_dir, keep=keep, exclude=exclude)
    gen_kmer_dict(dict_file, records, k, overwrite=overwrite)
