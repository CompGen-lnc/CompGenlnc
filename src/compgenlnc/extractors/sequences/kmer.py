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
    load_seq_ss,
    load_fasta,
    read_sequence,
)
from compgenlnc.fileman.log_manager import send_msg, track_progress
from compgenlnc.structs.seq_record import SeqRecord
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.molecules import SeqLike

if TYPE_CHECKING:
    from _typeshed import StrPath


def get_kmer(seq: SeqLike, /, k: int) -> NDArray[np.float32]:
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


def get_kmer_by_name(folder: StrPath, id_: str, k: int) -> NDArray[np.float32]:
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

    seq = read_sequence(folder, id_)
    return get_kmer(seq, k)


def gen_kmer_dict(
    dict_file: StrPath,
    records: Collection[SeqSSRecord],
    k: int,
    *,
    overwrite: bool = False,
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
    records = load_fasta(
        fasta_file,
        mode=SEQ_MODE,
        keep=keep,
        exclude=exclude,
    )
    gen_kmer_dict(dict_file, records, k, overwrite=overwrite)


def gen_kmer_dict_from_folder(
    dict_file: StrPath,
    seq_ss_folder: StrPath,
    k: int,
    *,
    keep: Collection[str] | None = None,
    exclude: Collection[str] | None = None,
    overwrite: bool = False,
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
    records = load_seq_ss(seq_ss_folder, keep=keep, exclude=exclude)
    gen_kmer_dict(dict_file, records, k, overwrite=overwrite)
