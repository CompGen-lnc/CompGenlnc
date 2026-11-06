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


def get_ctd(seq: SeqLike, /) -> NDArray[np.float32]:
    """Get the CTD of a RNA sequences from the given sequence.

        seq: The RNA sequence string on ATCG alphabet.

    Returns:
    Returns
    -------
    NDArray
        A concatenated numpy array with the CTD of the sequence.
    """
    if isinstance(seq, SeqSSRecord):
        seq = str(seq.seq)
    if isinstance(seq, SeqRecord):
        seq = str(seq)
    if not seq:
        return np.array([0.0] * 30, np.float32)

    seq = seq.upper().replace("T", "U")
    codes = {"A": 0, "U": 1, "G": 2, "C": 3}
    n = len(seq)

    nums = Counter(seq)
    num_A, num_U, num_G, num_C = nums["A"], nums["U"], nums["G"], nums["C"]
    nums_array = np.array([num_A, num_U, num_G, num_C], np.float32) / n

    trans = Counter(seq[i : i + 2] for i in range(n - 1))
    AU_trans = trans["AU"] + trans["UA"]
    AG_trans = trans["AG"] + trans["GA"]
    AC_trans = trans["AC"] + trans["CA"]
    UG_trans = trans["UG"] + trans["GU"]
    UC_trans = trans["UC"] + trans["CU"]
    GC_trans = trans["GC"] + trans["CG"]
    trans_array = np.array(
        [AU_trans, AG_trans, AC_trans, UG_trans, UC_trans, GC_trans],
        np.float32,
    ) / (n - 1)

    count = [0] * 4
    dist = [[0.0] * 5 for i in range(4)]
    for i in range(n):
        actual_num = nums[seq[i]]
        indx = codes[seq[i]]
        count[indx] += 1
        actual_count = count[indx]
        if actual_count == 1:
            dist[indx][0] = ((i * 1.0) + 1) / n
        if actual_count == int(round(actual_num / 4.0)):
            dist[indx][1] = ((i * 1.0) + 1) / n
        if actual_count == int(round(actual_num / 2.0)):
            dist[indx][2] = ((i * 1.0) + 1) / n
        if actual_count == int(round((actual_num * 3 / 4.0))):
            dist[indx][3] = ((i * 1.0) + 1) / n
        if actual_count == actual_num:
            dist[indx][4] = ((i * 1.0) + 1) / n

    for i, j in product(range(4), range(1, 5)):
        if dist[i][j] < dist[i][j - 1]:
            dist[i][j] = dist[i][j - 1]
    dist_array = np.array(dist, np.float32).flatten()

    return np.concatenate((nums_array, trans_array, dist_array))


def get_ctd_by_name(folder: StrPath, id_: str) -> NDArray[np.float32]:
    """Get the CTD of a RNA sequences from the seq+ss file of the given
    molecule.

    Args:
        id_: Name of the RNA sequence.
        folder: Path of the folder where the seq+ss file is.
        mature_only: Only saves the mature section of the sequence if
    Parameters
    ----------
    folder :
        Path of the folder where the seq+ss file is.
    id_: Name of the RNA sequence.

    Returns:
        NDArray: A concatenated numpy array with the CTD of the
            sequence.
    """
    folder = Path(folder).resolve()
    seq = read_sequence(folder, id_)
    return get_ctd(seq)


def gen_ctd_dict(
    dict_file: StrPath,
    records: Collection[SeqSSRecord],
    *,
    overwrite: bool = False,
) -> None:
    """Get the CTD of a RNA sequences from the seq+ss file of the given
    molecule.

    Args:
        record_list: Iterable that contains all the sequence records to
            save into the dictionary.
        dict_file: Path of the file where the dictionary will be saved.
        keep: List of the RNA names to save. Save all if None.
        mature_only: Only saves the mature section of the sequence if
            True.
    """
    send_msg("Generating CTD dictionary\n")
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
    desc = "  Calculating CTD"
    records = track_progress(records, desc, total, "seq", initial)

    ctd_dict = {record.id: get_ctd(record.seq) for record in records}
    save_dict(dict_file, ctd_dict, overwrite=overwrite)

    send_msg("  CTD dictionary generated\n")


def gen_ctd_dict_from_fasta(
    dict_file: StrPath,
    fasta_file: StrPath,
    *,
    keep: Collection[str] | None = None,
    exclude: Collection[str] | None = None,
    overwrite: bool = False,
) -> None:
    """Get the ctd of a RNA sequences from 1 to k from a FASTA file with RNA
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
    records = load_fasta(fasta_file, mode=SEQ_MODE, keep=keep, exclude=exclude)
    gen_ctd_dict(dict_file, records, overwrite=overwrite)


def gen_ctd_dict_from_folder(
    dict_file: StrPath,
    seq_ss_folder: StrPath,
    *,
    keep: Collection[str] | None = None,
    exclude: Collection[str] | None = None,
    overwrite: bool = False,
) -> None:
    """Get the CTD of a RNA sequences from the seq+ss files.

    Args:
        folder: Path of the folder where the seq+ss file is.
        dict_file: Path of the file where the dictionary will be saved.
        keep: List of the RNA names to save. Save all if None.
        mature_only: Only saves the mature section of the sequence if
            True.
    """
    records = load_seq_ss(seq_ss_folder, keep=keep, exclude=exclude)
    gen_ctd_dict(dict_file, records, overwrite=overwrite)
