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


def get_ctd(seq: SeqLike, /) -> NDArray[np.float32]:
    """Get the CTD of a RNA sequence.

    Parameters
    ----------
    seq : SeqLike
        The RNA sequence string on AUCG alphabet.

    Returns
    -------
    numpy.ndarray of float32
        A 1D array with the CTD of the sequence.
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


def get_ctd_by_name(dir: StrPath, id_: str) -> NDArray[np.float32]:
    """Get the CTD of a RNA sequence from the seq+ss file.

    Parameters
    ----------
    dir : str or os.PathLike
        Path of the directory where the seq+ss file is.
    id_ : str
        ID of the RNA sequence to read.

    Returns
    -------
    numpy.ndarray of float32
        A 1D array with the CTD of the sequence.
    """
    dir = Path(dir).resolve()
    seq = read_sequence(dir, id_)
    return get_ctd(seq)


def gen_ctd_dict(
    dict_file: StrPath,
    records: Collection[SeqSSRecord],
    *,
    overwrite: bool = False,
) -> None:
    """Generate a dictionary of CTD from an RNA sequence set.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    records : collection of SeqSSRecord
        Set containing all the sequence records to save into the
        dictionary.
    overwrite : bool, optional
        Replace the existing dictionary if `True` in case there is.
        Otherwise, append the items to the existing dictionary, ignoring
        those that have already been processed. Default: `False`.
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
    """Generate a dictionary of CTD from a FASTA file.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    fasta_file : str or os.PathLike
        Path of the FASTA file with the sequences.
    keep : collection of str, optional
        Set of IDs of the RNAs to save. Do not apply if `None` or empty.
        Default: `None`.
    exclude : collection of str, optional
        Set of IDs of the RNAs to ignore. Do not apply if `None` or
        empty. Default: `None`.
    overwrite : bool, optional
        Replace the existing dictionary if `True` in case there is.
        Otherwise, append the items to the existing dictionary, ignoring
        those that have already been processed. Default: `False`.
    """
    records = load_fasta(fasta_file, mode=SEQ_MODE, keep=keep, exclude=exclude)
    gen_ctd_dict(dict_file, records, overwrite=overwrite)


def gen_ctd_dict_from_dir(
    dict_file: StrPath,
    seq_ss_dir: StrPath,
    *,
    keep: Collection[str] | None = None,
    exclude: Collection[str] | None = None,
    overwrite: bool = False,
) -> None:
    """Generate a dictionary of CTD from seq+ss files.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    seq_ss_dir : str or os.PathLike
        Path of the directory where the seq+ss files are.
    keep : collection of str, optional
        Set of IDs of the RNAs to save. Do not apply if `None` or empty.
        Default: `None`.
    exclude : collection of str, optional
        Set of IDs of the RNAs to ignore. Do not apply if `None` or
        empty. Default: `None`.
    overwrite : bool, optional
        Replace the existing dictionary if `True` in case there is.
        Otherwise, append the items to the existing dictionary, ignoring
        those that have already been processed. Default: `False`.
    """
    records = load_seq_ss(seq_ss_dir, keep=keep, exclude=exclude)
    gen_ctd_dict(dict_file, records, overwrite=overwrite)
