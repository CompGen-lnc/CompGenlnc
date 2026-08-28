from collections import Counter
from itertools import product
from pathlib import Path
from typing import Iterable
import os

import numpy as np

from compgenlnc.config.paths import SEQ_SS_PREFIX
from compgenlnc.structs import SeqRecord, SeqSSRecord
from compgenlnc.utils.dict_manager import save_dict
from compgenlnc.utils.fasta_manager import read_sequence, iter_seq_ss

def get_ctd(
        seq: str | SeqRecord | SeqSSRecord
) -> np.typing.NDArray[np.float32]:
    """Get the CTD of a RNA sequences from the given sequence.

    Args:
        seq: The RNA sequence string on ATCG alphabet.

    Returns:
        NDArray: A concatenated numpy array with the CTD of the
            sequence.
    """
    if isinstance(seq, SeqSSRecord):
        seq = str(seq.seq)
    if isinstance(seq, SeqRecord):
        seq = str(seq)
    if not seq:
        return np.array([0.0] * 30, np.float32)
    
    seq = seq.upper().replace('U', 'T')
    codes = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    n = len(seq)

    nums = Counter(seq)
    num_A, num_T, num_G, num_C = nums['A'], nums['T'], nums['G'], nums['C']
    trans = Counter(seq[i : i + 2] for i in range(n - 1))
    AT_trans = trans['AT'] + trans['TA']
    AG_trans = trans['AG'] + trans['GA']
    AC_trans = trans['AC'] + trans['CA']
    TG_trans = trans['TG'] + trans['GT']
    TC_trans = trans['TC'] + trans['CT']
    GC_trans = trans['GC'] + trans['CG']

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
    [A0_dis, A1_dis, A2_dis, A3_dis, A4_dis] = dist[0]
    [T0_dis, T1_dis, T2_dis, T3_dis, T4_dis] = dist[1]
    [G0_dis, G1_dis, G2_dis, G3_dis, G4_dis] = dist[2]
    [C0_dis, C1_dis, C2_dis, C3_dis, C4_dis] = dist[3]

    return np.array([
        num_A / n, num_T / n, num_G / n, num_C / n,
        AT_trans / (n - 1), AG_trans / (n - 1), AC_trans / (n - 1),
        TG_trans / (n - 1), TC_trans / (n - 1), GC_trans / (n - 1),
        A0_dis, A1_dis, A2_dis, A3_dis, A4_dis,
        T0_dis, T1_dis, T2_dis, T3_dis, T4_dis,
        G0_dis, G1_dis, G2_dis, G3_dis, G4_dis,
        C0_dis, C1_dis, C2_dis, C3_dis, C4_dis,
    ], np.float32)


def get_ctd_by_name(
        id_: str,
        folder: str | os.PathLike,
        /,
        mature_only: bool = False
) -> np.typing.NDArray[np.float32]:
    """Get the CTD of a RNA sequences from the seq+ss file of the given
    molecule.

    Args:
        id_: Name of the RNA sequence.
        folder: Path of the folder where the seq+ss file is.
        mature_only: Only saves the mature section of the sequence if
            True.

    Returns:
        NDArray: A concatenated numpy array with the CTD of the
            sequence.
    """
    filename = folder / f'{SEQ_SS_PREFIX}{id_}.dat'
    seq = read_sequence(filename, mature_only=mature_only)
    return get_ctd(seq)


def gen_ctd_dict(
        record_list: Iterable[SeqSSRecord],
        dict_file: str | os.PathLike,
        /,
        keep: Iterable[str] | None = None,
        mature_only: bool = False
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
    dict_file = Path(dict_file).resolve()
    
    dict_folder = dict_file.parent
    dict_folder.mkdir(parents=True, exist_ok=True)

    ctd_dict = {record.id: get_ctd(record.seq) for record in record_list}
    save_dict(dict_file, ctd_dict)


def gen_ctd_dict_from_folder(
        folder: str | os.PathLike,
        dict_file: str | os.PathLike,
        /,
        keep: Iterable[str] | None = None,
        mature_only: bool = False
) -> None:
    """Get the CTD of a RNA sequences from the seq+ss file of the given
    molecule.

    Args:
        folder: Path of the folder where the seq+ss file is.
        dict_file: Path of the file where the dictionary will be saved.
        keep: List of the RNA names to save. Save all if None.
        mature_only: Only saves the mature section of the sequence if
            True.
    """
    folder = Path(folder).resolve()
    dict_file = Path(dict_file).resolve()

    dict_folder = dict_file.parent
    dict_folder.mkdir(parents=True, exist_ok=True)

    ctd_dict = {
        record.id: get_ctd(record.seq)
        for record in iter_seq_ss(folder, keep=keep, mature_only=mature_only)
    }
    save_dict(dict_file, ctd_dict)
