import os
from pathlib import Path
from typing import Iterable

import numpy as np
from Bio.SeqRecord import SeqRecord

from compgenlnc.config.paths import SEQ_SS_PREFIX
from compgenlnc.utils.dict_manager import save_dict
from compgenlnc.utils.fasta_manager import read_sequence, iter_seq_ss

def get_ctd(seq: str) -> np.typing.NDArray[np.float64]:
    """Get the CTD of a RNA sequences from the given sequence.

    Args:
        seq: The RNA sequence string on ATCG alphabet.

    Returns:
        NDArray: A concatenated numpy array with the CTD of the
            sequence.
    """
    n = len(seq)
    num_A = num_T = num_G = num_C = 0.0
    AT_trans = AG_trans = AC_trans = TG_trans = TC_trans = GC_trans = 0.0
    for i in range(len(seq) - 1):
        if seq[i] == "A":
            num_A = num_A + 1
        if seq[i] == "T":
            num_T = num_T + 1
        if seq[i] == "G":
            num_G = num_G + 1
        if seq[i] == "C":
            num_C = num_C + 1
        if ((seq[i] == "A" and seq[i + 1] == "T") or
            (seq[i] == "T" and seq[i + 1] == "A")):
            AT_trans = AT_trans + 1
        if ((seq[i] == "A" and seq[i + 1] == "G") or
            (seq[i] == "G" and seq[i + 1] == "A")):
            AG_trans = AG_trans + 1
        if ((seq[i] == "A" and seq[i + 1] == "C") or
            (seq[i] == "C" and seq[i + 1] == "A")):
            AC_trans = AC_trans + 1
        if ((seq[i] == "T" and seq[i + 1] == "G") or
            (seq[i] == "G" and seq[i + 1] == "T")):
            TG_trans = TG_trans + 1
        if ((seq[i] == "T" and seq[i + 1] == "C") or
            (seq[i] == "C" and seq[i + 1] == "T")):
            TC_trans = TC_trans + 1
        if ((seq[i] == "G" and seq[i + 1] == "C") or
            (seq[i] == "C" and seq[i + 1] == "G")):
            GC_trans = GC_trans + 1

    a, t, g, c = 0, 0, 0, 0
    A0_dis, A1_dis, A2_dis, A3_dis, A4_dis = 0.0, 0.0, 0.0, 0.0, 0.0
    T0_dis, T1_dis, T2_dis, T3_dis, T4_dis = 0.0, 0.0, 0.0, 0.0, 0.0
    G0_dis, G1_dis, G2_dis, G3_dis, G4_dis = 0.0, 0.0, 0.0, 0.0, 0.0
    C0_dis, C1_dis, C2_dis, C3_dis, C4_dis = 0.0, 0.0, 0.0, 0.0, 0.0
    for i in range(len(seq) - 1):
        if seq[i] == "A":
            a = a + 1
            if a == 1:
                A0_dis = ((i * 1.0) + 1) / n
            if a == int(round(num_A / 4.0)):
                A1_dis = ((i * 1.0) + 1) / n
            if a == int(round(num_A / 2.0)):
                A2_dis = ((i * 1.0) + 1) / n
            if a == int(round((num_A * 3 / 4.0))):
                A3_dis = ((i * 1.0) + 1) / n
            if a == num_A:
                A4_dis = ((i * 1.0) + 1) / n
        if seq[i] == "T":
            t = t + 1
            if t == 1:
                T0_dis = ((i * 1.0) + 1) / n
            if t == int(round(num_T / 4.0)):
                T1_dis = ((i * 1.0) + 1) / n
            if t == int(round((num_T / 2.0))):
                T2_dis = ((i * 1.0) + 1) / n
            if t == int(round((num_T * 3 / 4.0))):
                T3_dis = ((i * 1.0) + 1) / n
            if t == num_T:
                T4_dis = ((i * 1.0) + 1) / n
        if seq[i] == "G":
            g = g + 1
            if g == 1:
                G0_dis = ((i * 1.0) + 1) / n
            if g == int(round(num_G / 4.0)):
                G1_dis = ((i * 1.0) + 1) / n
            if g == int(round(num_G / 2.0)):
                G2_dis = ((i * 1.0) + 1) / n
            if g == int(round(num_G * 3 / 4.0)):
                G3_dis = ((i * 1.0) + 1) / n
            if g == num_G:
                G4_dis = ((i * 1.0) + 1) / n
        if seq[i] == "C":
            c = c + 1
            if c == 1:
                C0_dis = ((i * 1.0) + 1) / n
            if c == int(round(num_C / 4.0)):
                C1_dis = ((i * 1.0) + 1) / n
            if c == int(round(num_C / 2.0)):
                C2_dis = ((i * 1.0) + 1) / n
            if c == int(round(num_C * 3 / 4.0)):
                C3_dis = ((i * 1.0) + 1) / n
            if c == num_C:
                C4_dis = ((i * 1.0) + 1) / n
    return (num_A / n, num_T / n, num_G / n, num_C / n,
            AT_trans / (n - 1), AG_trans / (n - 1), AC_trans / (n - 1),
            TG_trans / (n - 1), TC_trans / (n - 1), GC_trans / (n - 1),
            A0_dis, A1_dis, A2_dis, A3_dis, A4_dis,
            T0_dis, T1_dis, T2_dis, T3_dis, T4_dis,
            G0_dis, G1_dis, G2_dis, G3_dis, G4_dis,
            C0_dis, C1_dis, C2_dis, C3_dis, C4_dis)


def get_ctd_by_name(
        id_: str,
        folder: str | os.PathLike,
        /,
        mature_only: bool = False
) -> np.typing.NDArray[np.float64]:
    folder = Path(folder).resolve()
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
        record_list: Iterable[SeqRecord],
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
