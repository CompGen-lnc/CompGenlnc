from __future__ import annotations

import os
from collections.abc import Collection, Iterable
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import RNA
from numpy.typing import NDArray

from compgenlnc.consts.paths import LOOPS_PREFIX
from compgenlnc.fileman.dict_manager import load_dict, save_dict
from compgenlnc.fileman.fasta_manager import iter_seq_ss
from compgenlnc.fileman.log_manager import send_msg, track_progress
from compgenlnc.fileman.loops_manager import (
    read_2d_structure_loops,
    read_2d_structure_loops_from_dir,
)
from compgenlnc.structs.seq_ss_record import SeqSSRecord

if TYPE_CHECKING:
    from _typeshed import StrPath


def extract_2d_structure_loops(
    loops_dir: StrPath, record: SeqSSRecord
) -> float:
    """Calculate the loops of a RNA structure.

    Parameters
    ----------
    loops_dir : str or os.PathLike
        Path of the directory to save the loops file.
    record : SeqSSRecord
        The 2D structure record to calculate the loops.

    Returns
    -------
    float
        The structure MFE obtained from the loops.
    """
    loops_dir = Path(loops_dir).resolve()
    if loops_dir.exists():
        files = os.listdir(loops_dir)
        if f"{LOOPS_PREFIX}{record.id}.dat" in files:
            return read_2d_structure_loops(loops_dir, record.id).energy
    else:
        loops_dir.mkdir(parents=True, exist_ok=True)

    filename = loops_dir / f"{LOOPS_PREFIX}{record.id}.dat"
    with open(filename, "w") as out_file:
        try:
            fc = RNA.fold_compound(str(record.seq), None, RNA.OPTION_EVAL_ONLY)
            energy = fc.eval_structure_verbose(str(record.ss), out_file)
        except Exception:
            energy = 0
        out_file.write(f"{'Energy':<40}: {energy:>8.2f}")
    return energy


def extract_2d_structure_loops_from_list(
    loops_dir: StrPath,
    records: Iterable[SeqSSRecord],
) -> NDArray[np.float32]:
    """Calculate the loops of an RNA structure set.

    Parameters
    ----------
    loops_dir : str or os.PathLike
        Path of the directory to save the loops files.
    records : iterable object of SeqSSRecord
        Group of 2D structure records to calculate the loops.

    Returns
    -------
    numpy.ndarray of float
        The structure MFEs obtained from the loops.
    """
    energy_list = []
    for record in records:
        energy_list.append(extract_2d_structure_loops(loops_dir, record))

    return np.array(energy_list, dtype=np.float32)


def extract_2d_structure_loops_from_dir(
    loops_dir: StrPath,
    seq_ss_dir: StrPath,
    *,
    keep: Collection[str] | None = None,
    exclude: Collection[str] | None = None,
) -> NDArray[np.float32]:
    """Calculate the loops of an RNA structure set in seq+ss files.

    Parameters
    ----------
    loops_dir : str or os.PathLike
        Path of the directory to save the loops files.
    seq_ss_dir : str or os.PathLike
        Path of the directory where the seq+ss files are.
    keep : collection of str, optional
        Set of IDs of the RNAs to save. Do not apply if `None` or empty.
        Default: `None`.
    exclude : collection of str, optional
        Set of IDs of the RNAs to ignore. Do not apply if `None` or
        empty. Default: `None`.

    Returns
    -------
    numpy.ndarray of float
        The structure MFEs obtained from the loops.
    """
    return extract_2d_structure_loops_from_list(
        loops_dir, iter_seq_ss(seq_ss_dir, keep=keep, exclude=exclude)
    )


def gen_loops_dict(
    dict_file: StrPath,
    loops_dir: StrPath,
    records: Collection[SeqSSRecord],
) -> None:
    """Generate a dictionary of loops from a 2D structure set.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    loops_dir : str or os.PathLike
        Path of the directory to save the loops files.
    records : collection of SeqSSRecord
        Set containing all the 2D structure records to save into the
        dictionary.
    """
    send_msg("Generating loops dictionary\n")
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
    desc = "  Calculating loops"
    records = track_progress(records, desc, total, "ss", initial)

    extract_2d_structure_loops_from_list(loops_dir, records)
    id_list = set(rec.id for rec in records)
    loops_list = list(
        read_2d_structure_loops_from_dir(loops_dir, keep=id_list)
    )
    keys = (counter.rna for counter in loops_list)
    values = (counter.loops_vectorized for counter in loops_list)
    save_dict(dict_file, keys=keys, values=values)

    send_msg("  Loops dictionary generated\n")


def gen_loops_dict_from_dir(
    dict_file: StrPath,
    loops_dir: StrPath,
    *,
    keep: Collection[str] | None = None,
    exclude: Collection[str] | None = None,
) -> None:
    """Generate a dictionary of loops from loops files.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    loops_dir : str or os.PathLike
        Path of the directory containing the loops files.
    keep : collection of str, optional
        Set of IDs of the RNAs to save. Do not apply if `None` or empty.
        Default: `None`.
    exclude : collection of str, optional
        Set of IDs of the RNAs to ignore. Do not apply if `None` or
        empty. Default: `None`.
    """
    dict_file = Path(dict_file).resolve()
    if dict_file.exists():
        keys = load_dict(dict_file).keys()
        if exclude is None:
            exclude = keys
        else:
            exclude = set(exclude) | keys

    loops_list = list(
        read_2d_structure_loops_from_dir(loops_dir, keep=keep, exclude=exclude)
    )
    values = (counter.loops_vectorized for counter in loops_list)
    if keep:
        keys = sorted(keep)
    else:
        keys = (loop.rna for loop in loops_list)
    save_dict(dict_file, keys=keys, values=values)
