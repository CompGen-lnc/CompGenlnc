import os
import sys
from collections.abc import Iterable
from itertools import tee
from pathlib import Path

import numpy as np
import RNA

from compgenlnc.consts.paths import LOOPS_PREFIX
from compgenlnc.fileman.fasta_manager import iter_seq_ss
from compgenlnc.fileman.loops_manager import (
    read_2d_structure_loops,
    read_2d_structure_loops_from_folder,
)
from compgenlnc.fileman.dict_manager import load_dict, save_dict
from compgenlnc.structs.seq_ss_record import SeqSSRecord


def extract_2d_structure_loops(
    record: SeqSSRecord, loops_folder: str | os.PathLike
) -> float:
    loops_folder = Path(loops_folder).resolve()
    if loops_folder.exists():
        files = os.listdir(loops_folder)
        if f"{LOOPS_PREFIX}{record.id}.dat" in files:
            return read_2d_structure_loops(loops_folder, record.id).energy
    else:
        loops_folder.mkdir(parents=True, exist_ok=True)

    filename = loops_folder / f"{LOOPS_PREFIX}{record.id}.dat"
    with open(filename, "w") as out_file:
        try:
            fc = RNA.fold_compound(str(record.seq), None, RNA.OPTION_EVAL_ONLY)
            energy = fc.eval_structure_verbose(str(record.ss), out_file)
        except:
            energy = 0
        out_file.write(f"{'Energy':<40}: {energy:>8.2f}")
    return energy


def extract_2d_structure_loops_from_list(
    record_list: Iterable[SeqSSRecord],
    loops_folder: str | os.PathLike,
) -> np.typing.NDArray[np.float32]:
    energy_list = []
    for record in record_list:
        energy_list.append(extract_2d_structure_loops(record, loops_folder))

    return np.array(energy_list, dtype=np.float32)


def extract_2d_structure_loops_from_folder(
    seq_ss_folder: str | os.PathLike,
    loops_folder: str | os.PathLike,
    *,
    keep: Iterable[str] | None = None,
    exclude: Iterable[str] | None = None,
) -> np.typing.NDArray[np.float32]:
    return extract_2d_structure_loops_from_list(
        iter_seq_ss(seq_ss_folder, keep=keep, exclude=exclude), loops_folder
    )


def gen_loops_dict(
    dict_file: str | os.PathLike,
    record_list: Iterable[SeqSSRecord],
    loops_folder: str | os.PathLike,
) -> None:
    dict_file = Path(dict_file).resolve()
    if dict_file.exists():
        existing_dict = load_dict(dict_file)
        record_list = (
            record for record in record_list if record.id not in existing_dict
        )

    record_list, record_copy = tee(record_list)
    extract_2d_structure_loops_from_list(record_list, loops_folder)
    id_list = set(rec.id for rec in record_copy)
    loops_list = list(
        read_2d_structure_loops_from_folder(loops_folder, keep=id_list)
    )
    keys = (counter.rna for counter in loops_list)
    values = (counter.loops_vectorized for counter in loops_list)
    save_dict(dict_file, keys=keys, values=values)


def gen_loops_dict_from_folder(
    dict_file: str | os.PathLike,
    loops_folder: str | os.PathLike,
    *,
    keep: Iterable[str] | None = None,
    exclude: Iterable[str] | None = None,
) -> None:
    dict_file = Path(dict_file).resolve()
    if dict_file.exists():
        keys = load_dict(dict_file).keys()
        if exclude is None:
            exclude = keys
        else:
            exclude = set(exclude) | keys

    loops_list = list(
        read_2d_structure_loops_from_folder(
            loops_folder, keep=keep, exclude=exclude
        )
    )
    values = (counter.loops_vectorized for counter in loops_list)
    if keep:
        keys = sorted(keep)
    else:
        keys = (loop.rna for loop in loops_list)
    save_dict(dict_file, keys=keys, values=values)
