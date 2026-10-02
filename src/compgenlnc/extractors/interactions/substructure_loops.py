import os
from collections.abc import Iterable, Iterator
from pathlib import Path

import numpy as np

from compgenlnc.consts.fields import (
    HIGH_POSITION,
    IS_POSITIVE_FIELD,
    LNCRNA_FIELD,
    LOW_POSITION,
    MIRNA_FIELD,
    PAIR_FIELDS,
)
from compgenlnc.consts.params import LNCRNA, MIRNA
from compgenlnc.config.paths import BINDING_PREFIX, LOOPS_PREFIX
from compgenlnc.fileman.csv_manager import load_interactions
from compgenlnc.fileman.dict_manager import save_dict
from compgenlnc.fileman.loops_manager import get_2d_struture_loops
from compgenlnc.fileman.miranda_manager import get_binding_zone
from compgenlnc.structs.binding_zone import BindingZone
from compgenlnc.structs.loop_counter import LoopCounter
from compgenlnc.typing.numpy_dtypes import InteractionTuple


def get_substructure_loops(
    loop_counter: LoopCounter, binding_zone: BindingZone, molecule: str
) -> LoopCounter:
    if not binding_zone:
        return LoopCounter([], 0)
    loops_arr = loop_counter.loops
    if molecule == LNCRNA:
        start, end = binding_zone.lnc_start, binding_zone.lnc_end
    if molecule == MIRNA:
        start, end = binding_zone.mir_start, binding_zone.mir_end
    mask = (loops_arr[LOW_POSITION] >= start) & (
        loops_arr[HIGH_POSITION] <= end
    )
    new_loops = loops_arr[mask]
    return LoopCounter(new_loops)


def get_substructure_loops_from_list(
    loops_it: Iterable[LoopCounter],
    binding_it: Iterable[BindingZone],
    molecule: str,
) -> Iterator[LoopCounter]:
    for loop_counter, binding_zone in zip(loops_it, binding_it):
        yield get_substructure_loops(loop_counter, binding_zone, molecule)


def gen_substructure_loops_dict(
    dict_file: str | os.PathLike,
    interactions: Iterable[InteractionTuple],
    loops_it: Iterable[LoopCounter],
    binding_it: Iterable[BindingZone],
    molecule: str,
) -> None:
    interactions = list(interactions)
    new_loops = get_substructure_loops_from_list(
        loops_it, binding_it, molecule
    )
    keys = (
        f"{interaction[LNCRNA_FIELD]}_{interaction[MIRNA_FIELD]}"
        for interaction in interactions
    )
    values = (
        loop_counter.loops_vectorized
        if interaction[IS_POSITIVE_FIELD]
        else np.array([-1] * 3 + [0] * 4)
        for loop_counter, interaction in zip(new_loops, interactions)
    )
    save_dict(dict_file, keys=keys, values=values)


# Fix when precursors are mapped to mature miRNAs
def gen_substructure_loops_dict_from_files(
    dict_file: str | os.PathLike,
    interactions_file: str | os.PathLike,
    loops_folder: str | os.PathLike,
    binding_folder: str | os.PathLike,
    molecule: str,
) -> None:
    loops_folder = Path(loops_folder).resolve()
    binding_folder = Path(binding_folder).resolve()

    interactions = load_interactions(interactions_file)
    interactions.sort(order=PAIR_FIELDS)
    keys = [f"{lnc}_{mir}" for lnc, mir in interactions[PAIR_FIELDS]]

    def iter_interactions_loops():
        for key in keys:
            loops_filename = loops_folder / f"{LOOPS_PREFIX}{key}.dat"
            binding_filename = binding_folder / f"{BINDING_PREFIX}{key}.dat"
            if not loops_filename.exists() or not binding_filename.exists():
                continue
            binding_zone = get_binding_zone(binding_filename)
            if not binding_zone:
                yield [-1] * 3 + [0] * 4
            loop_counter = get_2d_struture_loops(loops_filename)
            yield get_substructure_loops(
                loop_counter, binding_zone, molecule
            ).loops_vectorized

    save_dict(dict_file, keys=keys, values=iter_interactions_loops())
