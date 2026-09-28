import os
from collections.abc import Iterable, Iterator
from pathlib import Path

from compgenlnc.config.paths import BINDING_PREFIX, LOOPS_PREFIX
from compgenlnc.consts.params import LNCRNA_TYPE, MIRNA_TYPE
from compgenlnc.structs.loop_counter import LoopCounter
from compgenlnc.structs.binding_zone import BindingZone
from compgenlnc.typing.numpy_dtypes import InteractionTuple
from compgenlnc.utils.dict_manager import save_dict
from compgenlnc.utils.interactions_manager import load_interactions
from compgenlnc.utils.loops_manager import get_2d_struture_loops
from compgenlnc.utils.miranda_manager import get_binding_zone


def get_substructure_loops(
    loop_counter: LoopCounter, binding_zone: BindingZone, molecule: str
) -> LoopCounter:
    if not binding_zone:
        return LoopCounter([], 0)
    loops_arr = loop_counter.loops
    if molecule == LNCRNA_TYPE:
        start, end = binding_zone.lnc_start, binding_zone.lnc_end
    if molecule == MIRNA_TYPE:
        start, end = binding_zone.mir_start, binding_zone.mir_end
    mask = loops_arr["low"] >= start & loops_arr["high"] <= end
    new_loops = loops_arr[mask]
    return LoopCounter(new_loops)


def get_substructure_loops_from_list(
    loops_list: Iterable[LoopCounter],
    binding_list: Iterable[BindingZone],
    molecule: str,
) -> Iterator[LoopCounter]:
    for loop_counter, binding_zone in zip(loops_list, binding_list):
        yield get_substructure_loops(loop_counter, binding_zone, molecule)


def gen_substructure_loops_dict(
    dict_file: str | os.PathLike,
    interactions: Iterable[InteractionTuple],
    loops_list: Iterable[LoopCounter],
    binding_list: Iterable[BindingZone],
    molecule: str,
) -> None:
    dict_file = Path(dict_file).resolve()
    dict_file.mkdir(parents=True, exist_ok=True)

    new_loops = get_substructure_loops_from_list(
        loops_list, binding_list, molecule
    )
    keys = (
        f"{interaction['lncRNA']}_{interaction['miRNA']}"
        for interaction in interactions
    )
    values = (
        loop_counter.loops_vectorized if binding_zone else [-1] * 3 + [0] * 4
        for loop_counter, binding_zone in zip(new_loops, binding_list)
    )
    save_dict(dict_file, keys=keys, values=values)


def gen_substructure_from_files(
    dict_file: str | os.PathLike,
    interactions_file: str | os.PathLike,
    loops_folder: str | os.PathLike,
    binding_folder: str | os.PathLike,
    molecule: str,
) -> None:
    loops_folder = Path(loops_folder).resolve()
    binding_folder = Path(binding_folder).resolve()

    interactions = load_interactions(interactions_file)
    interactions.sort(order=["lncRNA", "miRNA"])
    keys = [f"{lnc}_{mir}" for lnc, mir in interactions[["lncRNA", "miRNA"]]]

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
