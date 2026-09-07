import os
import sys
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np
from ViennaRNA import fold_compound

from compgenlnc.config.paths import LOOPS_PREFIX
from compgenlnc.consts.regex import EXTERNAL_LOOP_RE, LOOP_RE
from compgenlnc.structs import SeqSSRecord
from compgenlnc.utils import iter_seq_ss


loop_tuple = np.dtype([
    ('kind', 'U10'), 
    ('low', np.int32), 
    ('high', np.int32), 
    ('energy', np.int32),
])


def gen_2d_structure_loops(
        record: SeqSSRecord,
        loops_folder: str | os.PathLike
) -> float:
    loops_folder = Path(loops_folder).resolve()
    loops_folder.mkdir(parents=True, exist_ok=True)

    fc = fold_compound(str(record.seq))

    filename = loops_folder / f'{LOOPS_PREFIX}{record.id}.dat'
    with open(filename, 'w') as out_file:
        energy = fc.eval_structure_verbose(str(record.ss), out_file)
        out_file.write(f'{'Energy':<40}: {energy:>8.2f}')
    return energy


def gen_2d_structure_loops_from_list(
        record_list: Iterable[SeqSSRecord],
        loops_folder: str | os.PathLike,
) -> np.typing.NDArray[np.float32]:
    loops_folder = Path(loops_folder).resolve()
    loops_folder.mkdir(parents=True, exist_ok=True)
    
    energy_list = []
    for record in record_list:
        energy_list.append(gen_2d_structure_loops(record, loops_folder))

    return np.array(energy_list, dtype=np.float32)


def gen_2d_structure_loops_from_folder(
        seq_ss_folder: str | os.PathLike,
        loops_folder: str | os.PathLike,
        /, *,
        keep: Iterable[str] | None = None
) -> np.typing.NDArray[np.float32]:
    return gen_2d_structure_loops_from_list(
        iter_seq_ss(seq_ss_folder, keep=keep), loops_folder
    )



def get_2d_struture_loops(
        loops_file: str | os.PathLike
) -> tuple[float, np.typing.NDArray[np.void]]:
    def match_loop(line: str) -> tuple[bool, tuple[str, int, int, int]]:
        m = LOOP_RE.match(line)
        if m:
            return True, (
                m.group(1).split()[0],
                *[int(x) for x in m.group(2, 3, 4)]
            )
        m = EXTERNAL_LOOP_RE.match(line)
        if m:
            return True, (m.group(1).split()[0], 0, 0, int(m.group(2)))
        return False, ()

    energy = 0
    loops = []
    with open(loops_file) as file:
        for line in file.readlines():
            if not line.strip():
                continue
            is_loop, loop = match_loop(line)
            words = line.split()
            if is_loop:
                loops.append(loop)
            elif words[0] == 'Energy':
                energy = float(words[-1])
    return energy, np.array(loops, dtype=loop_tuple)


def get_2d_structure_loops_from_folder(
        loops_folder: str | os.PathLike,
        /, *,
        keep: Iterable[str] | None = None,
) -> Iterator[tuple[float, np.typing.NDArray[np.void]]]:
    loops_folder = Path(loops_folder).resolve()
    for filename in os.listdir(loops_folder):
        id_ = filename.split('.')[0][len(LOOPS_PREFIX):]
        if (keep is not None and id_ not in keep or
            not filename.startswith(LOOPS_PREFIX)):
            continue
        filename = loops_folder / filename
        yield get_2d_struture_loops(filename)



def count_2d_structure_loops(
        loops_file: str | os.PathLike
) -> tuple[float, dict[str, int], dict[str, float]]:
    keys = ['Hairpin', 'Interior', 'Multi']
    total_energy = 0
    energy_dict = dict.fromkeys(keys, 0)
    count_dict = dict.fromkeys(keys, 0)

    with open(loops_file) as file:
        for line in file.readlines():
            if not line.strip():
                continue
            words = line.split()
            kind, energy = words[0], float(words[-1])

            if kind == 'Energy':
                total_energy = energy
            elif kind in keys:
                count_dict[kind] += 1
                energy_dict[kind] += energy

    for kind in keys:
        energy_dict[kind] /= 100

    return total_energy, energy_dict, count_dict


def count_2d_structure_loops_from_folder(
        loops_folder: str | os.PathLike,
        /, *,
        keep: Iterable[str] | None = None,
) -> Iterator[tuple[float, dict[str, int], dict[str, float]]]:
    loops_folder = Path(loops_folder).resolve()
    for filename in os.listdir(loops_folder):
        id_ = filename.split('.')[0][len(LOOPS_PREFIX):]
        if (keep is not None and id_ not in keep or
            not filename.startswith(LOOPS_PREFIX)):
            continue
        filename = loops_folder / filename
        yield count_2d_structure_loops(filename)
    