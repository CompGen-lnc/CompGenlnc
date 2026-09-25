import os
from collections.abc import Iterable
from pathlib import Path

import numpy as np
from ViennaRNA import fold_compound

from compgenlnc.config.paths import LOOPS_PREFIX
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.utils.fasta_manager import iter_seq_ss

dtype_loop_tuple = np.dtype(
    [
        ("kind", "U10"),
        ("low", np.int32),
        ("high", np.int32),
        ("energy", np.int32),
    ]
)


def gen_2d_structure_loops(
    record: SeqSSRecord, loops_folder: str | os.PathLike
) -> float:
    loops_folder = Path(loops_folder).resolve()
    loops_folder.mkdir(parents=True, exist_ok=True)

    fc = fold_compound(str(record.seq))

    filename = loops_folder / f"{LOOPS_PREFIX}{record.id}.dat"
    with open(filename, "w") as out_file:
        energy = fc.eval_structure_verbose(str(record.ss), out_file)
        out_file.write(f"{'Energy':<40}: {energy:>8.2f}")
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
    /,
    keep: Iterable[str] | None = None,
) -> np.typing.NDArray[np.float32]:
    return gen_2d_structure_loops_from_list(
        iter_seq_ss(seq_ss_folder, keep=keep), loops_folder
    )
