import numpy as np
from constants import SS_LOOPS_FOLDER

from compgenlnc.consts import ENERGY_FIELD
from compgenlnc.config.paths import LOOPS_PREFIX
from compgenlnc.fileman import (
    get_2d_structure_loops_from_folder,
    get_2d_struture_loops,
)


def test_get_2d_structure_loops(seq_ss_record):
    loops_file = SS_LOOPS_FOLDER / f"{LOOPS_PREFIX}{seq_ss_record.id}.dat"
    loop_counter = get_2d_struture_loops(loops_file)
    total_energy = loop_counter.loops[ENERGY_FIELD].sum().item() / 100
    assert loop_counter.energy == total_energy


def test_get_2d_structure_loops_from_folder(example_id_filter):
    loops_iter = get_2d_structure_loops_from_folder(
        SS_LOOPS_FOLDER, keep=example_id_filter
    )
    for id_, loop_counter in zip(example_id_filter, loops_iter):
        loops_file = SS_LOOPS_FOLDER / f"{LOOPS_PREFIX}{id_}.dat"
        expected = get_2d_struture_loops(loops_file)
        assert loop_counter == expected
