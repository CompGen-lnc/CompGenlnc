import numpy as np

from compgenlnc.config.paths import LOOPS_PREFIX
from compgenlnc.utils import iter_seq_ss
from compgenlnc.generators.structures2d import (
    count_2d_structure_loops,
    count_2d_structure_loops_from_folder,
    gen_2d_structure_loops,
    gen_2d_structure_loops_from_folder,
    gen_2d_structure_loops_from_list,
    get_2d_struture_loops,
    get_2d_structure_loops_from_folder,
)

from constants import SEQ_SS_FOLDER, SS_LOOPS_FOLDER


def test_gen_2d_structure_loops(seq_ss_record, tmp_path):
    gen_2d_structure_loops(seq_ss_record, tmp_path)
    filename = tmp_path / f'{LOOPS_PREFIX}{seq_ss_record.id}.dat'
    assert filename.exists()


def test_gen_2d_structure_loops_from_list(seq_ss_filtered_list, tmp_path):
    listed_energy = gen_2d_structure_loops_from_list(
        seq_ss_filtered_list, tmp_path
    )
    for record in seq_ss_filtered_list:
        listed_file = tmp_path / f'{LOOPS_PREFIX}{record.id}.dat'
        assert listed_file.exists()
        expected_file = SS_LOOPS_FOLDER / f'{LOOPS_PREFIX}{record.id}.dat'
        assert listed_file.read_text() == expected_file.read_text()


def test_gen_2d_structure_loops_from_folder(mirna_id_filter, tmp_path):
    listed_energy = gen_2d_structure_loops_from_folder(
        SEQ_SS_FOLDER, tmp_path, keep=mirna_id_filter
    )
    for record in iter_seq_ss(SEQ_SS_FOLDER, keep=mirna_id_filter):
        listed_file = tmp_path / f'{LOOPS_PREFIX}{record.id}.dat'
        assert listed_file.exists()
        expected_file = SS_LOOPS_FOLDER / f'{LOOPS_PREFIX}{record.id}.dat'
        assert listed_file.read_text() == expected_file.read_text()


def test_get_2d_structure_loops(seq_ss_record):
    loops_file = SS_LOOPS_FOLDER / f'{LOOPS_PREFIX}{seq_ss_record.id}.dat'
    energy, loops = get_2d_struture_loops(loops_file)
    total_energy = loops['energy'].sum().item() / 100
    assert energy == total_energy


def test_get_2d_structure_loops_from_folder(mirna_id_filter):
    loops_iter = get_2d_structure_loops_from_folder(
        SS_LOOPS_FOLDER, keep=mirna_id_filter
    )
    for id_, structure_loops in zip(mirna_id_filter, loops_iter):
        energy, loops = structure_loops
        loops_file = SS_LOOPS_FOLDER / f'{LOOPS_PREFIX}{id_}.dat'
        expectd_energy, expected_loops = get_2d_struture_loops(loops_file)
        assert energy == expectd_energy
        assert np.array_equal(loops, expected_loops)


def test_count_2d_structure_loops(seq_ss_record):
    loops_file = SS_LOOPS_FOLDER / f'{LOOPS_PREFIX}{seq_ss_record.id}.dat'
    (
        total_energy,
        energy_dict,
        count_dict,
    ) = count_2d_structure_loops(loops_file)
    expected_energy, loops = get_2d_struture_loops(loops_file)
    assert total_energy == expected_energy

    for kind in energy_dict.keys():
        kind_list = loops[loops['kind'] == kind]
        kind_energy = kind_list['energy'].sum().item() / 100
        kind_count = len(kind_list)
        assert kind_energy == energy_dict[kind]
        assert kind_count == count_dict[kind]


def test_count_2d_structure_loops_from_folder(mirna_id_filter):
    count_iter = count_2d_structure_loops_from_folder(
        SEQ_SS_FOLDER, keep=mirna_id_filter
    )
    for id_, count_loops in zip(mirna_id_filter, count_iter):
        total_energy, energy_dict, count_dict = count_loops
        loops_file = SS_LOOPS_FOLDER / f'{LOOPS_PREFIX}{id_}.dat'
        (
            expected_total_energy,
            expected_energy_dict,
            expected_count_dict,
        ) = count_2d_structure_loops(loops_file)
        assert total_energy == expected_total_energy
        assert energy_dict == expected_energy_dict
        assert count_dict == expected_count_dict
