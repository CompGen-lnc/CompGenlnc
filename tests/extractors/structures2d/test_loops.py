import numpy as np

from compgenlnc.config.paths import LOOPS_PREFIX
from compgenlnc.extractors.structures2d import (
    extract_2d_structure_loops,
    extract_2d_structure_loops_from_folder,
    extract_2d_structure_loops_from_list,
    gen_loops_dict,
    gen_loops_dict_from_folder,
)
from compgenlnc.fileman import get_2d_struture_loops, iter_seq_ss, load_dict

from constants import SEQ_SS_FOLDER, SS_LOOPS_FOLDER


def test_extract_2d_structure_loops(seq_ss_record, tmp_path):
    extract_2d_structure_loops(seq_ss_record, tmp_path)
    filename = tmp_path / f"{LOOPS_PREFIX}{seq_ss_record.id}.dat"
    assert filename.exists()


def test_extract_2d_structure_loops_from_list(seq_ss_filtered_list, tmp_path):
    extract_2d_structure_loops_from_list(seq_ss_filtered_list, tmp_path)
    for record in seq_ss_filtered_list:
        listed_file = tmp_path / f"{LOOPS_PREFIX}{record.id}.dat"
        assert listed_file.exists()
        expected_file = SS_LOOPS_FOLDER / f"{LOOPS_PREFIX}{record.id}.dat"
        assert listed_file.read_text() == expected_file.read_text()


def test_extract_2d_structure_loops_from_folder(example_id_filter, tmp_path):
    extract_2d_structure_loops_from_folder(
        SEQ_SS_FOLDER, tmp_path, keep=example_id_filter
    )
    for record in iter_seq_ss(SEQ_SS_FOLDER, keep=example_id_filter):
        listed_file = tmp_path / f"{LOOPS_PREFIX}{record.id}.dat"
        assert listed_file.exists()
        expected_file = SS_LOOPS_FOLDER / f"{LOOPS_PREFIX}{record.id}.dat"
        assert listed_file.read_text() == expected_file.read_text()


def test_gen_loops_dict(seq_ss_list, tmp_path):
    dict_path = tmp_path / "loops.dict"
    gen_loops_dict(dict_path, seq_ss_list, SS_LOOPS_FOLDER)
    loops_dict = load_dict(dict_path)
    for key, value in loops_dict.items():
        expected_file = SS_LOOPS_FOLDER / f"{LOOPS_PREFIX}{key}.dat"
        expected = get_2d_struture_loops(expected_file).loops_vectorized
        assert np.array_equal(value, expected)


def test_gen_loops_dict_from_folder(seq_ss_list, tmp_path):
    dict_path = tmp_path / "loops.dict"
    expected_path = tmp_path / "expected.dict"
    gen_loops_dict_from_folder(dict_path, SS_LOOPS_FOLDER)
    gen_loops_dict(expected_path, seq_ss_list, SS_LOOPS_FOLDER)
    loops_dict, expected_dict = load_dict(dict_path), load_dict(expected_path)
    for key in expected_dict.keys():
        assert np.array_equal(loops_dict[key], expected_dict[key])
