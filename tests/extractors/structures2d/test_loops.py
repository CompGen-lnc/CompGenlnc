from constants import SEQ_SS_FOLDER, SS_LOOPS_FOLDER

from compgenlnc.config.paths import LOOPS_PREFIX
from compgenlnc.extractors.structures2d import (
    gen_2d_structure_loops,
    gen_2d_structure_loops_from_folder,
    gen_2d_structure_loops_from_list,
)
from compgenlnc.utils import iter_seq_ss


def test_gen_2d_structure_loops(seq_ss_record, tmp_path):
    gen_2d_structure_loops(seq_ss_record, tmp_path)
    filename = tmp_path / f"{LOOPS_PREFIX}{seq_ss_record.id}.dat"
    assert filename.exists()


def test_gen_2d_structure_loops_from_list(seq_ss_filtered_list, tmp_path):
    gen_2d_structure_loops_from_list(seq_ss_filtered_list, tmp_path)
    for record in seq_ss_filtered_list:
        listed_file = tmp_path / f"{LOOPS_PREFIX}{record.id}.dat"
        assert listed_file.exists()
        expected_file = SS_LOOPS_FOLDER / f"{LOOPS_PREFIX}{record.id}.dat"
        assert listed_file.read_text() == expected_file.read_text()


def test_gen_2d_structure_loops_from_folder(example_id_filter, tmp_path):
    gen_2d_structure_loops_from_folder(
        SEQ_SS_FOLDER, tmp_path, keep=example_id_filter
    )
    for record in iter_seq_ss(SEQ_SS_FOLDER, keep=example_id_filter):
        listed_file = tmp_path / f"{LOOPS_PREFIX}{record.id}.dat"
        assert listed_file.exists()
        expected_file = SS_LOOPS_FOLDER / f"{LOOPS_PREFIX}{record.id}.dat"
        assert listed_file.read_text() == expected_file.read_text()
