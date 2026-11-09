import numpy as np

import compgenlnc as cgl


def test_extract_2d_structure_loops(seq_ss_list, tmp_path):
    for record in seq_ss_list:
        cgl.extract_2d_structure_loops(tmp_path, record)
        filename = tmp_path / f"{cgl.paths.LOOPS_PREFIX}{record.id}.dat"
        assert filename.exists()


def test_extract_2d_structure_loops_from_list(seq_ss_list, tmp_path):
    cgl.extract_2d_structure_loops_from_list(tmp_path, seq_ss_list)
    for record in seq_ss_list:
        listed_file = tmp_path / f"{cgl.paths.LOOPS_PREFIX}{record.id}.dat"
        assert listed_file.exists()
        expected_file = (
            cgl.paths.MIRNA_LOOPS_DIR
            / f"{cgl.paths.LOOPS_PREFIX}{record.id}.dat"
        )
        assert listed_file.read_text() == expected_file.read_text()


def test_extract_2d_structure_loops_from_dir(example_id_filter, tmp_path):
    cgl.extract_2d_structure_loops_from_dir(
        tmp_path, cgl.paths.MIRNA_SEQ_SS_DIR, keep=example_id_filter
    )
    for record in cgl.iter_seq_ss(
        cgl.paths.MIRNA_SEQ_SS_DIR, keep=example_id_filter
    ):
        listed_file = tmp_path / f"{cgl.paths.LOOPS_PREFIX}{record.id}.dat"
        assert listed_file.exists()
        expected_file = (
            cgl.paths.MIRNA_LOOPS_DIR
            / f"{cgl.paths.LOOPS_PREFIX}{record.id}.dat"
        )
        assert listed_file.read_text() == expected_file.read_text()


def test_gen_loops_dict(seq_ss_list, tmp_path):
    dict_path = tmp_path / "loops.dict"
    cgl.gen_loops_dict(dict_path, cgl.paths.MIRNA_LOOPS_DIR, seq_ss_list)
    loops_dict = cgl.load_dict(dict_path)
    for key, value in loops_dict.items():
        expected = cgl.read_2d_structure_loops(
            cgl.paths.MIRNA_LOOPS_DIR, key
        ).loops_vectorized
        assert np.array_equal(value, expected)


def test_gen_loops_dict_from_dir(seq_ss_list, tmp_path):
    dict_path = tmp_path / "loops.dict"
    expected_path = tmp_path / "expected.dict"
    cgl.gen_loops_dict_from_dir(dict_path, cgl.paths.MIRNA_LOOPS_DIR)
    cgl.gen_loops_dict(expected_path, cgl.paths.MIRNA_LOOPS_DIR, seq_ss_list)
    loops_dict, expected_dict = (
        cgl.load_dict(dict_path),
        cgl.load_dict(expected_path),
    )
    for key in expected_dict.keys():
        assert np.array_equal(loops_dict[key], expected_dict[key])
