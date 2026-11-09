import numpy as np
import pytest

import compgenlnc as cgl


@pytest.mark.parametrize(
    "dict_gen, args, type_, expected_factory",
    [
        (
            cgl.gen_ctd_dict_from_dir,
            [],
            np.float32,
            lambda: {
                record.id: cgl.get_ctd(record.seq)
                for record in cgl.iter_seq_ss(cgl.paths.MIRNA_SEQ_SS_DIR)
            },
        ),
        (
            cgl.gen_kmer_dict_from_dir,
            [3],
            np.float32,
            lambda: {
                record.id: cgl.get_kmer(record.seq, 3)
                for record in cgl.iter_seq_ss(cgl.paths.MIRNA_SEQ_SS_DIR)
            },
        ),
        (
            cgl.gen_kmer_dict_from_dir,
            [4],
            np.float32,
            lambda: {
                record.id: cgl.get_kmer(record.seq, 4)
                for record in cgl.iter_seq_ss(cgl.paths.MIRNA_SEQ_SS_DIR)
            },
        ),
    ],
)
def test_load_dict(dict_gen, args, type_, expected_factory, tmp_path):
    expected = expected_factory()
    filename = tmp_path / "example.dict"
    dict_gen(filename, cgl.paths.MIRNA_SEQ_SS_DIR, *args)
    assert filename.exists()
    dict_ = cgl.load_dict(filename, type_)
    assert set(dict_.keys()) == set(expected.keys())
    assert all(
        np.array_equal(dict_[key], expected[key]) for key in list(dict_.keys())
    )


@pytest.mark.parametrize(
    "type_, expected",
    [
        (
            np.float32,
            {
                record.id: cgl.get_ctd(record.seq)
                for record in cgl.iter_fasta(
                    cgl.paths.PRE_MIRNA_FASTA, mode="seq"
                )
            },
        ),
        (
            np.float32,
            {
                record.id: cgl.get_kmer(record.seq, 3)
                for record in cgl.iter_fasta(
                    cgl.paths.PRE_MIRNA_FASTA, mode="seq"
                )
            },
        ),
        (
            np.float32,
            {
                record.id: cgl.get_kmer(record.seq, 4)
                for record in cgl.iter_fasta(
                    cgl.paths.PRE_MIRNA_FASTA, mode="seq"
                )
            },
        ),
    ],
)
def test_save_dict(type_, expected, tmp_path):
    new_dict = tmp_path / "dict.fa"
    cgl.save_dict(new_dict, expected)
    assert new_dict.exists()

    dict_ = cgl.load_dict(new_dict, type_)
    assert set(dict_.keys()) == set(expected.keys())
    assert all(
        np.array_equal(dict_[key], expected[key]) for key in list(dict_.keys())
    )
