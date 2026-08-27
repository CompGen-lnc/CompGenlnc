import numpy as np
import pytest

from compgenlnc.config.paths import SEQ_SS_PREFIX
from compgenlnc.generators.sequences.ctd import (
    gen_ctd_dict_from_folder, get_ctd
)
from compgenlnc.generators.sequences.kmer import (
    gen_kmer_dict_from_folder, get_kmer
)
from compgenlnc.structs import SeqSSRecord
from compgenlnc.utils.fasta_manager import load_fasta
from compgenlnc.utils.dict_manager import (
    load_dict,
    save_dict,
)
from constants import SEQ_FASTA_EXAMPLE, SEQ_SS_FOLDER


@pytest.mark.parametrize(
    'dict_gen, args, type_, expected',
    [
        (gen_ctd_dict_from_folder, [], np.float64, {
            record.id: get_ctd(record.seq)
            for record in load_fasta(SEQ_FASTA_EXAMPLE, mode='seq')
        }),
        (gen_kmer_dict_from_folder, [3], np.float64, {
            record.id: get_kmer(record.seq, 3)
            for record in load_fasta(SEQ_FASTA_EXAMPLE, mode='seq')
        }),
        (gen_kmer_dict_from_folder, [4], np.float64, {
            record.id: get_kmer(record.seq, 4)
            for record in load_fasta(SEQ_FASTA_EXAMPLE, mode='seq')
        }),
    ]
)
def test_load_dict(dict_gen, args, type_, expected, tmp_path):
    filename = tmp_path / 'example.dict'
    dict_gen(SEQ_SS_FOLDER, filename, *args)
    dict_ = load_dict(filename, type_)
    assert set(dict_.keys()) == set(expected.keys())
    assert all(
        np.array_equal(dict_[key], expected[key])
        for key in list(dict_.keys())
    )


@pytest.mark.parametrize(
    'type_, expected',
    [
        (np.float64, {
            record.id: get_ctd(record.seq)
            for record in load_fasta(SEQ_FASTA_EXAMPLE, mode='seq')
        }),
        (np.float64, {
            record.id: get_kmer(record.seq, 3)
            for record in load_fasta(SEQ_FASTA_EXAMPLE, mode='seq')
        }),
        (np.float64, {
            record.id: get_kmer(record.seq, 4)
            for record in load_fasta(SEQ_FASTA_EXAMPLE, mode='seq')
        }),
    ]
)
def test_save_dict(type_, expected, tmp_path):
    new_dict = tmp_path / 'dict.fa'
    save_dict(new_dict, expected)
    assert new_dict.exists()

    dict_= load_dict(new_dict, type_)
    assert set(dict_.keys()) == set(expected.keys())
    assert all(
        np.array_equal(dict_[key], expected[key])
        for key in list(dict_.keys())
    )