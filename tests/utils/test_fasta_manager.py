import pytest

from compgenlnc.config.paths import SEQ_SS_PREFIX
from compgenlnc.structs import SeqSSRecord
from compgenlnc.utils.fasta_manager import (
    iter_seq_fasta,
    iter_seq_ss,
    load_fasta,
    read_sequence,
    read_structure,
    save_fasta,
)

from constants import SEQ_FASTA_EXAMPLE, SEQ_SS_FOLDER
from params import seq_list, T_converted_seq_list


@pytest.mark.parametrize(
    'index, expected',
    list(enumerate(seq_list))
)
def test_load_fasta(index, expected, seq_fasta_loader):
    record = seq_fasta_loader[index]
    assert record.id == expected.id
    assert record.seq == expected.seq


def test_save_fasta(tmp_path):
    new_fasta = tmp_path / 'fasta.fa'
    molecules = load_fasta(SEQ_FASTA_EXAMPLE, mode='seq')
    save_fasta(new_fasta, molecules, mode='seq')
    assert new_fasta.exists()

    it = load_fasta(SEQ_FASTA_EXAMPLE, mode='seq')
    new_it = load_fasta(new_fasta, mode='seq')
    assert all(a == b for a, b in zip(it, new_it))


@pytest.mark.parametrize(
    'index, expected',
    list(enumerate(seq_list))
)
def test_iter_seq_fasta(index, expected, seq_fasta_list):
    record = seq_fasta_list[index]
    assert record.id == expected.id
    assert record.seq == expected.seq


@pytest.mark.parametrize(
    'filter, expected',
    [
        (['7a'], 3),
        (['hsa'], 10),
        (['-1', '-2'], 4),
        (['non'], 0),
    ]
)
def test_iter_seq_fasta_filtered(filter, expected):
    records_list = [
        record
        for record in iter_seq_fasta(SEQ_FASTA_EXAMPLE, filter)
    ]
    assert len(records_list) == expected


@pytest.mark.parametrize(
    'index, expected',
    list(enumerate(T_converted_seq_list))
)
def test_iter_seq_ss(index, expected, seq_ss_list):
    record = seq_ss_list[index]
    assert record.id == expected.id
    assert record.seq == expected.seq


@pytest.mark.parametrize(
    'keep, expected',
    [
        (['hsa-let-7a-1', 'hsa-let-7a-2', 'hsa-let-7b'], 3),
        (['hsa-let-7f-1', 'hsa-let-7f-2'], 2),
        ([], 10),
        (['non1'], 0),
    ]
)
def test_iter_seq_ss_filtered(keep, expected):
    records_list = [
        record
        for record in iter_seq_ss(SEQ_SS_FOLDER, keep=keep)
    ]
    assert len(records_list) == expected


@pytest.mark.parametrize(
    'record',
    T_converted_seq_list
)
def test_read_sequence(record):
    filename = SEQ_SS_FOLDER / f'{SEQ_SS_PREFIX}{record.id}.dat'
    assert filename.exists()

    seq = read_sequence(filename)
    assert seq == record.seq


def test_read_structure(ss_record_identified):
    filename = SEQ_SS_FOLDER / f'{SEQ_SS_PREFIX}{ss_record_identified.id}.dat'
    assert filename.exists()

    ss = read_structure(filename)
    assert ss == ss_record_identified.ss
