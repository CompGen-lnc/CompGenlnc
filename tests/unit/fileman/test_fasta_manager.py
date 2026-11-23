import pytest

import compgenlnc as cgl
import compgenlnc.fileman.fasta_manager as fasta_manager


def test_load_fasta():
    records = cgl.iter_fasta(cgl.paths.PRE_MIRNA_FASTA, mode=cgl.SEQ_MODE)
    for record in records:
        assert isinstance(record, cgl.structs.SeqSSRecord)
        assert len(record.seq) > len(record.ss) or len(record.seq) == 0

    records = cgl.iter_fasta(cgl.paths.PRE_MIRNA_FASTA, mode=cgl.SS_MODE)
    for record in records:
        assert isinstance(record, cgl.structs.SeqSSRecord)
        assert len(record.ss) > len(record.seq) or len(record.ss) == 0


def test_save_fasta(tmp_path):
    new_fasta = tmp_path / "fasta.fa"
    molecules = cgl.iter_fasta(cgl.paths.PRE_MIRNA_FASTA, mode="seq")
    cgl.save_fasta(new_fasta, molecules, mode="seq")
    assert new_fasta.exists()

    it = cgl.iter_fasta(cgl.paths.PRE_MIRNA_FASTA, mode="seq")
    new_it = cgl.iter_fasta(new_fasta, mode="seq")
    assert all(a == b for a, b in zip(it, new_it))


def test_iter_seq_fasta():
    records = fasta_manager.iter_seq_fasta(cgl.paths.PRE_MIRNA_FASTA)
    for record in records:
        assert isinstance(record, cgl.structs.SeqSSRecord)
        assert len(record.seq) > len(record.ss) or len(record.seq) == 0


@pytest.mark.parametrize(
    "filter, expected",
    [
        (["7a"], 3),
        (["hsa"], 387),
        (["-1", "-2"], 127),
        (["non"], 0),
    ],
)
def test_iter_seq_fasta_filtered(filter, expected):
    records_list = [
        record
        for record in fasta_manager.iter_seq_fasta(
            cgl.paths.PRE_MIRNA_FASTA, filter
        )
    ]
    assert len(records_list) == expected


def test_iter_seq_ss():
    records = cgl.iter_seq_ss(cgl.paths.MIRNA_SEQ_SS_DIR)
    for record in records:
        assert isinstance(record, cgl.structs.SeqSSRecord)
        assert len(record.seq) == len(record.ss)


@pytest.mark.parametrize(
    "keep, expected",
    [
        (["hsa-let-7a-1", "hsa-let-7a-2", "hsa-let-7b"], 2),
        (["hsa-let-7f-1", "hsa-let-7f-2"], 2),
        ([], 387),
        (["non1"], 0),
    ],
)
def test_iter_seq_ss_filtered(keep, expected):
    records = cgl.load_seq_ss(cgl.paths.MIRNA_SEQ_SS_DIR, keep=keep)
    assert len(records) == expected


def test_read_sequence():
    records = cgl.iter_seq_ss(cgl.paths.MIRNA_SEQ_SS_DIR)
    for record in records:
        seq = cgl.read_sequence(cgl.paths.MIRNA_SEQ_SS_DIR, record.id)
        assert str(record.seq) == seq


def test_read_structure():
    records = cgl.iter_seq_ss(cgl.paths.MIRNA_SEQ_SS_DIR)
    for record in records:
        ss = cgl.read_structure(cgl.paths.MIRNA_SEQ_SS_DIR, record.id)
        assert str(record.ss) == ss
