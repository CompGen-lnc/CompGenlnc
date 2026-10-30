import pytest

import compgenlnc as cgl


@pytest.mark.parametrize(
    "filter, expected",
    [
        (["7a"], ["hsa-let-7a-1", "hsa-mir-4477a", "hsa-mir-517a"]),
        (["7f", "7c"], ["hsa-let-7c", "hsa-let-7f-1", "hsa-let-7f-2"]),
        (
            ["hsa"],
            [
                record.id
                for record in cgl.load_fasta(
                    cgl.paths.PRE_MIRNA_FASTA, mode="seq"
                )
            ],
        ),
        (["non"], []),
    ],
)
def test_filter_seq_fasta(filter, expected, tmp_path):
    new_fasta = tmp_path / "new_fasta.fa"
    cgl.filter_fasta(cgl.paths.PRE_MIRNA_FASTA, new_fasta, filter)
    assert new_fasta.exists()

    ids = [record.id for record in cgl.load_fasta(new_fasta, mode="seq")]
    assert ids == expected


def test_gen_fasta_2d(example_id_filter, tmp_path):
    ss_fasta = tmp_path / "fasta2d.fa"
    cgl.gen_fasta_2d(
        cgl.paths.PRE_MIRNA_FASTA, ss_fasta, keep=example_id_filter
    )
    seq_list = sorted(
        cgl.load_fasta(
            cgl.paths.PRE_MIRNA_FASTA, mode="seq", keep=example_id_filter
        ),
        key=lambda x: x.id,
    )
    ss_list = sorted(cgl.load_fasta(ss_fasta, mode="ss"), key=lambda x: x.id)

    for seq, expected in zip(seq_list, ss_list):
        ss = cgl.extract_2d_structure(seq)
        assert ss == expected.ss


def test_seq_ss_to_fasta(tmp_path):
    new_seq_fasta_file = tmp_path / "new_seq_fasta.fa"
    new_ss_fasta_file = tmp_path / "new_ss_fasta.fa"
    cgl.seq_ss_to_fasta(
        cgl.paths.MIRNA_SEQ_SS_FOLDER, new_seq_fasta_file, new_ss_fasta_file
    )
    assert new_seq_fasta_file.exists()

    new_seq_fasta = sorted(
        cgl.load_fasta(new_seq_fasta_file), key=lambda x: x.id
    )
    seq_fasta = sorted(
        cgl.load_fasta(cgl.paths.PRE_MIRNA_FASTA), key=lambda x: x.id
    )
    assert new_seq_fasta == seq_fasta
    new_ss_fasta = sorted(
        cgl.load_fasta(new_ss_fasta_file), key=lambda x: x.id
    )
    ss_fasta = sorted(
        cgl.load_fasta(cgl.paths.MIRNA_SS_FASTA), key=lambda x: x.id
    )
    assert new_ss_fasta == ss_fasta
