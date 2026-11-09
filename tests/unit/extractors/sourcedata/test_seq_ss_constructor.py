import compgenlnc as cgl


def test_fasta_to_seq_ss(tmp_path):
    dir = tmp_path / "seq+ss_join"
    total = cgl.fasta_to_seq_ss(
        cgl.paths.PRE_MIRNA_FASTA, cgl.paths.MIRNA_SS_FASTA, dir
    )
    assert total == 387


def test_join_seq_ss(tmp_path):
    records = cgl.iter_seq_fasta(cgl.paths.PRE_MIRNA_FASTA)

    for record in records:
        seq = str(record.seq)
        ss = ""
        id_ = record.id
        dir = tmp_path / "seq+ss_join"
        dir.mkdir(parents=True, exist_ok=True)
        cgl.join_seq_ss(dir, id_, seq, ss)

        path = dir / f"{cgl.paths.SEQ_SS_PREFIX}{id_}.dat"
        assert path.exists()

        new_seq = cgl.read_sequence(dir, id_)
        expected = cgl.read_sequence(cgl.paths.MIRNA_SEQ_SS_DIR, id_)
        assert new_seq == expected
