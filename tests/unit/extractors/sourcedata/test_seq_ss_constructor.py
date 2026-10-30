import compgenlnc as cgl


def test_fasta_to_seq_ss(tmp_path):
    folder = tmp_path / "seq+ss_join"
    total = cgl.fasta_to_seq_ss(
        cgl.paths.PRE_MIRNA_FASTA, cgl.paths.MIRNA_SS_FASTA, folder
    )
    assert total == 387


def test_join_seq_ss(tmp_path):
    records = cgl.iter_seq_fasta(cgl.paths.PRE_MIRNA_FASTA)

    for record in records:
        seq = str(record.seq)
        ss = ""
        id_ = record.id
        folder = tmp_path / "seq+ss_join"
        folder.mkdir(parents=True, exist_ok=True)
        cgl.join_seq_ss(seq, ss, id_, folder)

        path = folder / f"{cgl.paths.SEQ_SS_PREFIX}{id_}.dat"
        assert path.exists()

        new_seq = cgl.read_sequence(folder, id_)
        expected = cgl.read_sequence(cgl.paths.MIRNA_SEQ_SS_FOLDER, id_)
        assert new_seq == expected
