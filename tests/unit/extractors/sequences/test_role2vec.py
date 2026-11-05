import numpy as np

import compgenlnc as cgl


def test_gen_role2vec_dict(doc2vec_model_deterministic, tmp_path):
    sequences = tuple(cgl.load_fasta(cgl.paths.PRE_MIRNA_FASTA, mode="seq"))
    doc2vec_model = doc2vec_model_deterministic

    kmer_dict = {
        record.id: cgl.get_kmer(record.seq, 3) for record in sequences
    }
    ctd_dict = {record.id: cgl.get_ctd(record.seq) for record in sequences}
    doc2vec_dict = {
        record.id: cgl.get_doc2vec(record.seq, model=doc2vec_model)
        for record in sequences
    }

    filename = tmp_path / "role2vec.dict"
    cgl.gen_role2vec_dict(filename, kmer_dict, ctd_dict, doc2vec_dict)
    assert filename.exists()

    role2vec_dict = cgl.load_dict(filename)
    assert all(
        len(vector) == 256 and np.issubdtype(vector.dtype, np.float32)
        for vector in role2vec_dict.values()
    )


def test_gen_role2vec_dict_from_files(doc2vec_model_deterministic, tmp_path):
    sequences = tuple(cgl.load_fasta(cgl.paths.PRE_MIRNA_FASTA, mode="seq"))
    doc2vec_model = doc2vec_model_deterministic

    kmer_file = tmp_path / "kmer.dict"
    ctd_file = tmp_path / "ctd.dict"
    doc2vec_file = tmp_path / "doc2vec.dict"

    cgl.gen_kmer_dict(kmer_file, sequences, 3)
    cgl.gen_ctd_dict(ctd_file, sequences)
    cgl.gen_doc2vec_dict(doc2vec_file, sequences, model=doc2vec_model)

    filename = tmp_path / "role2vec.dict"
    cgl.gen_role2vec_dict_from_files(
        filename, kmer_file, ctd_file, doc2vec_file
    )
    assert filename.exists()

    role2vec_dict = cgl.load_dict(filename)
    assert all(
        len(vector) == 256 and np.issubdtype(vector.dtype, np.float32)
        for vector in role2vec_dict.values()
    )


def test_gen_role2vec_dict_from_folder(tmp_path):
    filename = tmp_path / "role2vec.dict"
    cgl.gen_role2vec_dict_from_folder(
        filename, cgl.paths.MIRNA_SEQ_SS_FOLDER, 3
    )
    assert filename.exists()

    role2vec_dict = cgl.load_dict(filename)
    assert all(
        len(vector) == 256 and np.issubdtype(vector.dtype, np.float32)
        for vector in role2vec_dict.values()
    )
