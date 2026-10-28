import numpy as np
from constants import SEQ_FASTA_EXAMPLE, SEQ_SS_FOLDER

from compgenlnc.extractors.sequences import (
    gen_ctd_dict,
    gen_doc2vec_dict,
    gen_kmer_dict,
    gen_role2vec_dict,
    gen_role2vec_dict_from_files,
    gen_role2vec_dict_from_folder,
    get_ctd,
    get_doc2vec,
    get_kmer,
)
from compgenlnc.fileman import load_dict, load_fasta


def test_gen_role2vec_dict(doc2vec_model_deterministic, tmp_path):
    sequences = tuple(load_fasta(SEQ_FASTA_EXAMPLE, mode="seq"))
    doc2vec_model = doc2vec_model_deterministic

    kmer_dict = {record.id: get_kmer(record.seq, 3) for record in sequences}
    ctd_dict = {record.id: get_ctd(record.seq) for record in sequences}
    doc2vec_dict = {
        record.id: get_doc2vec(record.seq, model=doc2vec_model)
        for record in sequences
    }

    filename = tmp_path / "role2vec.dict"
    gen_role2vec_dict(kmer_dict, ctd_dict, doc2vec_dict, filename)
    assert filename.exists()

    role2vec_dict = load_dict(filename)
    assert all(
        len(vector) == 256 and np.issubdtype(vector.dtype, np.float32)
        for vector in role2vec_dict.values()
    )


def test_gen_role2vec_dict_from_files(doc2vec_model_deterministic, tmp_path):
    sequences = tuple(load_fasta(SEQ_FASTA_EXAMPLE, mode="seq"))
    doc2vec_model = doc2vec_model_deterministic

    kmer_file = tmp_path / "kmer.dict"
    ctd_file = tmp_path / "ctd.dict"
    doc2vec_file = tmp_path / "doc2vec.dict"

    gen_kmer_dict(sequences, kmer_file, 3)
    gen_ctd_dict(sequences, ctd_file)
    gen_doc2vec_dict(sequences, doc2vec_file, model=doc2vec_model)

    filename = tmp_path / "role2vec.dict"
    gen_role2vec_dict_from_files(kmer_file, ctd_file, doc2vec_file, filename)
    assert filename.exists()

    role2vec_dict = load_dict(filename)
    assert all(
        len(vector) == 256 and np.issubdtype(vector.dtype, np.float32)
        for vector in role2vec_dict.values()
    )


def test_gen_role2vec_dict_from_folder(tmp_path):
    filename = tmp_path / "role2vec.dict"
    gen_role2vec_dict_from_folder(SEQ_SS_FOLDER, filename, 3)
    assert filename.exists()

    role2vec_dict = load_dict(filename)
    assert all(
        len(vector) == 256 and np.issubdtype(vector.dtype, np.float32)
        for vector in role2vec_dict.values()
    )
