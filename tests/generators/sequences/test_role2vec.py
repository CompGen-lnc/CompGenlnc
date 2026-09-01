import time

import numpy as np

from compgenlnc.generators.sequences.ctd import gen_ctd_dict, get_ctd
from compgenlnc.generators.sequences.doc2vec import (
    gen_doc2vec_dict,
    get_doc2vec,
)
from compgenlnc.generators.sequences.kmer import gen_kmer_dict, get_kmer
from compgenlnc.generators.sequences.role2vec import (
    gen_role2vec_dict,
    gen_role2vec_dict_from_files,
    gen_role2vec_dict_from_folder,
)
from compgenlnc.utils.dict_manager import load_dict
from compgenlnc.utils.fasta_manager import load_fasta
from constants import SEQ_FASTA_EXAMPLE, SEQ_SS_FOLDER


def test_gen_role2vec_dict(doc2vec_model_deterministic, tmp_path):
    sequences = tuple(load_fasta(SEQ_FASTA_EXAMPLE, 'seq'))
    doc2vec_model = doc2vec_model_deterministic

    ctd_dict = {
        record.id: get_ctd(record.seq)
        for record in sequences
    }
    doc2vec_dict = {
        record.id: get_doc2vec(record.seq, model=doc2vec_model)
        for record in sequences
    }
    kmer_dict = {
        record.id: get_kmer(record.seq, 3)
        for record in sequences
    }
    
    filename = tmp_path / 'role2vec.dict'
    gen_role2vec_dict(ctd_dict, doc2vec_dict, kmer_dict, filename)
    assert filename.exists()

    role2vec_dict = load_dict(filename)
    assert all(
        len(vector) == 256 and np.issubdtype(vector.dtype, np.float32)
        for vector in role2vec_dict.values()
    )


def test_gen_role2vec_dict_from_files(doc2vec_model_deterministic, tmp_path):
    sequences = tuple(load_fasta(SEQ_FASTA_EXAMPLE, 'seq'))
    doc2vec_model = doc2vec_model_deterministic

    ctd_file = tmp_path / 'ctd.dict'
    doc2vec_file = tmp_path / 'doc2vec.dict'
    kmer_file = tmp_path / 'kmer.dict'

    gen_ctd_dict(sequences, ctd_file)
    gen_doc2vec_dict(sequences, doc2vec_file, model=doc2vec_model)
    gen_kmer_dict(sequences, kmer_file, 3)
    
    filename = tmp_path / 'role2vec.dict'
    gen_role2vec_dict_from_files(ctd_file, doc2vec_file, kmer_file, filename)
    assert filename.exists()

    role2vec_dict = load_dict(filename)
    assert all(
        len(vector) == 256 and np.issubdtype(vector.dtype, np.float32)
        for vector in role2vec_dict.values()
    )

def test_gen_role2vec_dict_from_folder(tmp_path):
    sequences = tuple(load_fasta(SEQ_FASTA_EXAMPLE, 'seq'))
    filename = tmp_path / 'role2vec.dict'
    gen_role2vec_dict_from_folder(SEQ_SS_FOLDER, filename, 3)
    assert filename.exists()

    role2vec_dict = load_dict(filename)
    assert all(
        len(vector) == 256 and np.issubdtype(vector.dtype, np.float32)
        for vector in role2vec_dict.values()
    )
