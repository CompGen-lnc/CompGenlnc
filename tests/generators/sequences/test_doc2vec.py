import numpy as np
import pytest

from compgenlnc.generators.sequences.doc2vec import (
    get_doc2vec,
    get_doc2vec_by_name,
    gen_doc2vec_dict,
    gen_doc2vec_dict_from_folder,
)
from compgenlnc.utils.dict_manager import load_dict
from compgenlnc.utils.fasta_manager import iter_seq_fasta, iter_seq_ss
from constants import SEQ_FASTA_EXAMPLE, SEQ_SS_FOLDER


def test_train_doc2vec_model(doc2vec_model_deterministic):
    model = doc2vec_model_deterministic
    assert model is not None
    vector = model.infer_vector(['ACG', 'CGT', 'GTA'])
    assert len(vector) == model.vector_size

@pytest.mark.parametrize(
    'seq, segments',
    [
        ('AAAAAA', ['AAA'] * 4),
        ('ACGTACGTACGT', ['ACG', 'CGT', 'GTA', 'TAC'] * 2 + ['ACG', 'CGT']),
    ]
)
def test_get_doc2vec(seq, segments, doc2vec_model_deterministic):
    model = doc2vec_model_deterministic
    vector = get_doc2vec(seq, model=model)
    expected = model.infer_vector(segments)
    assert np.array_equal(vector, expected)


@pytest.mark.parametrize(
    'id_, segments',
    [
        ('hsa1', []),
        ('hsa2', ['ACG', 'CGT', 'GTG', 'TGC', 'GCA']),
        ('mol1', ['ACG', 'CGT', 'GTA', 'TAC'] * 2 + ['ACG', 'CGT']),
        ('mol2', ['AAA'] * 2),
        ('mol3', ['GCG', 'CGC'] * 3),
    ]
)
def test_get_doc2vec_by_name(id_, segments, doc2vec_model_deterministic):
    model = doc2vec_model_deterministic
    vector = get_doc2vec_by_name(id_, SEQ_SS_FOLDER, model=model)
    expected = model.infer_vector(segments) if segments else [0.0] * 256
    assert np.array_equal(vector, expected)


@pytest.mark.parametrize(
    'record_list',
    [
        (iter_seq_ss(SEQ_SS_FOLDER, keep=['hsa1', 'hsa2', 'mol1'])),
        (iter_seq_ss(SEQ_SS_FOLDER, keep=['hsa1', 'mol3', 'mol2'])),
        (iter_seq_ss(SEQ_SS_FOLDER, keep=['hsa2', 'mol2'])),
        (iter_seq_ss(SEQ_SS_FOLDER, keep=['hsa1', 'hsa2', 'mol1'])),
    ]
)
def test_gen_doc2vec_dict(record_list, doc2vec_model_deterministic, tmp_path):
    model = doc2vec_model_deterministic
    dict_file = tmp_path / 'kmer.dict'
    gen_doc2vec_dict(record_list, dict_file, model=model)
    assert dict_file.exists()

    kmer_dict = load_dict(dict_file)
    assert all([
        np.array_equal(
            kmer_dict[record.id],
            get_doc2vec(record.seq, model=model)
        )
        for record in record_list
    ])


@pytest.mark.parametrize(
    'seq_list',
    [
        (['hsa1', 'hsa2', 'mol1']),
        (['hsa1', 'mol3', 'mol2']),
        (['hsa2', 'mol2']),
        (['hsa1', 'hsa2', 'mol1']),
    ]
)
def test_gen_doc2vec_dict_from_folder(
        seq_list, doc2vec_model_deterministic, tmp_path
):
    model = doc2vec_model_deterministic
    dict_file = tmp_path / 'doc2vec.dict'
    gen_doc2vec_dict_from_folder(
        SEQ_SS_FOLDER, dict_file, model=model, keep=seq_list
    )
    assert dict_file.exists()

    doc2vec_dict = load_dict(dict_file)
    assert all([
        np.array_equal(
            doc2vec_dict[id_],
            get_doc2vec_by_name(id_, SEQ_SS_FOLDER, model=model)
        )
        for id_ in seq_list
    ])
