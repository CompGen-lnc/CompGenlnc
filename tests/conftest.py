import pytest

from compgenlnc.config.seeds import DOC2VEC_MODEL_SEED
from compgenlnc.generators.sequences.doc2vec import train_doc2vec_model
from compgenlnc.utils.fasta_manager import (
    load_fasta,
    iter_seq_fasta,
    iter_seq_ss
)

from constants import SEQ_FASTA_EXAMPLE, SEQ_SS_FOLDER
from params import id_list


@pytest.fixture(scope='session')
def doc2vec_model():
    return train_doc2vec_model(iter_seq_fasta(SEQ_FASTA_EXAMPLE))


@pytest.fixture
def doc2vec_model_deterministic(doc2vec_model):
    model = doc2vec_model
    infer_original = model.infer_vector

    def fixed_infer(doc_words, **kwargs):
        model.random.seed(0)
        return infer_original(doc_words, **kwargs)

    model.infer_vector = fixed_infer

    yield model


@pytest.fixture(scope='session')
def seq_fasta_loader():
    return list(load_fasta(SEQ_FASTA_EXAMPLE, mode='seq'))


@pytest.fixture(scope='session')
def seq_fasta_list():
    return list(iter_seq_fasta(SEQ_FASTA_EXAMPLE))


@pytest.fixture
def seq_ss_list(scope='session'):
    return list(iter_seq_ss(SEQ_SS_FOLDER))


@pytest.fixture(
    params=id_list
)
def mirna_id(request):
    return request.param


@pytest.fixture(
    params=[
        iter_seq_ss(SEQ_SS_FOLDER),
        iter_seq_ss(SEQ_SS_FOLDER, keep=[id_list[i] for i in [0, 1, 2]]),
        iter_seq_ss(SEQ_SS_FOLDER, keep=[id_list[i] for i in [2, 4, 6, 8]]),
        iter_seq_ss(SEQ_SS_FOLDER, keep=[id_list[i] for i in [7, 8]]),
    ]
)
def seq_ss_filtered_list(request):
    return request.param


@pytest.fixture(
    params=[
        [],
        [id_list[i] for i in [0, 1, 2]],
        [id_list[i] for i in [2, 4, 6, 8]],
        [id_list[i] for i in [7, 8]],
    ]
)
def mirna_id_filter(request):
    return request.param
