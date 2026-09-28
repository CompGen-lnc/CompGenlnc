import pytest
from constants import SEQ_FASTA_EXAMPLE

from compgenlnc.extractors.sequences import train_doc2vec_model
from compgenlnc.fileman import iter_seq_fasta


@pytest.fixture(scope="session")
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
