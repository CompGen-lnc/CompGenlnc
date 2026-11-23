import pytest

import compgenlnc as cgl


@pytest.fixture(scope="module")
def doc2vec_model():
    return cgl.train_doc2vec_model(
        cgl.iter_fasta(cgl.paths.PRE_MIRNA_FASTA, mode = cgl.SEQ_MODE)
    )


@pytest.fixture
def doc2vec_model_deterministic(doc2vec_model):
    model = doc2vec_model
    infer_original = model.infer_vector

    def fixed_infer(doc_words, **kwargs):
        model.random.seed(0)
        return infer_original(doc_words, **kwargs)

    model.infer_vector = fixed_infer

    yield model
