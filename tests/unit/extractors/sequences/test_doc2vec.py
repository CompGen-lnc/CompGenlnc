import numpy as np
import pytest

import compgenlnc as cgl


def test_train_doc2vec_model(doc2vec_model_deterministic):
    model = doc2vec_model_deterministic
    assert model is not None
    vector = model.infer_vector(["ACG", "CGT", "GTA"])
    assert len(vector) == 256


@pytest.mark.parametrize(
    "seq, segments",
    [
        ("AAAAAA", ["AAA"] * 4),
        ("ACGTACGTACGT", ["ACG", "CGT", "GTA", "TAC"] * 2 + ["ACG", "CGT"]),
    ],
)
def test_get_doc2vec(seq, segments, doc2vec_model_deterministic):
    model = doc2vec_model_deterministic
    vector = cgl.get_doc2vec(seq, model=model)
    expected = model.infer_vector(segments)
    assert np.array_equal(vector, expected)


def test_get_doc2vec_by_name(example_id_list, doc2vec_model_deterministic):
    model = doc2vec_model_deterministic
    for id_ in example_id_list:
        seq = cgl.read_sequence(cgl.paths.MIRNA_SEQ_SS_FOLDER, id_)
        segments = [seq[i : i + 3] for i in range(len(seq) - 2)]
        vector = cgl.get_doc2vec_by_name(
            cgl.paths.MIRNA_SEQ_SS_FOLDER, id_, model=model
        )
        expected = model.infer_vector(segments) if segments else [0.0] * 256
        assert np.array_equal(vector, expected)


def test_gen_doc2vec_dict(seq_ss_list, doc2vec_model_deterministic, tmp_path):
    model = doc2vec_model_deterministic
    dict_file = tmp_path / "kmer.dict"
    cgl.gen_doc2vec_dict(dict_file, seq_ss_list, model=model)
    assert dict_file.exists()

    kmer_dict = cgl.load_dict(dict_file)
    assert all(
        [
            np.array_equal(
                kmer_dict[record.id], cgl.get_doc2vec(record.seq, model=model)
            )
            for record in seq_ss_list
        ]
    )


def test_gen_doc2vec_dict_from_folder(
    example_id_filter, doc2vec_model_deterministic, tmp_path
):
    model = doc2vec_model_deterministic
    dict_file = tmp_path / "doc2vec.dict"
    cgl.gen_doc2vec_dict_from_folder(
        dict_file,
        cgl.paths.MIRNA_SEQ_SS_FOLDER,
        model=model,
        keep=example_id_filter,
    )
    assert dict_file.exists()

    doc2vec_dict = cgl.load_dict(dict_file)
    assert all(
        [
            np.array_equal(
                doc2vec_dict[id_],
                cgl.get_doc2vec_by_name(
                    cgl.paths.MIRNA_SEQ_SS_FOLDER, id_, model=model
                ),
            )
            for id_ in example_id_filter
        ]
    )
