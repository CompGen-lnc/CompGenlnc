import numpy as np
import pytest
from constants import SEQ_SS_FOLDER

from compgenlnc.config.paths import SEQ_SS_PREFIX
from compgenlnc.generators.sequences.doc2vec import (
    gen_doc2vec_dict,
    gen_doc2vec_dict_from_folder,
    get_doc2vec,
    get_doc2vec_by_name,
)
from compgenlnc.utils.dict_manager import load_dict
from compgenlnc.utils.fasta_manager import read_sequence


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
    vector = get_doc2vec(seq, model=model)
    expected = model.infer_vector(segments)
    assert np.array_equal(vector, expected)


def test_get_doc2vec_by_name(example_id, doc2vec_model_deterministic):
    model = doc2vec_model_deterministic
    seq = read_sequence(SEQ_SS_FOLDER / f"{SEQ_SS_PREFIX}{example_id}.dat")
    segments = [seq[i : i + 3] for i in range(len(seq) - 2)]
    vector = get_doc2vec_by_name(example_id, SEQ_SS_FOLDER, model=model)
    expected = model.infer_vector(segments) if segments else [0.0] * 256
    assert np.array_equal(vector, expected)


def test_gen_doc2vec_dict(
    seq_ss_filtered_list, doc2vec_model_deterministic, tmp_path
):
    model = doc2vec_model_deterministic
    dict_file = tmp_path / "kmer.dict"
    gen_doc2vec_dict(seq_ss_filtered_list, dict_file, model=model)
    assert dict_file.exists()

    kmer_dict = load_dict(dict_file)
    assert all(
        [
            np.array_equal(
                kmer_dict[record.id], get_doc2vec(record.seq, model=model)
            )
            for record in seq_ss_filtered_list
        ]
    )


def test_gen_doc2vec_dict_from_folder(
    example_id_filter, doc2vec_model_deterministic, tmp_path
):
    model = doc2vec_model_deterministic
    dict_file = tmp_path / "doc2vec.dict"
    gen_doc2vec_dict_from_folder(
        SEQ_SS_FOLDER, dict_file, model=model, keep=example_id_filter
    )
    assert dict_file.exists()

    doc2vec_dict = load_dict(dict_file)
    assert all(
        [
            np.array_equal(
                doc2vec_dict[id_],
                get_doc2vec_by_name(id_, SEQ_SS_FOLDER, model=model),
            )
            for id_ in example_id_filter
        ]
    )
