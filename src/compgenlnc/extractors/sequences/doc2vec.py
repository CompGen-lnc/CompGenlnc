import os
from collections.abc import Iterable
from pathlib import Path

import numpy as np
from gensim.models.doc2vec import Doc2Vec, TaggedDocument

from compgenlnc.config.paths import SEQ_SS_PREFIX
from compgenlnc.config.seeds import DOC2VEC_MODEL_SEED
from compgenlnc.fileman.dict_manager import save_dict
from compgenlnc.fileman.fasta_manager import (
    iter_seq_ss,
    load_fasta,
    read_sequence,
)
from compgenlnc.structs.seq_record import SeqRecord
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.molecules import SeqLike


def train_doc2vec_model(
    sequences: Iterable[SeqSSRecord],
    model_file: str | os.PathLike = "",
) -> Doc2Vec:
    tokens = [
        TaggedDocument(
            [record.seq[j : j + 3] for j in range(len(record.seq) - 2)], [i]
        )
        for i, record in enumerate(sequences)
    ]
    model = Doc2Vec(vector_size=256, min_count=3, epochs=100, workers=12)
    model.build_vocab(tokens)
    model.train(tokens, total_examples=model.corpus_count, epochs=model.epochs)
    if model_file:
        model_file = Path(model_file).resolve()
        model_file.parent.mkdir(parents=True, exist_ok=True)
        model.save(str(model_file))
    return model


def train_doc2vec_model_from_fasta(
    filename: str | os.PathLike,
    model_file: str | os.PathLike = "",
) -> Doc2Vec:
    filename = Path(filename).resolve()
    model_file = Path(model_file).resolve()
    return train_doc2vec_model(load_fasta(filename, mode="seq"), model_file)


def get_doc2vec(
    seq: SeqLike,
    *,
    model: Doc2Vec | None = None,
    model_file: str | os.PathLike = "",
) -> np.typing.NDArray[np.float32]:
    if isinstance(seq, SeqSSRecord):
        seq = str(seq.seq)
    if isinstance(seq, SeqRecord):
        seq = str(seq)
    if not seq:
        return np.array([0.0] * 256, np.float32)
    doc = [seq[i : i + 3] for i in range(len(seq) - 2)]

    if model is None:
        model = Doc2Vec.load(model_file)

    model.random.seed(DOC2VEC_MODEL_SEED)
    return model.infer_vector(doc)


def get_doc2vec_by_name(
    id_: str,
    folder: str | os.PathLike,
    /,
    *,
    mature_only: bool = False,
    model: Doc2Vec | None = None,
    model_file: str | os.PathLike = "",
) -> np.typing.NDArray[np.float32]:
    folder = Path(folder).resolve()
    seq = read_sequence(folder, id_, mature_only=mature_only)
    return get_doc2vec(seq, model=model, model_file=model_file)


def gen_doc2vec_dict(
    record_list: Iterable[SeqSSRecord],
    dict_file: str | os.PathLike,
    /,
    *,
    model: Doc2Vec | None = None,
    model_file: str | os.PathLike = "",
    keep: Iterable[str] | None = None,
    mature_only: bool = False,
) -> None:
    dict_file = Path(dict_file).resolve()

    dict_folder = dict_file.parent
    dict_folder.mkdir(parents=True, exist_ok=True)

    doc2vec_dict = {
        record.id: get_doc2vec(record.seq, model=model, model_file=model_file)
        for record in record_list
    }
    save_dict(dict_file, doc2vec_dict)


def gen_doc2vec_dict_from_fasta(
    filename: str | os.PathLike,
    dict_file: str | os.PathLike,
    /,
    *,
    model: Doc2Vec | None = None,
    model_file: str | os.PathLike = "",
    keep: Iterable[str] | None = None,
    mature_only: bool = False,
) -> None:
    dict_file = Path(dict_file).resolve()
    dict_folder = dict_file.parent
    dict_folder.mkdir(parents=True, exist_ok=True)

    if model is None:
        model = Doc2Vec.load(str(model_file))

    doc2vec_dict = {
        record.id: get_doc2vec(record.seq, model=model)
        for record in load_fasta(filename, mode="seq", keep=keep)
    }
    save_dict(dict_file, doc2vec_dict)


def gen_doc2vec_dict_from_folder(
    folder: str | os.PathLike,
    dict_file: str | os.PathLike,
    /,
    *,
    model: Doc2Vec | None = None,
    model_file: str | os.PathLike = "",
    keep: Iterable[str] | None = None,
    mature_only: bool = False,
) -> None:
    folder = Path(folder).resolve()
    dict_file = Path(dict_file).resolve()

    dict_folder = dict_file.parent
    dict_folder.mkdir(parents=True, exist_ok=True)

    doc2vec_dict = {
        record.id: get_doc2vec(record.seq, model=model, model_file=model_file)
        for record in iter_seq_ss(folder, keep=keep, mature_only=mature_only)
    }
    save_dict(dict_file, doc2vec_dict)
