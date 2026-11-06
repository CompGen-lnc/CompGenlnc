from __future__ import annotations

from collections.abc import Collection, Iterable
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from gensim.models.doc2vec import Doc2Vec, TaggedDocument
from numpy.typing import NDArray

from compgenlnc.config import DOC2VEC_MODEL_SEED
from compgenlnc.consts.params import SEQ_MODE
from compgenlnc.fileman.dict_manager import save_dict
from compgenlnc.fileman.fasta_manager import (
    iter_seq_ss,
    load_fasta,
    read_sequence,
)
from compgenlnc.fileman.log_manager import send_msg, track_progress
from compgenlnc.structs.seq_record import SeqRecord
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.molecules import SeqLike

if TYPE_CHECKING:
    from _typeshed import StrPath


def train_doc2vec_model(
    sequences: Iterable[SeqSSRecord],
    /,
    model_file: StrPath | None = None,
) -> Doc2Vec:
    tokens = [
        TaggedDocument(
            [record.seq[j : j + 3] for j in range(len(record.seq) - 2)], [i]
        )
        for i, record in enumerate(sequences)
    ]
    model = Doc2Vec(vector_size=256, min_count=10, epochs=16, workers=12)
    model.build_vocab(tokens)
    model.train(tokens, total_examples=model.corpus_count, epochs=model.epochs)
    if model_file is not None:
        model_file = Path(model_file).resolve()
        model_file.parent.mkdir(parents=True, exist_ok=True)
        model.save(str(model_file))
    return model


def train_doc2vec_model_from_fasta(
    fasta_file: StrPath,
    /,
    model_file: StrPath | None = None,
) -> Doc2Vec:
    fasta_file = Path(fasta_file).resolve()
    return train_doc2vec_model(
        load_fasta(fasta_file, mode=SEQ_MODE), model_file
    )


def get_doc2vec(
    seq: SeqLike,
    /,
    *,
    model: Doc2Vec | None = None,
    model_file: StrPath | None = None,
) -> NDArray[np.float32]:
    if isinstance(seq, SeqSSRecord):
        seq = str(seq.seq)
    if isinstance(seq, SeqRecord):
        seq = str(seq)
    if not seq:
        return np.array([0.0] * 256, np.float32)
    doc = [seq[i : i + 3] for i in range(len(seq) - 2)]

    if model is None:
        model = Doc2Vec.load(str(model_file))

    model.random.seed(DOC2VEC_MODEL_SEED)
    return model.infer_vector(doc)


def get_doc2vec_by_name(
    folder: StrPath,
    id_: str,
    *,
    model: Doc2Vec | None = None,
    model_file: StrPath | None = None,
) -> NDArray[np.float32]:
    folder = Path(folder).resolve()
    seq = read_sequence(folder, id_)
    return get_doc2vec(seq, model=model, model_file=model_file)


def gen_doc2vec_dict(
    dict_file: StrPath,
    records: Collection[SeqSSRecord],
    *,
    model: Doc2Vec | None = None,
    model_file: StrPath | None = None,
) -> None:
    send_msg("Generating Doc2Vec dictionary\n")

    if model is None:
        if model_file is not None and Path(model_file).resolve().exists():
            send_msg("  Loading trained model")
            model = Doc2Vec.load(str(model_file))
            complete_msg = "\r  Model loaded in memory\n"
        else:
            send_msg("  Training data model")
            model = train_doc2vec_model(records, model_file)
            complete_msg = "\r  Data model trained \n"
        send_msg(complete_msg)

    desc = "  Calculating Doc2Vec"
    records = track_progress(records, desc, unit="seq")

    doc2vec_dict = {
        record.id: get_doc2vec(record.seq, model=model, model_file=model_file)
        for record in records
    }
    save_dict(dict_file, doc2vec_dict, overwrite=True)

    send_msg("  Doc2Vec dictionary generated\n")


def gen_doc2vec_dict_from_fasta(
    dict_file: StrPath,
    fasta_file: StrPath,
    *,
    model: Doc2Vec | None = None,
    model_file: StrPath | None = None,
    keep: Collection[str] | None = None,
    exclude: Collection[str] | None = None,
) -> None:
    records = load_fasta(fasta_file, mode=SEQ_MODE, keep=keep, exclude=exclude)
    gen_doc2vec_dict(dict_file, records, model=model, model_file=model_file)


def gen_doc2vec_dict_from_folder(
    dict_file: StrPath,
    seq_ss_folder: StrPath,
    *,
    model: Doc2Vec | None = None,
    model_file: StrPath | None = None,
    keep: Collection[str] | None = None,
    exclude: Collection[str] | None = None,
) -> None:
    records = iter_seq_ss(seq_ss_folder, keep=keep, exclude=exclude)
    gen_doc2vec_dict(dict_file, records, model=model, model_file=model_file)
