from __future__ import annotations

from collections.abc import Collection, Iterable
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from gensim.models.doc2vec import Doc2Vec, TaggedDocument
from numpy.typing import NDArray

import compgenlnc.config as config
from compgenlnc.consts.params import SEQ_MODE
from compgenlnc.fileman.dict_manager import save_dict
from compgenlnc.fileman.fasta_manager import (
    iter_fasta,
    load_fasta,
    load_seq_ss,
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
    *,
    model_file: StrPath | None = None,
) -> Doc2Vec:
    """Train a Doc2Vec model from an RNA sequence set.

    Parameters
    ----------
    sequences : iterable object of SeqSSRecord
        Group of RNA sequences to train the model with.
    model_file : str or os.PathLike, optional
        Path of the file to save the model. Default: `None`.

    Returns
    -------
    gensim.models.Doc2Vec
        The model to extract the Doc2Vec vectors from the sequences.
    """
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
    *,
    model_file: StrPath | None = None,
) -> Doc2Vec:
    """Train a Doc2Vec model from an RNA sequence set.

    Parameters
    ----------
    fasta_file : str or os.PathLike
        Path of the FASTA file containing the sequences.
    model_file : str or os.PathLike, optional
        Path of the file to save the model. Default: `None`.

    Returns
    -------
    gensim.models.Doc2Vec
        The model to extract the Doc2Vec from the sequences.
    """
    fasta_file = Path(fasta_file).resolve()
    return train_doc2vec_model(
        iter_fasta(fasta_file, mode=SEQ_MODE), model_file
    )


def get_doc2vec(
    seq: SeqLike,
    /,
    *,
    model: Doc2Vec | None = None,
    model_file: StrPath | None = None,
) -> NDArray[np.float32]:
    """Get the Doc2Vec of a RNA sequence.

    Parameters
    ----------
    seq : SeqLike
        The RNA sequence string on AUCG alphabet.
    model : Doc2Vec, optional
        The model to extract the Doc2Vec from `seq`. Cannot be `None` if
        model_file is `None`. Default: `None`.
    model_file : str or os.PathLike, optional
        Path of the file where the model is stored. Default: `None`.

    Returns
    -------
    numpy.ndarray of float32
        A 1D array with the Doc2Vec of the sequence.
    """
    if isinstance(seq, SeqSSRecord):
        seq = str(seq.seq)
    if isinstance(seq, SeqRecord):
        seq = str(seq)
    if not seq:
        return np.array([0.0] * 256, np.float32)
    doc = [seq[i : i + 3] for i in range(len(seq) - 2)]

    if model is None:
        model = Doc2Vec.load(str(model_file))

    model.random.seed(config.doc2vec_seed)
    return model.infer_vector(doc)


def get_doc2vec_by_name(
    dir: StrPath,
    id_: str,
    *,
    model: Doc2Vec | None = None,
    model_file: StrPath | None = None,
) -> NDArray[np.float32]:
    """Get the Doc2Vec of a RNA sequence.

    Parameters
    ----------
    dir : str or os.PathLike
        Path of the directory where the seq+ss file is.
    id_ : str
        ID of the RNA sequence to read.
    model : Doc2Vec, optional
        The model to extract the Doc2Vec from the sequence. Cannot be
        `None` if model_file is `None`. Default: `None`.
    model_file : str or os.PathLike, optional
        Path of the file where the model is stored. Default: `None`.

    Returns
    -------
    numpy.ndarray of float32
        A 1D array with the Doc2Vec of the sequence.
    """
    dir = Path(dir).resolve()
    seq = read_sequence(dir, id_)
    return get_doc2vec(seq, model=model, model_file=model_file)


def gen_doc2vec_dict(
    dict_file: StrPath,
    records: Collection[SeqSSRecord],
    *,
    model: Doc2Vec | None = None,
    model_file: StrPath | None = None,
) -> None:
    """Generate a dictionary of Doc2Vec from an RNA sequence set.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    records : collection of SeqSSRecord
        Set containing all the sequence records to save into the
        dictionary.
    model : Doc2Vec, optional
        The model to extract the Doc2Vec from the sequences. Cannot be
        `None` if model_file is `None`. Default: `None`.
    model_file : str or os.PathLike, optional
        Path of the file where the model is stored. Default: `None`.
    """
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
    """Generate a dictionary of Doc2Vec from a FASTA file.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    fasta_file : str or os.PathLike
        Path of the FASTA file with the sequences.
    model : Doc2Vec, optional
        The model to extract the Doc2Vec from the sequences. Cannot be
        `None` if model_file is `None`. Default: `None`.
    model_file : str or os.PathLike, optional
        Path of the file where the model is stored. Default: `None`.
    keep : collection of str, optional
        Set of IDs of the RNAs to save. Do not apply if `None` or empty.
        Default: `None`.
    exclude : collection of str, optional
        Set of IDs of the RNAs to ignore. Do not apply if `None` or
        empty. Default: `None`.
    """
    records = load_fasta(fasta_file, mode=SEQ_MODE, keep=keep, exclude=exclude)
    gen_doc2vec_dict(dict_file, records, model=model, model_file=model_file)


def gen_doc2vec_dict_from_dir(
    dict_file: StrPath,
    seq_ss_dir: StrPath,
    *,
    model: Doc2Vec | None = None,
    model_file: StrPath | None = None,
    keep: Collection[str] | None = None,
    exclude: Collection[str] | None = None,
) -> None:
    """Generate a dictionary of Doc2Vec from seq+ss files.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    seq_ss_dir : str or os.PathLike
        Path of the directory where the seq+ss files are.
    model : Doc2Vec, optional
        The model to extract the Doc2Vec from the sequences. Cannot be
        `None` if model_file is `None`. Default: `None`.
    model_file : str or os.PathLike, optional
        Path of the file where the model is stored. Default: `None`.
    keep : collection of str, optional
        Set of IDs of the RNAs to save. Do not apply if `None` or empty.
        Default: `None`.
    exclude : collection of str, optional
        Set of IDs of the RNAs to ignore. Do not apply if `None` or
        empty. Default: `None`.
    """
    records = load_seq_ss(seq_ss_dir, keep=keep, exclude=exclude)
    gen_doc2vec_dict(dict_file, records, model=model, model_file=model_file)
