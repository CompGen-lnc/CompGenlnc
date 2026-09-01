import os
from pathlib import Path
from typing import Iterable

import networkx as nx
import numpy as np
from karateclub.node_embedding.structural import Role2Vec
from multiprocessing import cpu_count
from sklearn.neighbors import KDTree

from compgenlnc.generators.sequences.ctd import get_ctd
from compgenlnc.generators.sequences.doc2vec import (
    get_doc2vec, train_doc2vec_model_from_fasta, train_doc2vec_model
)
from compgenlnc.generators.sequences.kmer import get_kmer
from compgenlnc.utils.fasta_manager import load_fasta, iter_seq_ss
from compgenlnc.utils.dict_manager import load_dict, save_dict


def gen_role2vec_dict(
        ctd_dict: dict[str, np.typing.NDArray],
        doc2vec_dict: dict[str, np.typing.NDArray],
        kmer_dict: dict[str, np.typing.NDArray],
        filename: str | os.PathLike,
) -> None:
    keys = sorted(kmer_dict.keys() & doc2vec_dict.keys() & ctd_dict.keys())
    vectors = np.array([
        np.concatenate((ctd_dict[key], doc2vec_dict[key], kmer_dict[key]), dtype=np.float32)
        for key in keys
    ])
    kdt = KDTree(vectors, leaf_size=30, metric='euclidean')
    k_near = kdt.query(vectors, k=10, return_distance=False)

    g = nx.Graph()
    for key in keys:
        near = k_near[keys.index(key)][1:]
        for x in near:
            g.add_edge(x,keys.index(key))

    model = Role2Vec(dimensions=256, workers=cpu_count(),epochs=16)
    model.fit(g)
    embedding = model.get_embedding()
    embedding_dict = {}
    for i in range(len(g.nodes)):
        embedding_dict[keys[i]] = embedding[i].astype(np.float32)

    filename = Path(filename).resolve()
    save_dict(filename, embedding_dict)


def gen_role2vec_dict_from_files(
        ctd_file: str | os.PathLike,
        doc2vec_file: str | os.PathLike,
        kmer_file: str | os.PathLike,
        filename: str | os.PathLike,
) -> None:
    ctd_file = Path(ctd_file).resolve()
    doc2vec_file = Path(doc2vec_file).resolve()
    kmer_file = Path(kmer_file).resolve()
    filename = Path(filename).resolve()

    ctd_dict = load_dict(ctd_file)
    doc2vec_dict = load_dict(doc2vec_file)
    kmer_dict = load_dict(kmer_file)

    gen_role2vec_dict(ctd_dict, doc2vec_dict, kmer_dict, filename)


def gen_role2vec_dict_from_fasta(
        filename: str | os.PathLike,
        dict_file: str | os.PathLike,
        k: int,
        /, *,
        keep: Iterable[str] | None = None,
        mature_only: bool = False,
) -> None:
    filename = Path(filename).resolve()
    dict_file = Path(dict_file).resolve()

    dict_folder = dict_file.parent
    dict_folder.mkdir(parents=True, exist_ok=True)

    doc2vec_model = train_doc2vec_model_from_fasta(filename)

    ctd_dict = {
        record.id: get_ctd(record.seq)
        for record in load_fasta(filename, 'seq') if record.id in keep
    }
    doc2vec_dict = {
        record.id: get_doc2vec(record.seq, model=doc2vec_model)
        for record in load_fasta(filename, 'seq') if record.id in keep
    }
    kmer_dict = {
        record.id: get_kmer(record.seq, k)
        for record in load_fasta(filename, 'seq') if record.id in keep
    }

    gen_role2vec_dict(ctd_dict, doc2vec_dict, kmer_dict, dict_file)


def gen_role2vec_dict_from_folder(
        folder: str | os.PathLike,
        dict_file: str | os.PathLike,
        k: int,
        /, *,
        keep: Iterable[str] | None = None,
        mature_only: bool = False,
) -> None:
    folder = Path(folder).resolve()
    dict_file = Path(dict_file).resolve()

    dict_folder = dict_file.parent
    dict_folder.mkdir(parents=True, exist_ok=True)

    doc2vec_model = train_doc2vec_model(
        iter_seq_ss(folder, keep=keep, mature_only=mature_only)
    )

    ctd_dict = {
        record.id: get_ctd(record.seq)
        for record in iter_seq_ss(folder, keep=keep, mature_only=mature_only)
    }
    doc2vec_dict = {
        record.id: get_doc2vec(record.seq, model=doc2vec_model)
        for record in iter_seq_ss(folder, keep=keep, mature_only=mature_only)
    }
    kmer_dict = {
        record.id: get_kmer(record.seq, k)
        for record in iter_seq_ss(folder, keep=keep, mature_only=mature_only)
    }

    gen_role2vec_dict(ctd_dict, doc2vec_dict, kmer_dict, dict_file)
