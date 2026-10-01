import os
from collections.abc import Iterable
from multiprocessing import cpu_count
from pathlib import Path

import numpy as np
import torch
from torch_geometric.data import Data
from torch_geometric.nn import knn_graph, Node2Vec

from compgenlnc.consts.params import WORKERS
from compgenlnc.extractors.sequences.ctd import get_ctd
from compgenlnc.extractors.sequences.doc2vec import (
    get_doc2vec,
    train_doc2vec_model,
    train_doc2vec_model_from_fasta,
)
from compgenlnc.extractors.sequences.kmer import get_kmer
from compgenlnc.fileman.dict_manager import load_dict, save_dict
from compgenlnc.fileman.fasta_manager import iter_seq_ss, load_fasta


def gen_role2vec_dict(
    kmer_dict: dict[str, np.typing.NDArray],
    ctd_dict: dict[str, np.typing.NDArray],
    doc2vec_dict: dict[str, np.typing.NDArray],
    filename: str | os.PathLike,
    /,
) -> None:
    keys = sorted(kmer_dict.keys() & doc2vec_dict.keys() & ctd_dict.keys())
    vectors = [
        np.concatenate(
            (ctd_dict[key], doc2vec_dict[key], kmer_dict[key]),
            dtype=np.float32,
        )
        for key in keys
    ]
    tensor = torch.tensor(np.stack(vectors, dtype=np.float32))

    g = knn_graph(tensor, k=10, num_workers=WORKERS)
    data = Data(edge_index=g, num_nodes=tensor.size(0))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = Node2Vec(
        edge_index=data.edge_index,
        embedding_dim=256,
        walk_length=20,
        context_size=10,
        sparse=True,
    ).to(device)
    loader = model.loader(batch_size=128, shuffle=True, num_workers=WORKERS)
    optimizer = torch.optim.SparseAdam(model.parameters())

    for epoch in range(16):
        model.train()
        total_loss = 0
        for pos_rw, neg_rw in loader:
            optimizer.zero_grad()
            loss = model.loss(pos_rw.to(device), neg_rw.to(device))
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        # return total_loss / len(loader)

    embedding = model.embedding.weight.data.numpy(force=True)
    save_dict(filename, keys=keys, values=embedding)


def gen_role2vec_dict_from_files(
    kmer_file: str | os.PathLike,
    ctd_file: str | os.PathLike,
    doc2vec_file: str | os.PathLike,
    filename: str | os.PathLike,
    /,
) -> None:
    kmer_file = Path(kmer_file).resolve()
    ctd_file = Path(ctd_file).resolve()
    doc2vec_file = Path(doc2vec_file).resolve()
    filename = Path(filename).resolve()

    kmer_dict = load_dict(kmer_file)
    ctd_dict = load_dict(ctd_file)
    doc2vec_dict = load_dict(doc2vec_file)

    gen_role2vec_dict(kmer_dict, ctd_dict, doc2vec_dict, filename)


def gen_role2vec_dict_from_fasta(
    filename: str | os.PathLike,
    dict_file: str | os.PathLike,
    k: int,
    /,
    *,
    keep: Iterable[str] | None = None,
    mature_only: bool = False,
) -> None:
    filename = Path(filename).resolve()
    dict_file = Path(dict_file).resolve()

    dict_folder = dict_file.parent
    dict_folder.mkdir(parents=True, exist_ok=True)

    doc2vec_model = train_doc2vec_model_from_fasta(filename)

    kmer_dict = {
        record.id: get_kmer(record.seq, k)
        for record in load_fasta(filename, "seq")
        if record.id in keep
    }
    ctd_dict = {
        record.id: get_ctd(record.seq)
        for record in load_fasta(filename, "seq")
        if record.id in keep
    }
    doc2vec_dict = {
        record.id: get_doc2vec(record.seq, model=doc2vec_model)
        for record in load_fasta(filename, "seq")
        if record.id in keep
    }

    gen_role2vec_dict(kmer_dict, ctd_dict, doc2vec_dict, dict_file)


def gen_role2vec_dict_from_folder(
    folder: str | os.PathLike,
    dict_file: str | os.PathLike,
    k: int,
    /,
    *,
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

    kmer_dict = {
        record.id: get_kmer(record.seq, k)
        for record in iter_seq_ss(folder, keep=keep, mature_only=mature_only)
    }
    ctd_dict = {
        record.id: get_ctd(record.seq)
        for record in iter_seq_ss(folder, keep=keep, mature_only=mature_only)
    }
    doc2vec_dict = {
        record.id: get_doc2vec(record.seq, model=doc2vec_model)
        for record in iter_seq_ss(folder, keep=keep, mature_only=mature_only)
    }

    gen_role2vec_dict(kmer_dict, ctd_dict, doc2vec_dict, dict_file)
