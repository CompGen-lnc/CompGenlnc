from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

from compgenlnc.consts.params import LNCRNA, MATURE, PRECURSOR
from compgenlnc.consts.paths import (
    CTD_DICT,
    DATASET_DIR,
    DOC2VEC_DICT,
    KMER_DICT,
    LNCRNA_DICTS_DIR,
    LNCRNA_DOC2VEC_MODEL,
    MAT_MIRNA_DICTS_DIR,
    MAT_MIRNA_DOC2VEC_MODEL,
    NORMAL_SEQ_DICT,
    PRE_MIRNA_DICTS_DIR,
    PRE_MIRNA_DOC2VEC_MODEL,
    ROLE2VEC_DICT,
)
from compgenlnc.dict_indexer import dict_to_numpy
from compgenlnc.extractors.sequences.ctd import gen_ctd_dict
from compgenlnc.extractors.sequences.doc2vec import gen_doc2vec_dict
from compgenlnc.extractors.sequences.kmer import gen_kmer_dict
from compgenlnc.extractors.sequences.normal_sequence import (
    gen_normalized_sequence_dict,
)
from compgenlnc.extractors.sequences.role2vec import (
    gen_role2vec_dict_from_files,
)
from compgenlnc.fileman.dict_manager import load_dict
from compgenlnc.structs.seq_ss_record import SeqSSRecord

if TYPE_CHECKING:
    from _typeshed import StrPath


def gen_seq_features(
    records: Sequence[SeqSSRecord],
    dicts_dir: StrPath,
    k: int,
    dv2_model_file: StrPath,
    *,
    rna: str,
) -> None:
    dicts_dir = Path(dicts_dir).resolve()
    dicts_dir.parent.mkdir(parents=True, exist_ok=True)

    kmer_path = dicts_dir / KMER_DICT
    ctd_path = dicts_dir / CTD_DICT
    d2v_path = dicts_dir / DOC2VEC_DICT
    r2v_path = dicts_dir / ROLE2VEC_DICT
    normal_path = dicts_dir / NORMAL_SEQ_DICT

    gen_kmer_dict(kmer_path, records, k)
    gen_ctd_dict(ctd_path, records)
    gen_doc2vec_dict(d2v_path, records, model_file=dv2_model_file)
    gen_role2vec_dict_from_files(r2v_path, kmer_path, ctd_path, d2v_path)
    gen_normalized_sequence_dict(normal_path, records)

    def path_to_npy(filename):
        return (DATASET_DIR / f"{rna}_{filename}").with_suffix(".npy")

    kmer_npy = path_to_npy(KMER_DICT)
    ctd_npy = path_to_npy(CTD_DICT)
    d2v_npy = path_to_npy(DOC2VEC_DICT)
    r2v_npy = path_to_npy(ROLE2VEC_DICT)
    normal_npy = path_to_npy(NORMAL_SEQ_DICT)
    dict_to_numpy(kmer_npy, load_dict(kmer_path), rna=rna)
    dict_to_numpy(ctd_npy, load_dict(ctd_path), rna=rna)
    dict_to_numpy(d2v_npy, load_dict(d2v_path), rna=rna)
    dict_to_numpy(r2v_npy, load_dict(r2v_path), rna=rna)
    dict_to_numpy(normal_npy, load_dict(normal_path), rna=rna)


def gen_lncrna_seq_features(records: Sequence[SeqSSRecord]) -> None:
    gen_seq_features(
        records,
        LNCRNA_DICTS_DIR,
        k=4,
        dv2_model_file=LNCRNA_DOC2VEC_MODEL,
        rna=LNCRNA,
    )


def gen_pre_mirna_seq_features(records: Sequence[SeqSSRecord]) -> None:
    gen_seq_features(
        records,
        PRE_MIRNA_DICTS_DIR,
        k=3,
        dv2_model_file=PRE_MIRNA_DOC2VEC_MODEL,
        rna=PRECURSOR,
    )


def gen_mat_mirna_seq_features(records: Sequence[SeqSSRecord]) -> None:
    gen_seq_features(
        records,
        MAT_MIRNA_DICTS_DIR,
        k=3,
        dv2_model_file=MAT_MIRNA_DOC2VEC_MODEL,
        rna=MATURE,
    )


def gen_seq_database(
    lnc_records: Sequence[SeqSSRecord],
    pre_mir_records: Sequence[SeqSSRecord],
    mat_mir_records: Sequence[SeqSSRecord],
) -> None:
    gen_lncrna_seq_features(lnc_records)
    gen_pre_mirna_seq_features(pre_mir_records)
    gen_mat_mirna_seq_features(mat_mir_records)
