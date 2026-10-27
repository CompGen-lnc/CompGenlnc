import os
from collections.abc import Sequence
from pathlib import Path

from compgenlnc.consts.params import LNCRNA, MATURE, PRECURSOR
from compgenlnc.consts.paths import (
    CTD_DICT,
    DATASET_FOLDER,
    DOC2VEC_DICT,
    KMER_DICT,
    LNCRNA_DICTS_FOLDER,
    LNCRNA_DOC2VEC_MODEL,
    MAT_MIRNA_DICTS_FOLDER,
    MAT_MIRNA_DOC2VEC_MODEL,
    NORMAL_SEQ_DICT,
    PRE_MIRNA_DICTS_FOLDER,
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


def gen_seq_features(
    record_list: Sequence[SeqSSRecord],
    dicts_folder: str | os.PathLike,
    k: int,
    dv2_model_file: str | os.PathLike,
    *,
    rna: str,
) -> None:
    dicts_folder = Path(dicts_folder).resolve()
    dicts_folder.parent.mkdir(parents=True, exist_ok=True)

    kmer_path = dicts_folder / KMER_DICT
    ctd_path = dicts_folder / CTD_DICT
    d2v_path = dicts_folder / DOC2VEC_DICT
    r2v_path = dicts_folder / ROLE2VEC_DICT
    normal_path = dicts_folder / NORMAL_SEQ_DICT

    gen_kmer_dict(record_list, kmer_path, k)
    gen_ctd_dict(record_list, ctd_path)
    gen_doc2vec_dict(record_list, d2v_path, model_file=dv2_model_file)
    gen_role2vec_dict_from_files(kmer_path, ctd_path, d2v_path, r2v_path)
    gen_normalized_sequence_dict(record_list, normal_path)

    def path_to_npy(filename):
        return (DATASET_FOLDER / f"{rna}_{filename}").with_suffix(".npy")

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


def gen_lncrna_seq_features(record_list: Sequence[SeqSSRecord]) -> None:
    gen_seq_features(
        record_list,
        LNCRNA_DICTS_FOLDER,
        k=4,
        dv2_model_file=LNCRNA_DOC2VEC_MODEL,
        rna=LNCRNA,
    )


def gen_pre_mirna_seq_features(record_list: Sequence[SeqSSRecord]) -> None:
    gen_seq_features(
        record_list,
        PRE_MIRNA_DICTS_FOLDER,
        k=3,
        dv2_model_file=PRE_MIRNA_DOC2VEC_MODEL,
        rna=PRECURSOR,
    )


def gen_mat_mirna_seq_features(record_list: Sequence[SeqSSRecord]) -> None:
    gen_seq_features(
        record_list,
        MAT_MIRNA_DICTS_FOLDER,
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
