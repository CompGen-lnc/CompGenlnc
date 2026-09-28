import os
import sys
from pathlib import Path

from compgenlnc.config.paths import (
    CTD_DICT,
    DOC2VEC_DICT,
    KMER_DICT,
    NORMAL_SEQ_DICT,
    ROLE2VEC_DICT,
)
from compgenlnc.extractors.sequences import (
    gen_ctd_dict_from_fasta,
    gen_doc2vec_dict_from_fasta,
    gen_kmer_dict_from_fasta,
    gen_normalized_sequence_dict_from_fasta,
    gen_role2vec_dict,
    train_doc2vec_model_from_fasta,
)
from compgenlnc.utils import load_dict

raiz_proyecto = str(Path(__file__).resolve().parents[2])
if raiz_proyecto not in sys.path:
    sys.path.append(raiz_proyecto)

from tests.constants import (
    EMBEDDING_FOLDER,
    LNCRNA_FASTA,
    LNCRNA_FOLDER,
    MAT_MIRNA_FASTA,
    MIRNA_FOLDER,
    PRE_MIRNA_FASTA,
)


def generate_lnc_seq_characterization():
    folder = LNCRNA_FOLDER / "dict"
    kmer_dict = folder / KMER_DICT
    ctd_dict = folder / CTD_DICT
    doc2vec_dict = folder / DOC2VEC_DICT
    normal_dict = folder / NORMAL_SEQ_DICT
    role2vec_dict = folder / ROLE2VEC_DICT

    if not kmer_dict.exists():
        gen_kmer_dict_from_fasta(LNCRNA_FASTA, kmer_dict, 4)
    if not ctd_dict.exists():
        gen_ctd_dict_from_fasta(LNCRNA_FASTA, ctd_dict)
    if not doc2vec_dict.exists():
        model_path = EMBEDDING_FOLDER / "model" / "lnc_doc2vec.model"
        if not model_path.exists():
            train_doc2vec_model_from_fasta(LNCRNA_FASTA, model_path)
        gen_doc2vec_dict_from_fasta(
            LNCRNA_FASTA, doc2vec_dict, model_file=model_path
        )
    if not normal_dict.exists():
        gen_normalized_sequence_dict_from_fasta(LNCRNA_FASTA, normal_dict)
    if not role2vec_dict.exists():
        gen_role2vec_dict(
            load_dict(kmer_dict),
            load_dict(ctd_dict),
            load_dict(doc2vec_dict),
            role2vec_dict,
        )


def generate_pre_mir_seq_characterization():
    folder = MIRNA_FOLDER / "dict" / "precursor"
    kmer_dict = folder / KMER_DICT
    ctd_dict = folder / CTD_DICT
    doc2vec_dict = folder / DOC2VEC_DICT
    normal_dict = folder / NORMAL_SEQ_DICT
    role2vec_dict = folder / ROLE2VEC_DICT

    if not kmer_dict.exists():
        gen_kmer_dict_from_fasta(PRE_MIRNA_FASTA, kmer_dict, 4)
    if not ctd_dict.exists():
        gen_ctd_dict_from_fasta(PRE_MIRNA_FASTA, ctd_dict)
    if not doc2vec_dict.exists():
        model_path = EMBEDDING_FOLDER / "model" / "pre_mir_doc2vec.model"
        if not model_path.exists():
            train_doc2vec_model_from_fasta(PRE_MIRNA_FASTA, model_path)
        gen_doc2vec_dict_from_fasta(
            PRE_MIRNA_FASTA, doc2vec_dict, model_file=model_path
        )
    if not normal_dict.exists():
        gen_normalized_sequence_dict_from_fasta(PRE_MIRNA_FASTA, normal_dict)
    if not role2vec_dict.exists():
        gen_role2vec_dict(
            load_dict(kmer_dict),
            load_dict(ctd_dict),
            load_dict(doc2vec_dict),
            role2vec_dict,
        )


def generate_mat_mir_seq_characterization():
    folder = MIRNA_FOLDER / "dict" / "mature"
    kmer_dict = folder / KMER_DICT
    ctd_dict = folder / CTD_DICT
    doc2vec_dict = folder / DOC2VEC_DICT
    normal_dict = folder / NORMAL_SEQ_DICT
    role2vec_dict = folder / ROLE2VEC_DICT

    if not kmer_dict.exists():
        gen_kmer_dict_from_fasta(MAT_MIRNA_FASTA, kmer_dict, 4)
    if not ctd_dict.exists():
        gen_ctd_dict_from_fasta(MAT_MIRNA_FASTA, ctd_dict)
    if not doc2vec_dict.exists():
        model_path = EMBEDDING_FOLDER / "model" / "mat_mir_doc2vec.model"
        if not model_path.exists():
            train_doc2vec_model_from_fasta(MAT_MIRNA_FASTA, model_path)
        gen_doc2vec_dict_from_fasta(
            MAT_MIRNA_FASTA, doc2vec_dict, model_file=model_path
        )
    if not normal_dict.exists():
        gen_normalized_sequence_dict_from_fasta(MAT_MIRNA_FASTA, normal_dict)
    if not role2vec_dict.exists():
        gen_role2vec_dict(
            load_dict(kmer_dict),
            load_dict(ctd_dict),
            load_dict(doc2vec_dict),
            role2vec_dict,
        )


if __name__ == "__main__":
    generate_lnc_seq_characterization()
    generate_mat_mir_seq_characterization()
    generate_pre_mir_seq_characterization()
