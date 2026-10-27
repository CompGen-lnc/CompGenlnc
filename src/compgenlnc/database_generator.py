import math
import os

import numpy as np
import pandas as pd

from compgenlnc.consts.fields import (
    INTERACTION_FIELD,
    PAIR_FIELDS,
    PRECURSOR_FIELD,
)
from compgenlnc.consts.params import SEQ_MODE, SS_MODE
from compgenlnc.consts.paths import (
    BINDING_FOLDER,
    BINDING_PREFIX,
    DATA_MIRNA,
    LNCRNA_FASTA,
    LNCRNA_SEQ_SS_FOLDER,
    LNCRNA_SS_FASTA,
    MAT_MIRNA_FASTA,
    MIRANDA_FOLDER,
    MIRANDA_PREFIX,
    MIRNA_SEQ_SS_FOLDER,
    MIRNA_SS_FASTA,
    PAIRS_FILE,
    PAIRS_FILE_NEGATIVE,
    PAIRS_FILE_POSITIVE,
    PRE_MIRNA_FASTA,
)
from compgenlnc.dict_indexer import (
    gen_indexes_from_files,
    index_lncrna,
    index_mature_mirna,
    index_precursor_mirna,
    load_indexes,
)
from compgenlnc.extractors.interactions.binding_zone import (
    extract_binding_zone_from_folder,
    predict_miranda_from_files,
)
from compgenlnc.extractors.sourcedata.fasta_constructor import gen_fasta_2d
from compgenlnc.extractors.sourcedata.interaction_extractor import (
    extract_interactions,
)
from compgenlnc.extractors.sourcedata.mature_to_precursor import (
    gen_precursors_csv_from_files,
)
from compgenlnc.extractors.sourcedata.seq_ss_constructor import fasta_to_seq_ss
from compgenlnc.fileman.csv_manager import load_interactions, load_precursors
from compgenlnc.fileman.fasta_manager import (
    iter_seq_ss,
    load_fasta,
    save_fasta,
)
from compgenlnc.fileman.miranda_manager import read_binding_zone
from compgenlnc.generators.interaction_features import gen_interaction_database
from compgenlnc.generators.seq_features import gen_seq_database
from compgenlnc.generators.ss_features import gen_ss_database


def load_source_data(
    *,
    lnc_seq_fasta: str | os.PathLike | None = None,
    lnc_ss_fasta: str | os.PathLike | None = None,
    pre_mir_seq_fasta: str | os.PathLike | None = None,
    mat_mir_seq_fasta: str | os.PathLike | None = None,
    mir_ss_fasta: str | os.PathLike | None = None,
    data_mirna_csv: str | os.PathLike | None = None,
    interactions_csv: str | os.PathLike | None = None,
    pos_interactions_csv: str | os.PathLike | None = None,
    neg_interactions_csv: str | os.PathLike | None = None,
) -> None:
    if lnc_seq_fasta is not None:
        record_it = load_fasta(lnc_seq_fasta, mode=SEQ_MODE)
        save_fasta(LNCRNA_FASTA, record_it, mode=SEQ_MODE, overwrite=True)
    if lnc_ss_fasta is not None:
        record_it = load_fasta(lnc_ss_fasta, mode=SS_MODE)
        save_fasta(LNCRNA_SS_FASTA, record_it, mode=SS_MODE, overwrite=True)
    if pre_mir_seq_fasta is not None:
        record_it = load_fasta(pre_mir_seq_fasta, mode=SEQ_MODE)
        save_fasta(PRE_MIRNA_FASTA, record_it, mode=SEQ_MODE, overwrite=True)
    if mat_mir_seq_fasta is not None:
        record_it = load_fasta(mat_mir_seq_fasta, mode=SEQ_MODE)
        save_fasta(MAT_MIRNA_FASTA, record_it, mode=SEQ_MODE, overwrite=True)
    if mir_ss_fasta is not None:
        record_it = load_fasta(mir_ss_fasta, mode=SS_MODE)
        save_fasta(MIRNA_SS_FASTA, record_it, mode=SS_MODE, overwrite=True)
    if data_mirna_csv is not None:
        df = pd.read_csv(data_mirna_csv)
        df.to_csv(DATA_MIRNA, index=False)
    if interactions_csv is not None:
        df = pd.read_csv(interactions_csv)
        df.to_csv(PAIRS_FILE, index=False)
    if pos_interactions_csv is not None:
        df = pd.read_csv(pos_interactions_csv)
        df.to_csv(PAIRS_FILE_POSITIVE, index=False)
    if neg_interactions_csv is not None:
        df = pd.read_csv(neg_interactions_csv)
        df.to_csv(PAIRS_FILE_NEGATIVE, index=False)


def gen_source_data(
    filter_indexes: bool = True,
    keep_data: bool = True,
    max_size: int = math.inf,
):
    if not keep_data or not DATA_MIRNA.exists():
        gen_precursors_csv_from_files(
            DATA_MIRNA, MAT_MIRNA_FASTA, PRE_MIRNA_FASTA
        )
        lnc_list = list(
            rec.id
            for rec in load_fasta(LNCRNA_FASTA)
            if len(rec.seq) <= max_size
        )
    mir_list = list(rec.id for rec in load_fasta(MAT_MIRNA_FASTA))
    if not keep_data or not PAIRS_FILE.exists():
        extract_interactions(
            PAIRS_FILE,
            PAIRS_FILE_POSITIVE,
            PAIRS_FILE_NEGATIVE,
            separator=",",
            lnc_list=lnc_list,
            mir_list=mir_list,
        )
    if not keep_data or not PAIRS_FILE.exists():
        gen_indexes_from_files(
            PAIRS_FILE, DATA_MIRNA, lnc_list=lnc_list, mir_list=mir_list
        )

    lnc_keep = pre_mir_keep = mir_keep = []
    if filter_indexes:
        indexes = load_indexes()
        lnc_keep, mir_keep = zip(
            *(inter.split("_") for inter in indexes[INTERACTION_FIELD])
        )
        lnc_keep = set(lnc_keep)
        mir_keep = set(mir_keep)
        pre_mir_keep = set(indexes[PRECURSOR_FIELD])
    gen_fasta_2d(LNCRNA_FASTA, LNCRNA_SS_FASTA, keep=lnc_keep)
    gen_fasta_2d(PRE_MIRNA_FASTA, MIRNA_SS_FASTA, keep=pre_mir_keep)
    fasta_to_seq_ss(
        LNCRNA_FASTA, LNCRNA_SS_FASTA, LNCRNA_SEQ_SS_FOLDER, keep=lnc_keep
    )
    fasta_to_seq_ss(
        PRE_MIRNA_FASTA, MIRNA_SS_FASTA, MIRNA_SEQ_SS_FOLDER, keep=pre_mir_keep
    )

    interactions = load_interactions(PAIRS_FILE)
    miranda_exclude = binding_exclude = []
    str_interactions = ["_".join(inter[PAIR_FIELDS]) for inter in interactions]
    if MIRANDA_FOLDER.exists():
        miranda_files = [
            filename.rsplit(".", 1)[0][len(MIRANDA_PREFIX) :]
            for filename in os.listdir(MIRANDA_FOLDER)
        ]
        mask = np.isin(str_interactions, miranda_files)
        miranda_exclude = interactions[mask]
    if BINDING_FOLDER.exists():
        binding_files = [
            filename[len(BINDING_PREFIX) : -4]
            for filename in os.listdir(BINDING_FOLDER)
        ]
        mask = np.isin(str_interactions, binding_files)
        binding_exclude = interactions[mask]
    predict_miranda_from_files(
        PAIRS_FILE,
        LNCRNA_FASTA,
        PRE_MIRNA_FASTA,
        MAT_MIRNA_FASTA,
        MIRANDA_FOLDER,
        exclude=miranda_exclude,
    )
    extract_binding_zone_from_folder(
        MIRANDA_FOLDER, BINDING_FOLDER, exclude=binding_exclude
    )


def gen_database():
    lnc_records = list(iter_seq_ss(LNCRNA_SEQ_SS_FOLDER))
    pre_mir_records = list(iter_seq_ss(MIRNA_SEQ_SS_FOLDER))
    mat_mir_records = list(load_fasta(MAT_MIRNA_FASTA))
    gen_seq_database(lnc_records, pre_mir_records, mat_mir_records)
    gen_ss_database(lnc_records, pre_mir_records)

    lnc_records = index_lncrna()
    pre_mir_records = index_precursor_mirna()
    mat_mir_records = index_mature_mirna()

    lnc_set = list(rec.id for rec in iter_seq_ss(LNCRNA_SEQ_SS_FOLDER))
    mat_mir_set = list(rec.id for rec in load_fasta(MAT_MIRNA_FASTA))
    interactions = load_interactions(
        PAIRS_FILE, lnc_keep=lnc_set, mir_keep=mat_mir_set
    )

    binding_list = [
        read_binding_zone(BINDING_FOLDER, inter) for inter in interactions
    ]
    precursor_map = load_precursors(DATA_MIRNA)
    gen_interaction_database(
        lnc_records,
        mat_mir_records,
        pre_mir_records,
        interactions,
        precursor_map,
        binding_list,
    )
