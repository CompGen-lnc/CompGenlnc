from __future__ import annotations

import math
import os
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from compgenlnc.consts.fields import (
    INTERACTION_FIELD,
    PAIR_FIELDS,
    PRECURSOR_FIELD,
)
from compgenlnc.consts.params import SEQ_MODE, SS_MODE
from compgenlnc.consts.paths import (
    BINDING_DIR,
    BINDING_PREFIX,
    DATA_MIRNA,
    INDEX_NPY,
    LNCRNA_FASTA,
    LNCRNA_SEQ_SS_DIR,
    LNCRNA_SS_FASTA,
    MAT_MIRNA_FASTA,
    MIRANDA_DIR,
    MIRANDA_PREFIX,
    MIRNA_SEQ_SS_DIR,
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
    extract_binding_zone_from_dir,
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
    iter_fasta,
    iter_seq_ss,
    load_fasta,
    load_seq_ss,
    save_fasta,
)
from compgenlnc.fileman.miranda_manager import read_binding_zone
from compgenlnc.generators.interaction_features import gen_interaction_database
from compgenlnc.generators.seq_features import gen_seq_database
from compgenlnc.generators.ss_features import gen_ss_database

if TYPE_CHECKING:
    from _typeshed import StrPath


def load_source_data(
    base_dir: StrPath | None = None,
    *,
    lnc_seq_fasta: StrPath | None = None,
    lnc_ss_fasta: StrPath | None = None,
    pre_mir_seq_fasta: StrPath | None = None,
    mat_mir_seq_fasta: StrPath | None = None,
    mir_ss_fasta: StrPath | None = None,
    data_mirna_csv: StrPath | None = None,
    interactions_csv: StrPath | None = None,
    pos_interactions_csv: StrPath | None = None,
    neg_interactions_csv: StrPath | None = None,
) -> None:
    """Load the base files to the data directory.

    Parameters
    ----------
    base_dir : str or os.PathLike, optional
        Path of the directory containing the base files to load. These
        files do not overwrite paths defined below. Default: `None`.
    lnc_seq_fasta : str or os.PathLike, optional
        Path of the FASTA file with the lncRNA sequences.
        Default: `None`.
    lnc_ss_fasta : str or os.PathLike, optional
        Path of the FASTA file with the lncRNA 2D structures.
        Default: `None`.
    pre_mir_seq_fasta : str or os.PathLike, optional
        Path of the FASTA file with the precursor miRNA sequences.
        Default: `None`.
    mat_mirna_seq_fasta : str or os.PathLike, optional
        Path of the FASTA file with the mature miRNA sequences.
        Default: `None`.
    mirna_ss_fasta : str or os.PathLike, optional
        Path of the FASTA file with the precursor miRNA 2D structures.
        Default: `None`.
    data_mirna_csv : str or os.PathLike, optional
        Path of the CSV file with the mature-precursor miRNA pairs.
        Default: `None`.
    interactions_csv : str or os.PathLike, optional
        Path of the CSV file with the interaction tuples. Default:
        `None`.
    positive_interactions_csv : str or os.PathLike, optional
        Path of the CSV file with the positive interaction pairs.
        Default: `None`.
    negative_interactions_csv : str or os.PathLike, optional
        Path of the CSV file with the negative interaction pairs.
        Default: `None`.
    """
    if base_dir is not None:
        base_dir = Path(base_dir).resolve()

        def get_path(name, default):
            """Get the path from `base_dir` if not specified."""
            if default is not None:
                return default
            return base_dir / name

        lnc_seq_fasta = get_path(LNCRNA_FASTA.name, lnc_seq_fasta)
        lnc_ss_fasta = get_path(LNCRNA_SS_FASTA.name, lnc_ss_fasta)
        pre_mir_seq_fasta = get_path(PRE_MIRNA_FASTA.name, pre_mir_seq_fasta)
        mat_mir_seq_fasta = get_path(MAT_MIRNA_FASTA.name, mat_mir_seq_fasta)
        mir_ss_fasta = get_path(MIRNA_SS_FASTA.name, mir_ss_fasta)
        data_mirna_csv = get_path(DATA_MIRNA.name, data_mirna_csv)
        interactions_csv = get_path(PAIRS_FILE.name, interactions_csv)
        pos_interactions_csv = get_path(
            PAIRS_FILE_POSITIVE.name, pos_interactions_csv
        )
        neg_interactions_csv = get_path(
            PAIRS_FILE_NEGATIVE.name, neg_interactions_csv
        )

    def load_fasta_if_exists(old_fasta, new_fasta, mode):
        """Load a FASTA file if is specified."""
        if old_fasta is not None and Path(old_fasta).resolve().exists():
            records = iter_fasta(old_fasta, mode=mode)
            save_fasta(new_fasta, records, mode=mode, overwrite=True)

    load_fasta_if_exists(lnc_seq_fasta, LNCRNA_FASTA, mode=SEQ_MODE)
    load_fasta_if_exists(lnc_ss_fasta, LNCRNA_SS_FASTA, mode=SS_MODE)
    load_fasta_if_exists(pre_mir_seq_fasta, PRE_MIRNA_FASTA, mode=SEQ_MODE)
    load_fasta_if_exists(mat_mir_seq_fasta, MAT_MIRNA_FASTA, mode=SEQ_MODE)
    load_fasta_if_exists(mir_ss_fasta, MIRNA_SS_FASTA, mode=SS_MODE)

    def load_csv_if_exists(old_csv, new_csv):
        """Load a CSV file if is specified."""
        if old_csv is not None and Path(old_csv).resolve().exists():
            df = pd.read_csv(old_csv)
            Path(new_csv).resolve().parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(new_csv, index=False)

    load_csv_if_exists(data_mirna_csv, DATA_MIRNA)
    load_csv_if_exists(interactions_csv, PAIRS_FILE)
    load_csv_if_exists(pos_interactions_csv, PAIRS_FILE_POSITIVE)
    load_csv_if_exists(neg_interactions_csv, PAIRS_FILE_NEGATIVE)


def gen_source_data(
    filter_indexes: bool = True,
    keep_data: bool = True,
    max_size: int = math.inf,
) -> None:
    """Generate the base files that were not specified.

    Parameters
    ----------
    filter_indexes : bool, optional
        Whether the data should contain only the RNAs and interactions
        in the indexes. Default: `True`.
    keep_data : bool, optional
        Whether the data calculated should not be recalculated. Default:
        `True`.
    max_size : int, optional
        The maximum size of RNAs to process. Default: infinite.
    """
    # Generate the CSV files that were not specified
    if not keep_data or not DATA_MIRNA.exists():
        gen_precursors_csv_from_files(
            DATA_MIRNA, MAT_MIRNA_FASTA, PRE_MIRNA_FASTA
        )
    lnc_list = [
        rec.id for rec in iter_fasta(LNCRNA_FASTA) if len(rec.seq) <= max_size
    ]
    mir_list = [rec.id for rec in iter_fasta(MAT_MIRNA_FASTA)]
    if not keep_data or not PAIRS_FILE.exists():
        extract_interactions(
            PAIRS_FILE,
            PAIRS_FILE_POSITIVE,
            PAIRS_FILE_NEGATIVE,
            separator=",",
            lnc_list=lnc_list,
            mir_list=mir_list,
        )

    # Generate the indexes
    if not keep_data or not INDEX_NPY.exists():
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

    # Generate the 2D structure FASTA and seq+ss files
    gen_fasta_2d(LNCRNA_FASTA, LNCRNA_SS_FASTA, keep=lnc_keep)
    gen_fasta_2d(PRE_MIRNA_FASTA, MIRNA_SS_FASTA, keep=pre_mir_keep)
    fasta_to_seq_ss(
        LNCRNA_FASTA, LNCRNA_SS_FASTA, LNCRNA_SEQ_SS_DIR, keep=lnc_keep
    )
    fasta_to_seq_ss(
        PRE_MIRNA_FASTA, MIRNA_SS_FASTA, MIRNA_SEQ_SS_DIR, keep=pre_mir_keep
    )

    # Generate the binding zone files
    interactions = load_interactions(PAIRS_FILE)
    miranda_exclude = binding_exclude = []
    str_interactions = ["_".join(inter[PAIR_FIELDS]) for inter in interactions]
    if MIRANDA_DIR.exists():
        miranda_files = [
            filename.rsplit(".", 1)[0][len(MIRANDA_PREFIX) :]
            for filename in os.listdir(MIRANDA_DIR)
        ]
        mask = np.isin(str_interactions, miranda_files)
        miranda_exclude = interactions[mask]
    if BINDING_DIR.exists():
        binding_files = [
            filename[len(BINDING_PREFIX) : -4]
            for filename in os.listdir(BINDING_DIR)
        ]
        mask = np.isin(str_interactions, binding_files)
        binding_exclude = interactions[mask]
    predict_miranda_from_files(
        PAIRS_FILE,
        LNCRNA_FASTA,
        PRE_MIRNA_FASTA,
        MAT_MIRNA_FASTA,
        MIRANDA_DIR,
        exclude=miranda_exclude,
    )
    extract_binding_zone_from_dir(
        MIRANDA_DIR, BINDING_DIR, exclude=binding_exclude
    )


def gen_database() -> None:
    """Generate the lncRNA-miRNA interaction database."""
    # Generate the sequence and 2D structure features
    lnc_records = load_seq_ss(LNCRNA_SEQ_SS_DIR)
    pre_mir_records = load_seq_ss(MIRNA_SEQ_SS_DIR)
    mat_mir_records = load_fasta(MAT_MIRNA_FASTA)
    gen_seq_database(lnc_records, pre_mir_records, mat_mir_records)
    gen_ss_database(lnc_records, pre_mir_records)

    lnc_records = index_lncrna()
    pre_mir_records = index_precursor_mirna()
    mat_mir_records = index_mature_mirna()

    # Generate the interaction features
    lnc_list = [rec.id for rec in iter_seq_ss(LNCRNA_SEQ_SS_DIR)]
    mat_mir_list = [rec.id for rec in iter_fasta(MAT_MIRNA_FASTA)]
    interactions = load_interactions(
        PAIRS_FILE, lnc_keep=lnc_list, mir_keep=mat_mir_list
    )

    binding_list = [
        read_binding_zone(BINDING_DIR, inter) for inter in interactions
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
