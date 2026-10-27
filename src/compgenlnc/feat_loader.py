import os
from collections.abc import Collection
from typing import Annotated

import numpy as np
import pandas as pd

from compgenlnc.consts.params import LNCRNA, MATURE, MIRNA, PRECURSOR
from compgenlnc.consts.paths import (
    BINDING_FOLDER,
    CTD_DICT,
    DATASET_FOLDER,
    DOC2VEC_DICT,
    KMER_DICT,
    LOOPS_DICT,
    NORMAL_SEQ_DICT,
    ROLE2VEC_DICT,
)
from compgenlnc.dict_indexer import (
    index_lncrna,
    index_mature_mirna,
    index_precursor_mirna,
    interactions_from_index,
    load_indexes,
)
from compgenlnc.extractors.sequences.normal_sequence import (
    recode_sequence_list,
)
from compgenlnc.extractors.structures2d.adjacency_matrix import (
    get_adjacency_matrix_list,
    normalize_matrix_list,
)
from compgenlnc.extractors.structures2d.recoded_structure import (
    recode_2d_structure_list,
)
from compgenlnc.fileman.miranda_manager import (
    get_subrecord_list,
    get_subrecord_list_from_folder,
    iter_binding_zones,
)
from compgenlnc.typing.numpy_dtypes import Index


def load_feature_set(
    set_path: str | os.PathLike,
    *,
    keep: Collection[Index] | None = None,
    exclude: Collection[Index] | None = None,
) -> np.typing.NDArray[np.float32]:
    feature_set = np.load(set_path)
    mask = np.ones(len(feature_set), np.bool)
    indexes = load_indexes()
    if keep is not None and len(keep):
        mask &= np.isin(indexes, keep)
    if exclude is not None and len(exclude):
        mask &= ~np.isin(indexes, exclude)
    if mask.all():
        return feature_set
    return feature_set[mask]


def get_seq_features(
    prefix: str,
    *,
    from_interactions: bool = False,
    rna: str = LNCRNA,
) -> dict[str, np.typing.NDArray]:
    if not prefix.endswith("_"):
        prefix += "_"

    if rna == LNCRNA:
        records = index_lncrna()
    elif rna == PRECURSOR:
        records = index_precursor_mirna()
    elif rna == MATURE or rna == MIRNA:
        records = index_mature_mirna()
    if from_interactions:
        interactions = interactions_from_index()
        bindings = list(iter_binding_zones(BINDING_FOLDER, interactions))
        records = get_subrecord_list(records, interactions, bindings, rna=rna)

    features_dict = {}
    sequences = [rec.seq for rec in records]

    def get_npy_path(dict_name: str):
        return DATASET_FOLDER / f"{prefix}{dict_name}".replace(".dict", ".npy")

    def map_sequences(func, dtype):
        it = (func(seq) for seq in sequences)
        return np.fromiter(it, dtype, len(sequences))

    features_dict["sequence"] = map_sequences(str, np.dtypes.StringDType)
    features_dict["size"] = map_sequences(len, np.float32)
    if from_interactions:
        it = (binding.energy for binding in bindings)
        features_dict["energy"] = np.fromiter(it, np.float32, len(bindings))
    features = [KMER_DICT, CTD_DICT, DOC2VEC_DICT, ROLE2VEC_DICT]
    for dictname in features:
        feat = dictname.split(".", 1)[0]
        features_dict[feat] = load_feature_set(get_npy_path(dictname))

    if rna != PRECURSOR:
        features_dict["recoded sequence"] = recode_sequence_list(sequences)
        normal_file = get_npy_path(NORMAL_SEQ_DICT)
        features_dict["normal sequence"] = load_feature_set(normal_file)

    return features_dict


def get_ss_features(
    prefix: str,
    *,
    calculate_matrix: bool = False,
    from_interactions: bool = False,
    rna: str = LNCRNA,
) -> dict[str, np.typing.NDArray]:
    if not prefix.endswith("_"):
        prefix += "_"

    if rna == LNCRNA:
        records = index_lncrna()
    elif rna == PRECURSOR or rna == MIRNA:
        records = index_precursor_mirna()
    if from_interactions:
        interactions = interactions_from_index()
        records = get_subrecord_list_from_folder(
            records, interactions, BINDING_FOLDER, rna=rna
        )

    features_dict = {}
    structures = [rec.ss for rec in records]

    def get_npy_path(dict_name: str):
        return DATASET_FOLDER / f"{prefix}{dict_name}".replace(".dict", ".npy")

    def map_structures(func, dtype):
        it = (func(seq) for seq in structures)
        return np.fromiter(it, dtype, len(structures))

    features_dict["structure"] = map_structures(str, np.dtypes.StringDType)
    features_dict["loops"] = load_feature_set(get_npy_path(LOOPS_DICT))
    features_dict["recoded structure"] = recode_2d_structure_list(structures)

    if calculate_matrix and not from_interactions:
        n = len(structures)
        matrix = np.fromiter(get_adjacency_matrix_list(structures), object, n)
        features_dict["adjacency matrix"] = matrix
        features_dict["normal matrix"] = np.array(
            list(normalize_matrix_list(matrix)), np.float32
        )

    return features_dict


type FeatDict5 = Annotated[tuple[tuple[dict, str]], 5]


def get_feature_dictionaries() -> FeatDict5:
    lnc_dict = (
        get_seq_features("lnc", rna=LNCRNA)
        | get_ss_features("lnc", rna=LNCRNA),
        "lnc",
    )
    pre_dict = (
        get_seq_features("pre", rna=PRECURSOR)
        | get_ss_features("pre", rna=PRECURSOR),
        "pre",
    )
    mat_dict = get_seq_features("mat", rna=MATURE), "mat"
    bin_lnc_dict = (
        get_seq_features("binding_lnc", from_interactions=True, rna=LNCRNA)
        | get_ss_features("binding_lnc", from_interactions=True, rna=LNCRNA),
        "binding_lnc",
    )
    bin_mir_dict = (
        get_seq_features("binding_mir", from_interactions=True, rna=MIRNA)
        | get_ss_features("binding_mir", from_interactions=True, rna=MIRNA),
        "binding_mir",
    )

    return (lnc_dict, pre_dict, mat_dict, bin_lnc_dict, bin_mir_dict)


def get_interactions_features() -> pd.DataFrame:
    features_dict = {}

    def add_features(features, prefix):
        for key, value in features.items():
            key = f"{prefix} {key}".replace("_", " ")
            features_dict[key] = value

    for tup in get_feature_dictionaries():
        add_features(*tup)

    return pd.DataFrame.from_dict(features_dict)


def get_feature_vectors():
    dicts, _ = zip(*get_feature_dictionaries())
    feat_names = [
        "size",
        "energy",
        "kmer",
        "ctd",
        "doc2vec",
        "role2vec",
        "normal sequence",
        "loops",
        "normal matrix",
    ]
    features = []

    for feat_dict in dicts:
        for key in feat_names:
            if key not in feat_dict:
                continue
            features.append(feat_dict[key])

    return np.vstack(
        [np.hstack(row, dtype=np.float32) for row in zip(*features)]
    )
