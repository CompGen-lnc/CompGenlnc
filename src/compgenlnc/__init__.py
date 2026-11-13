import importlib.util
from warnings import warn

from . import consts, extractors, fileman, generators, structs
from .config import (
    set_base_path,
    set_doc2vec_seed,
    set_log_file,
)
from .consts import *
from .database_generator import gen_database, gen_source_data, load_source_data
from .dict_indexer import (
    dict_to_numpy,
    gen_indexes,
    get_mature_index,
    index_lncrna,
    index_mature_mirna,
    index_precursor_mirna,
    interactions_from_index,
    load_indexes,
)
from .dtypes import (
    dtype_index,
    dtype_interaction_pair,
    dtype_interaction_tuple,
    dtype_loop_tuple,
)
from .extractors import *
from .feat_loader import (
    get_feature_vectors,
    get_interaction_features,
    get_seq_features,
    get_ss_features,
)
from .fileman import *
from .generators import *
from .version import __version__

HAS_PYG = importlib.util.find_spec("pyg_lib") is not None
if not HAS_PYG:
    warn(
        "pyg_lib not installed, some functions will not work. Please install "
        "it manually or run `compgenlnc-install-deps` instead",
        stacklevel=2,
    )

__all__ = [
    "dict_to_numpy",
    "dtype_index",
    "dtype_interaction_pair",
    "dtype_interaction_tuple",
    "dtype_loop_tuple",
    "get_feature_vectors",
    "get_interaction_features",
    "get_mature_index",
    "get_seq_features",
    "get_ss_features",
    "gen_database",
    "gen_indexes",
    "gen_source_data",
    "index_lncrna",
    "index_mature_mirna",
    "index_precursor_mirna",
    "interactions_from_index",
    "load_indexes",
    "load_source_data",
    "set_base_path",
    "set_doc2vec_seed",
    "set_log_file",
]
__all__.extend(consts.__all__)
__all__.extend(extractors.__all__)
__all__.extend(fileman.__all__)
__all__.extend(generators.__all__)
__all__.extend(structs.__all__)
