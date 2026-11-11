from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

from compgenlnc import config

SEQ_SS_PREFIX = "seq+ss_join_"
LOOPS_PREFIX = "loops_"
MIRANDA_PREFIX = "miranda_"
BINDING_PREFIX = "binding_"


class DyPath(Path):
    """A path that is relative to `config.base_path`."""

    def __new__(cls, relative_subpath=""):
        """Create a new `DyPath` instance.

        Parameters
        ----------
        relative_subpath : str, optional
            The relative path to the file or directory. Default: `""`.
        """
        return super().__new__(cls)

    def __init__(self, relative_subpath=""):
        """Initialize the `DyPath` object.

        Parameters
        ----------
        relative_subpath : str, optional
            The relative path to the file or directory. Default: `""`.
        """
        self._relative_subpath = str(relative_subpath)

    @property
    def _real_path(self) -> Path:
        """Get the absolute path of the file or directory.

        Returns
        -------
        pathlib.Path
            The absolute path.
        """
        return config.base_path / self._relative_subpath

    def __str__(self) -> str:
        """Get the string representation of the path.

        Returns
        -------
        str
            The string representation.
        """
        return str(self._real_path)

    def __repr__(self) -> str:
        """Get the printable representation of the object.

        Returns
        -------
        str
            The printable representation.
        """
        return repr(self._real_path)

    def __fspath__(self) -> str:
        """Get the file system path representation of the object.

        Returns
        -------
        str
            The file system path representation.
        """
        return str(self._real_path)

    def __truediv__(self, key) -> DyPath:
        """Get the path of given key relative to the current path.

        Parameters
        ----------
        key : str
            The key to be used to get the path.

        Returns
        -------
        DyPath
            The new path.
        """
        return DyPath(Path(self._relative_subpath) / key)

    def __getattr__(self, name) -> Any:
        """Provide attribute access to the underlying path object.

        Parameters
        ----------
        name : str
            Name of the attribute to get.

        Returns
        -------
        Any
            The value of the attribute.
        """
        return getattr(self._real_path, name)

    def __reduce__(self) -> tuple[Callable, tuple[str]]:
        """Reduce the path object to a tuple for pickling.

        Returns
        -------
        func : callable
            The function to call to create a new path object.
        args : tuple of str
            The arguments to pass to the function.
        """
        return Path, (str(self._real_path),)


# Main data dirs
DATA_DIR = DyPath("data")
EMBEDDING_DIR = DATA_DIR / "embedding"
INTERACTIONS_DIR = DATA_DIR / "interactions"
LNCRNA_DIR = DATA_DIR / "lncrna"
MIRNA_DIR = DATA_DIR / "mirna"
TEMP_DIR = DATA_DIR / "temp"

# Embedding
INDEX_NPY = EMBEDDING_DIR / "index.npy"
DATASET_DIR = EMBEDDING_DIR / "dataset"

MODELS_DIR = EMBEDDING_DIR / "models"
LNCRNA_DOC2VEC_MODEL = MODELS_DIR / "lnc_doc2vec.model"
PRE_MIRNA_DOC2VEC_MODEL = MODELS_DIR / "pre_mir_doc2vec.model"
MAT_MIRNA_DOC2VEC_MODEL = MODELS_DIR / "mat_mir_doc2vec.model"
LNCRNA_INTERACTION_DOV2VEC_MODEL = MODELS_DIR / "binding_lnc_doc2vec.model"
MIRNA_INTERACTION_DOV2VEC_MODEL = MODELS_DIR / "binding_mir_doc2vec.model"

# Interactions
PAIRS_FILE = INTERACTIONS_DIR / "interactions.csv"
PAIRS_FILE_NEGATIVE = INTERACTIONS_DIR / "pairs" / "negative_pairs.csv"
PAIRS_FILE_POSITIVE = INTERACTIONS_DIR / "pairs" / "positive_pairs.csv"

BINDING_DIR = INTERACTIONS_DIR / "binding"
MIRANDA_DIR = INTERACTIONS_DIR / "miranda"

# lncRNA
LNCRNA_FASTA = LNCRNA_DIR / "fasta" / "lncRNA.fa"
LNCRNA_SS_FASTA = LNCRNA_DIR / "fasta" / "lncrna_2d.fa"

LNCRNA_SEQ_SS_DIR = LNCRNA_DIR / "characterization" / "seq+ss_join"
LNCRNA_LOOPS_DIR = LNCRNA_DIR / "characterization" / "ss_loops"
LNCRNA_DICTS_DIR = LNCRNA_DIR / "dicts"

# miRNA
DATA_MIRNA = MIRNA_DIR / "data_miRNA.csv"
PRE_MIRNA_FASTA = MIRNA_DIR / "fasta" / "miRNA.fa"
MAT_MIRNA_FASTA = MIRNA_DIR / "fasta" / "mature_miRNA.fa"
MIRNA_SS_FASTA = MIRNA_DIR / "fasta" / "miRNA_2d.fa"

MIRNA_SEQ_SS_DIR = MIRNA_DIR / "characterization" / "seq+ss_join"
MIRNA_LOOPS_DIR = MIRNA_DIR / "characterization" / "ss_loops"
MIRNA_DICTS_DIR = MIRNA_DIR / "dicts"
PRE_MIRNA_DICTS_DIR = MIRNA_DICTS_DIR / "precursor"
MAT_MIRNA_DICTS_DIR = MIRNA_DICTS_DIR / "mature"

# Dictionaries
KMER_DICT = "kmer.dict"
CTD_DICT = "ctd.dict"
DOC2VEC_DICT = "doc2vec.dict"
ROLE2VEC_DICT = "role2vec.dict"
NORMAL_SEQ_DICT = "normal_seq.dict"

LOOPS_DICT = "loops.dict"
NORMAL_MAT_DICT = "normal_matrix.dict"
