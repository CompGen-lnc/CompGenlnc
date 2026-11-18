from __future__ import annotations

import sys
from io import TextIOBase
from pathlib import Path
from typing import TYPE_CHECKING, Final

from numpy._typing import _ArrayLikeFloat_co

if TYPE_CHECKING:
    from _typeshed import StrPath

# Files
log_file: TextIOBase = sys.stdout

# Seeds
doc2vec_seed: _ArrayLikeFloat_co | None = None

# Paths
PACKAGE_DIR: Final[Path] = Path(__file__).resolve().parent
base_path: Path = Path(__file__).resolve().parent

# Environments
conda_environment: str = "base"


def set_base_path(path: StrPath = PACKAGE_DIR) -> None:
    """Set the base path of the data directory.

    Parameters
    ----------
    path : str or os.PathLike, optional
        Path of the new base directory. Default: package directory.
    """
    global base_path
    base_path = Path(path).resolve()


def set_conda_env(env: str = "base") -> None:
    """Set conda virtual environment where miRanda is installed.

    Parameters
    ----------
    env : str, optional
        Conda virtual environment where miRanda is installed. Default:
        `base`.
    """
    global conda_environment
    conda_environment = env


def set_doc2vec_seed(seed: _ArrayLikeFloat_co | None = None) -> None:
    """Set the Doc2Vec seed for the vectors calculation.

    Parameters
    ----------
    seed : number or array of number, optional
        Seed for the Doc2Vec calculations. Default: `None`.
    """
    global doc2vec_seed
    doc2vec_seed = seed


def set_log_file(file: TextIOBase = sys.stdout) -> None:
    """Set the file to which the logs will be written.

    Parameters
    ----------
    file : TextIOBase, optional
        Buffer to which the logs will be written. Default: `sys.stdout`.
    """
    global log_file
    log_file = file
