import os
import sys
from io import TextIOBase
from pathlib import Path
from typing import Final

from numpy._typing import _ArrayLikeFloat_co

# Files
LOG_FILE = sys.stdout

# Seeds
DOC2VEC_MODEL_SEED = 0

# Paths
BASE_PATH: Final[Path] = Path(__file__).resolve().parent


def set_base_path(path: str | os.PathLike) -> None:
    global BASE_PATH
    BASE_PATH = Path(path).resolve()


def set_doc2vec_seed(seed: _ArrayLikeFloat_co | None) -> None:
    global DOC2VEC_MODEL_SEED
    DOC2VEC_MODEL_SEED = seed


def set_log_file(file: TextIOBase) -> None:
    global LOG_FILE
    LOG_FILE = file
