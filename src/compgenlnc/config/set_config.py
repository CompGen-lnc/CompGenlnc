import os
from io import TextIOBase
from pathlib import Path

from numpy._typing import _ArrayLikeFloat_co

from compgenlnc.config import var_config


def set_data_path(path: str | os.PathLike) -> None:
    var_config.BASE_PATH = Path(path).resolve()


def set_doc2vec_seed(seed: _ArrayLikeFloat_co | None) -> None:
    var_config.DOC2VEC_MODEL_SEED = seed


def set_log_file(file: TextIOBase) -> None:
    var_config.LOG_FILE = file
