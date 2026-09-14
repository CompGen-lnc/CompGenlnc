import os
from pathlib import Path
from typing import Iterable

import numpy as np


def load_dict(
    filename: str | os.PathLike, val_type: type = np.float32
) -> dict[str, np.typing.NDArray]:
    """Load a dicationary from a file with pairs of RNA names and feature
    vectors.

    Args:
        filename: Path of the dictionary file.
        val_type: Type of the feature vector.

    Returns:
        dict
    """
    filename = Path(filename).resolve()

    lines = open(filename, "r").readlines()
    res_dict = {}
    actual_id, actual_value = "", []

    for line in lines:
        if not line.strip():
            continue

        if line[0] != "\t":
            if actual_id:
                res_dict[actual_id] = np.array(
                    actual_value, np.float32
                ).squeeze()
            actual_id, actual_value = line.strip(), []
            continue
        value = line.strip().split(',')
        actual_value.append([val_type(x) for x in value])
    res_dict[actual_id] = np.array(actual_value, np.float32).squeeze()

    return res_dict


def save_dict(
    filename: str | os.PathLike,
    dict_: dict[str, Iterable]
) -> None:
    """Save a dicationary pairs of RNA names and feature vectors into a file.

    Args:
        filename: Path of the dictionary file.
        dict_: Dictionary with the RNA names and feature vectors.
    """
    filename = Path(filename).resolve()
    filename.parent.mkdir(parents=True, exist_ok=True)

    keys = list(dict_.keys())
    with open(filename, 'w') as out_file:
        for key in keys:
            out_file.write(f'{key}\n')
            out_file.write(f'\t{','.join(str(x) for x in dict_[key])}\n')
