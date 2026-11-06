from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from _typeshed import StrPath


def load_dict(
    filename: StrPath, val_type: type = np.float32
) -> dict[str, NDArray]:
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

    for line in lines + ["none"]:
        if not line.strip():
            continue

        if line[0] != "\t":
            if actual_id:
                res_dict[actual_id] = np.array(
                    actual_value, val_type
                ).squeeze()
            actual_id, actual_value = line.strip(), []
            continue
        value = line.strip().split(",")
        actual_value.append([val_type(x) for x in value])

    return res_dict


def save_dict(
    filename: StrPath,
    dict_: dict[str, NDArray] | None = None,
    *,
    keys: Iterable[str] | None = None,
    values: Iterable[NDArray] | None = None,
    overwrite: bool = False,
) -> None:
    """Save a dicationary pairs of RNA names and feature vectors into a file.

    Args:
        filename: Path of the dictionary file.
        dict_: Dictionary with the RNA names and feature vectors.
    """
    filename = Path(filename).resolve()
    filename.parent.mkdir(parents=True, exist_ok=True)

    if dict_ is not None:
        items = dict_.items()
    elif keys is not None and values is not None:
        items = zip(keys, values, strict=True)
    else:
        return

    open_mode = "w" if overwrite else "a"
    with open(filename, open_mode) as out_file:
        wrapper = csv.writer(out_file)
        for key, value in items:
            if value.ndim == 1:
                value = value.reshape(1, len(value))
            out_file.write(f"{key}\n\t")
            for row in value:
                out_file.write("\t")
                wrapper.writerow(row)
