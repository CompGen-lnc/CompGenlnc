from __future__ import annotations

import csv
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from _typeshed import StrPath


def load_dict(
    dict_file: StrPath, val_type: type = np.float32
) -> dict[str, NDArray]:
    """Load a dicationary from a file with pairs of RNA IDs and vectors.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file containing the dictionary.
    val_type : type, optional
        Type of the feature vectors contained in the dictionary.
        Default: `float32`.

    Returns
    -------
    dict of str to numpy.ndarray
        Dictionary containing all vectors associated to their RNA ID.
    """
    dict_file = Path(dict_file).resolve()

    lines = open(dict_file, "r").readlines()
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
    dict_file: StrPath,
    src_dict: dict[str, NDArray] | None = None,
    *,
    keys: Sequence[str] | None = None,
    values: Sequence[NDArray] | None = None,
    overwrite: bool = False,
) -> None:
    """Save a dictionary pairs of RNA names and feature vectors.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the dictionary file.
    src_dict : dict of str to numpy.ndarray, optional
        Dictionary with the feature vectors associated to their RNA ID.
        Should not be `None` if `keys` or `values` are `None`, if it is,
        then `save_dict` will do nothing. Default: `None`.
    keys : sequence of str, optional
        Sequence with the RNA IDs to store. Does not have effect if
        `src_dict` is provided. Default: `None`.
    values : sequence of numpy.ndarray, optional
        Sequence with the feature vectors to store. Must have the same
        length than `keys`. Does not have effect if `src_dict` is
        provided. Default: `None`.
    overwrite : bool, optional
        Replace the existing dictionary if `True` in case there is.
        Otherwise, append the items to the existing dictionary, ignoring
        those that have already been processed. Default: `False`.
    """
    dict_file = Path(dict_file).resolve()
    dict_file.parent.mkdir(parents=True, exist_ok=True)

    if src_dict is not None:
        items = src_dict.items()
    elif keys is not None and values is not None:
        items = zip(keys, values, strict=True)
    else:
        return

    open_mode = "w" if overwrite else "a"
    with open(dict_file, open_mode) as out_file:
        wrapper = csv.writer(out_file)
        for key, value in items:
            if value.ndim == 1:
                value = value.reshape(1, len(value))
            out_file.write(f"{key}\n\t")
            for row in value:
                out_file.write("\t")
                wrapper.writerow(row)
