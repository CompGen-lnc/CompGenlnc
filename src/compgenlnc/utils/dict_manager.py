import os
from pathlib import Path
from typing import Iterable

import numpy as np

def load_dict(
    filename: str | os.PathLike,
    val_type: type = np.float64
) -> dict[str, np.typing.NDArray]:
    filename = Path(filename).resolve()

    lines = open(filename, "r").readlines()
    res_dict = {}
    actual_id, actual_value = '', []

    for line in lines:
        if line[0] != '\t':
            res_dict[actual_id] = np.array(actual_value).squeeze()
            actual_id, actual_value = line.strip(), []
            print(actual_id)
            continue

        if not line.strip():
            continue
        value = line.strip().split(',')
        actual_value.append([val_type(x) for x in value])
    res_dict[actual_id] = np.array(actual_value).squeeze()

    return res_dict


def save_dict(
    filename: str | os.PathLike,
    dict_: dict[str, Iterable]
) -> None:
    filename = Path(filename).resolve()

    keys = list(dict_.keys())
    with open(filename, 'w') as out_file:
        for key in keys:
            out_file.write(f'{key}\n')
            out_file.write(f'\t{','.join(str(x) for x in dict_[key])}\n')
