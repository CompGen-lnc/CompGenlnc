from __future__ import annotations

import os
import textwrap
from collections.abc import Collection, Iterable, Iterator
from pathlib import Path
from typing import TYPE_CHECKING

from compgenlnc.consts.params import SEQ_MODE, SS_MODE
from compgenlnc.consts.paths import (
    SEQ_SS_PREFIX,
)
from compgenlnc.structs.seq_ss_record import SeqSSRecord

if TYPE_CHECKING:
    from _typeshed import StrPath


def iter_fasta(
    filename: StrPath,
    *,
    mode: str = SEQ_MODE,
    keep: Collection[str] | None = None,
    exclude: Collection[str] | None = None,
) -> Iterator[SeqSSRecord]:
    """Load a FASTA file and iterate the RNA molecules into them.

    Parameters
    ----------
    filename : str or os.PathLike
        Path of the FASTA file.
    mode : str
        Type of the data saved in the file. If `'seq'`, it will treat
        the data as sequences. If `'ss'`, it will treat the data as 2D
        structures.
    keep : collection of str, optional
        Set of the RNA IDs to save. Do not apply if `None` or empty.
        Default: `None`.
    exclude : collection of str, optional
        Set of the RNA IDs to ignore. Do not apply if `None` or empty.
        Default: `None`.

    Yields
    ------
    SeqSSRecord
        A record with the sequence or 2D structure with its ID.

    Raises
    ------
    ValueError
        If `mode` is any other than `'seq'` or `'ss'`.
    """
    if mode not in [SEQ_MODE, SS_MODE]:
        raise ValueError(f"Value {mode} for mode is not allowed")

    if keep is None:
        keep = []
    if exclude is None:
        exclude = []

    filename = Path(filename).resolve()
    id_ = seq = ss = value = ""
    with open(filename) as file:
        # Add the '>' at the end for capturing the last molecule
        for raw_line in file.readlines() + [">"]:
            if raw_line[0] != ">":
                if id_ and raw_line.strip():
                    value += raw_line.split()[0]
                continue

            if mode == SEQ_MODE:
                seq = value.upper().replace("T", "U")
            elif mode == SS_MODE:
                ss = value

            if id_:
                yield SeqSSRecord(id_, seq, ss)
                id_ = ""
                value = ""

            if not raw_line[1:].strip():
                continue

            line = raw_line[1:].split()[0]
            if len(keep) and line not in keep or line in exclude:
                continue
            id_ = line


def load_fasta(filename: StrPath, **kwargs) -> list[SeqSSRecord]:
    """Load a FASTA file and list the RNA molecules into them.

    Parameters
    ----------
    filename : str or os.PathLike
        Path of the FASTA file.
    mode : str
        Type of the data saved in the file. If `'seq'`, it will treat
        the data as sequences. If `'ss'`, it will treat the data as 2D
        structures.
    keep : collection of str, optional
        Set of the RNA IDs to save. Do not apply if `None` or empty.
        Default: `None`.
    exclude : collection of str, optional
        Set of the RNA IDs to ignore. Do not apply if `None` or empty.
        Default: `None`.

    Returns
    -------
    list of SeqSSRecord
        A list of records with the sequence or 2D structure with their
        ID.

    Raises
    ------
    ValueError
        If `mode` is any other than `'seq'` or `'ss'`.
    """
    records = list(iter_fasta(filename, **kwargs))
    return records


def save_fasta(
    filename: StrPath,
    records: Iterable[SeqSSRecord],
    *,
    mode: str,
    overwrite: bool = False,
) -> None:
    """Save sequences or 2D structures from a set of RNA molecules.

    Parameters
    ----------
    filename : str or os.PathLike
        Path of the FASTA file.
    records : iterable object of SeqSSRecord
        Group of RNA molecules to save.
    mode : str
        Type of the data saved in the file. If `'seq'`, it will save the
        sequences. If `'ss'`, it will save the 2Dstructures.
    overwrite : bool, optional
        Append the items to the existing dictionary if `True`, ignoring
        those that have already been processed. Otherwise, replace the
        existing dictionary. Default: `False`.

    Raises
    ------
    ValueError
        If `mode` is any other than `'seq'` or `'ss'`.
    """
    if mode not in [SEQ_MODE, SS_MODE]:
        raise ValueError(f"Value {mode} for mode is not allowed")

    filename = Path(filename).resolve()
    filename.parent.mkdir(parents=True, exist_ok=True)
    open_mode = "w" if overwrite else "a"
    with open(filename, open_mode) as out_file:
        for record in records:
            value = record.ss if mode == SS_MODE else record.seq
            out_file.write(f">{record.id}\n")
            if not value:
                continue
            out_file.write(f"{textwrap.fill(str(value), 60)}\n")


def iter_seq_fasta(
    filename: StrPath, filter: Collection[str] | None = None
) -> Iterator[SeqSSRecord]:
    """Iterate all the RNA sequences from a FASTA file.

    Parameters
    ----------
    filename : str or os.PathLike
        Path of the FASTA file.
    filter : collection of str
        A set with all substrings contained in the ID of the desired
        sequences. If None, it is the same as using `load_fasta` with
        mode as `'seq'`.

    Yields
    ------
    SeqSSRecord
        A record with the sequence with its ID.
    """
    filename = Path(filename).resolve()
    if filter is None:
        filter = []

    for record in load_fasta(filename, mode='seq'):
        if len(filter) and not any(prefix in record.id for prefix in filter):
            continue
        yield record


def read_sequence(dir: StrPath, id_: str) -> str:
    """Read the first line of a seq+ss file.

    Parameters
    ----------
    dir : str or os.PathLike
        Path of the directory that contains the seq+ss file.
    id_ : str
        ID of the RNA to read the sequence

    Returns
    -------
    str
        Sequence of the RNA searched.
    """
    dir = Path(dir).resolve()
    filename = dir / f"{SEQ_SS_PREFIX}{id_}.dat"

    with open(filename) as f:
        seq = f.readlines()[0].strip()
    return seq.upper().replace("T", "U")


def read_structure(dir: StrPath, id_: str) -> str:
    """Read the second line of a seq+ss file.

    Parameters
    ----------
    dir : str or os.PathLike
        Path of the directory that contains the seq+ss file.
    id_ : str
        ID of the RNA to read the 2D structure

    Returns
    -------
    str
        2D structure of the RNA searched.
    """
    dir = Path(dir).resolve()
    filename = dir / f"{SEQ_SS_PREFIX}{id_}.dat"

    with open(filename) as f:
        ss = f.readlines()[1].strip()
    return ss


def iter_seq_ss(
    dir: StrPath,
    *,
    keep: Collection[str] | None = None,
    exclude: Collection[str] | None = None,
) -> Iterator[SeqSSRecord]:
    """Iterate all the seq+ss files from the same directory.

    Parameters
    ----------
    dir : str or os.PathLike
        Path of the directory that contains the seq+ss file.
    keep : collection of str, optional
        Set of the RNA IDs to save. Do not apply if `None` or empty.
        Default: `None`.
    exclude : collection of str, optional
        Set of the RNA IDs to ignore. Do not apply if `None` or empty.
        Default: `None`.

    Yields
    ------
    SeqSSRecord
        A record with the sequence and 2D structure with its ID.
    """
    dir = Path(dir).resolve()
    if keep is None:
        keep = []
    if exclude is None:
        exclude = []

    for filename in sorted(os.listdir(dir)):
        if not filename.startswith(SEQ_SS_PREFIX):
            continue
        id_ = os.path.splitext(filename[len(SEQ_SS_PREFIX) :])[0]
        if len(keep) and id_ not in keep or id_ in exclude:
            continue
        with open(dir / filename) as file:
            seq = file.readline().upper().replace("T", "U")
            ss = file.readline()
        yield SeqSSRecord(id_, seq, ss)


def load_seq_ss(dir: StrPath, **kwargs) -> list[SeqSSRecord]:
    """List all the seq+ss files from the same directory.

    Parameters
    ----------
    dir : str or os.PathLike
        Path of the directory that contains the seq+ss file.
    keep : collection of str, optional
        Set of the RNA IDs to save. Do not apply if `None` or empty.
        Default: `None`.
    exclude : collection of str, optional
        Set of the RNA IDs to ignore. Do not apply if `None` or empty.
        Default: `None`.

    Returns
    -------
    list of SeqSSRecord
        A list of records with the sequence and 2D structure with their ID.
    """
    records = list(iter_seq_ss(dir, **kwargs))
    return records
