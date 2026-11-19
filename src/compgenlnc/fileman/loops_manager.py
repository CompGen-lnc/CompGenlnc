from __future__ import annotations

import os
from collections.abc import Collection, Iterator
from pathlib import Path
from typing import TYPE_CHECKING

from compgenlnc.consts.paths import LOOPS_PREFIX
from compgenlnc.consts.regex import EXTERNAL_LOOP_RE, LOOP_RE
from compgenlnc.structs.loop_counter import LoopCounter

if TYPE_CHECKING:
    from _typeshed import StrPath


def read_2d_structure_loops(dir: StrPath, ss_id: str) -> LoopCounter:
    """Read the 2D structure loops of a RNA.

    Parameters
    ----------
    dir : str or os.PathLike
        The directory containing the loops file.
    ss_id : str
        The ID of the RNA.

    Returns
    -------
    LoopCounter
        The loops for the 2D structure of the RNA.
    """

    def match_loop(line: str) -> tuple[bool, tuple[str, int, int, int]]:
        """Match a loop line in the loops file."""
        m = LOOP_RE.match(line)
        if m:
            return True, (
                m.group(1).split()[0],
                *(int(x) for x in m.group(2, 3, 4)),
            )
        m = EXTERNAL_LOOP_RE.match(line)
        if m:
            return True, (m.group(1).split()[0], 0, 0, int(m.group(2)))
        return False, ()

    loops_file = Path(dir).resolve() / f"{LOOPS_PREFIX}{ss_id}.dat"
    energy = 0
    loops = []
    with open(loops_file) as file:
        for line in file.readlines():
            if not line.strip():
                continue
            is_loop, loop = match_loop(line)
            words = line.split()
            if is_loop:
                loops.append(loop)
            elif words[0] == "Energy":
                energy = float(words[-1])
    return LoopCounter(ss_id, loops, energy)


def read_2d_structure_loops_from_dir(
    loops_dir: StrPath,
    *,
    keep: Collection[str] | None = None,
    exclude: Collection[str] | None = None,
) -> Iterator[LoopCounter]:
    """Read the 2D structure loops from loops files.

    Parameters
    ----------
    loops_dir : str or os.PathLike
        The directory containing the loops files.
    keep : collection of str, optional
        Set of IDs of the RNAs to save. Do not apply if `None` or empty.
        Default: `None`.
    exclude : collection of str, optional
        Set of IDs of the RNAs to ignore. Do not apply if `None` or
        empty. Default: `None`.

    Yields
    ------
    LoopCounter
        The loops for the 2D structure of the RNA.
    """
    for filename in sorted(os.listdir(loops_dir)):
        id_ = filename[len(LOOPS_PREFIX) : -4]
        filtered = keep and id_ not in keep or exclude and id_ in exclude
        if filtered:
            continue
        yield read_2d_structure_loops(loops_dir, id_)
