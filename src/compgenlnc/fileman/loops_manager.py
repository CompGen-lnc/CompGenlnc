import os
from collections.abc import Iterable, Iterator
from pathlib import Path

from compgenlnc.config.paths import LOOPS_PREFIX
from compgenlnc.consts.regex import EXTERNAL_LOOP_RE, LOOP_RE
from compgenlnc.fileman.fasta_manager import iter_seq_ss
from compgenlnc.structs.loop_counter import LoopCounter


def get_2d_struture_loops(
    loops_file: str | os.PathLike,
) -> LoopCounter:
    def match_loop(line: str) -> tuple[bool, tuple[str, int, int, int]]:
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
    return LoopCounter(loops, energy)


def get_2d_structure_loops_from_folder(
    loops_folder: str | os.PathLike,
    /,
    keep: Iterable[str] | None = None,
) -> Iterator[LoopCounter]:
    loops_folder = Path(loops_folder).resolve()
    for filename in os.listdir(loops_folder):
        id_ = filename.split(".")[0][len(LOOPS_PREFIX) :]
        if (
            keep is not None
            and id_ not in keep
            or not filename.startswith(LOOPS_PREFIX)
        ):
            continue
        filename = loops_folder / filename
        yield get_2d_struture_loops(filename)
