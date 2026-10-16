import io
import math
from collections import defaultdict
from collections.abc import Iterable, Iterator

import numpy as np
import RNA

from compgenlnc.structs.seq_record import SeqRecord
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.structs.structure_record import StructureRecord
from compgenlnc.typing.molecules import SeqLike


def extract_massive_2d_structure(seq: SeqLike) -> str:
    ss = io.StringIO("." * len(seq))
    sub_ss_list = []
    events = []
    length = 0

    def joint_ss(start, end, sub_ss, energy, data):
        nonlocal length
        sub_ss_list.append([sub_ss, range(start - 1, end), energy])
        events.append((start - 1, True, length))
        events.append((end, False, length))
        length += 1

    RNA.Lfold_cb(seq, 300, joint_ss, None)
    sub_ss_list = np.array(sub_ss_list, dtype=object)

    events.sort()
    active = set()
    first_appear, last_appear = {}, {}
    total_sets = 0
    ss_started = False

    for _, start, i in events:
        if start:
            active.add(i)
            ss_started = True
        else:
            if ss_started:
                for j in active:
                    if j not in first_appear.keys():
                        first_appear[j] = total_sets
                    last_appear[j] = total_sets
                total_sets += 1
                ss_started = False
            active.discard(i)

    set_beginnings = defaultdict(list)
    for index, collision_set in first_appear.items():
        set_beginnings[collision_set].append(index)

    def indexes_from(k):
        if k == total_sets:
            yield []
        else:
            for i in set_beginnings[k]:
                for r in indexes_from(last_appear[i] + 1):
                    yield [*r, i]

    final_indexes = None
    min_energy = math.inf
    for indexes in indexes_from(0):
        if not indexes:
            continue
        total_energy = np.sum(sub_ss_list[indexes, 2])
        if total_energy < min_energy:
            min_energy = total_energy
            final_indexes = indexes

    if final_indexes:
        sub_ss_list = sub_ss_list[final_indexes]

    for sub_ss, sub_range in sub_ss_list[:, :2]:
        ss.seek(sub_range.start)
        ss.write(sub_ss)
    return ss.getvalue()


def extract_2d_structure(seq: SeqLike) -> StructureRecord:
    if isinstance(seq, SeqSSRecord):
        seq = str(seq.seq)
    elif isinstance(seq, SeqRecord):
        seq = str(seq)
    if not seq:
        return StructureRecord("")

    if len(seq) < 3500:
        fc = RNA.fold_compound(seq)
        ss, _ = fc.mfe()
    else:
        ss = extract_massive_2d_structure(seq)

    return StructureRecord(ss)


def extract_2d_structure_list(
    seq_list: Iterable[SeqLike],
) -> Iterator[StructureRecord]:
    return (extract_2d_structure(seq) for seq in seq_list)


def extract_2d_structure_identified(seq: SeqSSRecord) -> SeqSSRecord:
    return SeqSSRecord(seq.id, seq.seq, extract_2d_structure(seq))
