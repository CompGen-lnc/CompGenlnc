from collections.abc import Iterable, Iterator
from itertools import combinations

import numpy as np
import RNA

from compgenlnc.structs.seq_record import SeqRecord
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.structs.structure_record import StructureRecord
from compgenlnc.typing.molecules import SeqLike


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
        raw_ss = bytearray(b"." * len(seq))
        ss_ranges = []
        substructures = []
        energies = []

        def joint_ss(start, end, sub_ss, energy, data):
            ss_ranges.append(range(start - 1, end))
            substructures.append(sub_ss.encode("ascii"))
            energies.append(energy)

        RNA.Lfold_cb(seq, 300, joint_ss, None)
        energies = np.array(energies, dtype=np.float32)

        length = len(ss_ranges)
        collisions = np.ones((length, length), dtype=np.uint8)
        for i, j in combinations(range(length)):
            ran1 = ss_ranges[i]
            ran2 = ss_ranges[j]
            if ran1.start in ran2 or ran1.stop - 1 in ran2:
                collisions[i, j] = collisions[j, i] = 0

        filtered_indexes = collisions[np.any(collisions == 0, 1)]
        total_energies = [
            np.sum(energies[indexes]) for indexes in filtered_indexes
        ]
        if total_energies.size:
            min_energy_idx = np.argsort(total_energies)[0]
            final_indexes = filtered_indexes[min_energy_idx]
        else:
            final_indexes = range(length)

        for i in final_indexes:
            sub_mfe = substructures[i]
            sub_range = ss_ranges[i]
            raw_ss[sub_range.start : sub_range.stop] = sub_mfe[0]
        ss = raw_ss.decode("ascii")

    return StructureRecord(ss)


def extract_2d_structure_list(
    seq_list: Iterable[SeqLike],
) -> Iterator[StructureRecord]:
    return (extract_2d_structure(seq) for seq in seq_list)


def extract_2d_structure_identified(seq: SeqSSRecord) -> SeqSSRecord:
    return SeqSSRecord(seq.id, seq.seq, extract_2d_structure(seq))
