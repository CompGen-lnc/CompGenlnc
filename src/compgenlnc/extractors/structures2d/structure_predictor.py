from collections.abc import Iterable, Iterator

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

    raw_ss = bytearray(b"." * len(seq))

    def joint_ss(start, end, sub_ss, energy, data):
        idx_start = start - 1
        for i, caracter in enumerate(sub_ss.encode("ascii")):
            if caracter == 40 or caracter == 41:
                pos = idx_start + i
                raw_ss[pos] = caracter
        return 0

    if len(seq) < 3500:
        fc = RNA.fold_compound(seq)
        (ss, _) = fc.mfe()
    else:
        RNA.Lfold_cb(seq, 300, joint_ss, None)
        ss = raw_ss.decode("ascii")

    return StructureRecord(ss)


def extract_2d_structure_list(
    seq_list: Iterable[SeqLike],
) -> Iterator[StructureRecord]:
    return (extract_2d_structure(seq) for seq in seq_list)


def extract_2d_structure_identified(seq: SeqSSRecord) -> SeqSSRecord:
    return SeqSSRecord(seq.id, seq.seq, extract_2d_structure(seq))
