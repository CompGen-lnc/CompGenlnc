from typing import Iterable

from ViennaRNA import fold

from compgenlnc.structs import SeqRecord, SeqSSRecord, StructureRecord
from compgenlnc.typing import SeqLike


def extract_2d_structure(seq: SeqLike) -> StructureRecord:
    if isinstance(seq, SeqSSRecord):
        seq = str(seq.seq)
    elif isinstance(seq, SeqRecord):
        seq = str(seq)

    (ss, _) = fold(seq)
    return StructureRecord(ss)


def extract_2d_structure_list(
        seq_list: Iterable[SeqLike]
) -> Iterable[StructureRecord]:
    return (extract_2d_structure(seq) for seq in seq_list)


def extract_2d_structure_identified(seq: SeqSSRecord) -> SeqSSRecord:
    return SeqSSRecord(seq.id, seq.seq, extract_2d_structure(seq))
    