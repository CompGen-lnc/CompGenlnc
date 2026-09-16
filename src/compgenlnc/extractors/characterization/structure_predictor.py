from collections.abc import Iterable

from ViennaRNA import fold_compound

from compgenlnc.structs import SeqRecord, SeqSSRecord, StructureRecord
from compgenlnc.typing import SeqLike


def extract_2d_structure(seq: SeqLike) -> StructureRecord:
    if isinstance(seq, SeqSSRecord):
        seq = str(seq.seq)
    elif isinstance(seq, SeqRecord):
        seq = str(seq)

    fc = fold_compound(seq)
    (ss, _) = fc.mfe()
    return StructureRecord(ss)


def extract_2d_structure_list(
    seq_list: Iterable[SeqLike],
) -> Iterable[StructureRecord]:
    return (extract_2d_structure(seq) for seq in seq_list)


def extract_2d_structure_identified(seq: SeqSSRecord) -> SeqSSRecord:
    return SeqSSRecord(seq.id, seq.seq, extract_2d_structure(seq))
