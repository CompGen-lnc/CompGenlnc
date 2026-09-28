from collections.abc import Iterable, Iterator

from ViennaRNA import fold_compound

from compgenlnc.structs.seq_record import SeqRecord
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.structs.structure_record import StructureRecord
from compgenlnc.typing.molecules import SeqLike


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
) -> Iterator[StructureRecord]:
    return (extract_2d_structure(seq) for seq in seq_list)


def extract_2d_structure_identified(seq: SeqSSRecord) -> SeqSSRecord:
    return SeqSSRecord(seq.id, seq.seq, extract_2d_structure(seq))
