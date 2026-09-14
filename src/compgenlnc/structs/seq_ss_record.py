from dataclasses import dataclass

from .seq_record import SeqRecord
from .structure_record import StructureRecord


@dataclass(slots=True)
class SeqSSRecord:
    id: str
    seq: SeqRecord
    ss: StructureRecord

    def __init__(
        self, id_: str, seq: str | SeqRecord, ss: str | StructureRecord, /
    ):
        self.id = id_
        self.seq = SeqRecord(seq)
        self.ss = StructureRecord(ss)
