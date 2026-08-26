from .seq_record import SeqRecord
from dataclasses import dataclass

@dataclass(slots=True)
class SeqSSRecord:
    id: str
    seq: SeqRecord
    ss: str

    def __init__(self, id_: str, seq: str | SeqRecord, ss: str, /):
        self.id = id_
        self.seq = SeqRecord(seq)
        self.ss = ss