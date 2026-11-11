from dataclasses import dataclass

from .seq_record import SeqRecord
from .structure_record import StructureRecord


@dataclass(slots=True)
class SeqSSRecord:
    """Record for the RNA sequence and 2D structure.

    Attributes
    ----------
    id : str
        The ID of the RNA.
    seq : SeqRecord
        The sequence of the RNA.
    ss : StructureRecord
        The 2D structure of the RNA.
    """

    id: str
    seq: SeqRecord
    ss: StructureRecord

    def __init__(
        self, id_: str, seq: str | SeqRecord, ss: str | StructureRecord, /
    ) -> None:
        """Initialize the SeqSSRecord object.

        Parameters
        ----------
        id_ : str
            The ID of the RNA.
        seq : SeqRecord
            The sequence of the RNA.
        ss : StructureRecord
            The 2D structure of the RNA.
        """
        self.id = id_
        self.seq = SeqRecord(seq)
        self.ss = StructureRecord(ss)
