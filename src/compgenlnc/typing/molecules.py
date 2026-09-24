from compgenlnc.structs.seq_record import SeqRecord
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.structs.structure_record import StructureRecord

type IDLike = str | SeqSSRecord
type SeqLike = str | SeqRecord | SeqSSRecord
type SSLike = str | StructureRecord | SeqSSRecord
