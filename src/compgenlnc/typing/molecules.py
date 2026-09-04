from compgenlnc.structs import SeqRecord, SeqSSRecord, StructureRecord


type SeqLike = str | SeqRecord | SeqSSRecord
type SSLike = str | StructureRecord | SeqSSRecord
