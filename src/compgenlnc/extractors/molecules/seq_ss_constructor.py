import os
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord

from compgenlnc.config.paths import SEQ_SS_PREFIX
from compgenlnc.extractors.molecules.fasta_constructor import iter_seq_fasta


def join_seq_ss(seq: str, ss: str, id_: str, folder: str | Path) -> None:
    if not isinstance(folder, Path):
        folder = Path(folder).resolve()

    path = folder / f'{SEQ_SS_PREFIX}{id_}.txt'
    with open(path, 'w') as out_file:
        out_file.write(seq + '\n')
        out_file.write(ss)


def fasta_to_seq_ss(filename: str | Path, folder: str | Path) -> int:
    if not isinstance(folder, Path):
        folder = Path(folder).resolve()

    seq_dict = {
        record.id: str(record.seq).replace('U', 'T')
        for record in iter_seq_fasta(filename)
    }
    id_list = list(seq_dict.keys())
    id_list.sort()
    folder.mkdir(parents=True, exist_ok=True)
    total = 0

    for id_ in id_list:
        try:
            seq = seq_dict[id_]
            join_seq_ss(seq, '', id_, folder)
            total += 1
        except:
            print(id_)

    return total
