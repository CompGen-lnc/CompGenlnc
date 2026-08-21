import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord

from compgenlnc.config.paths import *
from compgenlnc.extractors.molecules.fasta_constructor import iter_seq_fasta


def join_seq_ss(seq: str, ss: str, id: str, folder: str | Path) -> None:
    if not isinstance(folder, Path):
        folder = Path(folder).resolve()

    path = MIRNA_SEQ_SS_FOLDER / f'{SEQ_SS_PREFIX}{id}.txt'
    with open(path, 'w') as out_file:
        out_file.write(seq + '\n')
        out_file.write(ss)


def fasta_to_seq_ss(filename: str | Path, folder: str | Path) -> int:
    if not isinstance(folder, Path):
        folder = Path(folder).resolve()

    seq_dict = {
        record.id: record.seq.replace('U', 'T')
        for record in iter_seq_fasta(filename)
    }
    id_list = list(seq_dict.keys())
    folder.mkdir(parents=True, exist_ok=True)
    total = 0

    for id in id_list:
        try:
            seq = seq_dict[id]
            join_seq_ss(seq, '', id)
            total += 1
        except:
            print(id)

    return total
