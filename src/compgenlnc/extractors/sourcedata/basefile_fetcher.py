import math
import os
import requests

from compgenlnc.consts.params import (
    GIT_DATA_MIRNA,
    GIT_LNCRNA_FASTA,
    GIT_MAT_MIRNA_FASTA,
    GIT_POSITIVE_PAIRS,
    GIT_PRE_MIRNA_FASTA,
)
from compgenlnc.consts.paths import (
    DATA_MIRNA,
    LNCRNA_FASTA,
    MAT_MIRNA_FASTA,
    PAIRS_FILE_POSITIVE,
    PRE_MIRNA_FASTA,
)
from compgenlnc.fileman.log_manager import send_msg, track_progress

_funits = ["B", "kiB", "MiB", "GiB"]


def fetch_file(url: str, file: str | os.PathLike, file_desc: str = "") -> None:
    send_msg(f"Fetching {file_desc}")
    req = requests.get(url, stream=True)
    if req.status_code != 200:
        send_msg(f"\n{file_desc} not founded\n")
        return
    send_msg(f"\r{file_desc} fetched \n")
    total = int(req.headers.get("Content-Length"))

    exp = math.floor(math.log(total, 1024))
    conversion = 1024**exp
    total /= conversion
    steps = 4096 / conversion

    it = req.iter_content(4096)
    desc = f"Downloading {file_desc}"
    chunks = track_progress(it, desc, total, _funits[exp], steps=steps)
    with open(file, "wb") as out_file:
        for chunk in chunks:
            out_file.write(chunk)
    send_msg(f"{file_desc} downloaded\n")


def fetch_lncrna_files(seq_fasta: str | os.PathLike = LNCRNA_FASTA) -> None:
    fetch_file(GIT_LNCRNA_FASTA, seq_fasta, "lncRNA FASTA")


def fetch_mirna_files(
    mat_seq_fasta: str | os.PathLike = MAT_MIRNA_FASTA,
    pre_seq_fasta: str | os.PathLike = PRE_MIRNA_FASTA,
    precursor_file: str | os.PathLike = DATA_MIRNA,
) -> None:
    fetch_file(GIT_DATA_MIRNA, precursor_file, "Precursor miRNA FASTA")
    fetch_file(GIT_MAT_MIRNA_FASTA, mat_seq_fasta, "Mature miRNA FASTA")
    fetch_file(GIT_PRE_MIRNA_FASTA, pre_seq_fasta, "Precursor association CSV")


def fetch_interaction_files(
    positive_file: str | os.PathLike = PAIRS_FILE_POSITIVE,
) -> None:
    fetch_file(GIT_POSITIVE_PAIRS, positive_file, "Positive pairs CSV")
