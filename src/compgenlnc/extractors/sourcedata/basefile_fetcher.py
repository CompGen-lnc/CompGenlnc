from __future__ import annotations

import math
from pathlib import Path
from typing import TYPE_CHECKING

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

if TYPE_CHECKING:
    from _typeshed import StrPath
    from requests import _types as req_types

_funits = ["B", "kiB", "MiB", "GiB"]


def fetch_file(
    url: req_types.UriType, file: StrPath, *, file_desc: str = ""
) -> None:
    """Fetch and download the file from the URL.

    Parameters
    ----------
    url : UriType
        The URL where the file is stored.
    file : str or os.PathLike
        Path where the file will be saved.
    file_desc : str, optional
        A description of the file being fetched. Default: `""`.
    """
    send_msg(f"Fetching {file_desc}")
    req = requests.get(url, stream=True)
    if req.status_code != 200:
        send_msg(f"\n{file_desc} not founded\n")
        return
    send_msg(f"\r{file_desc} fetched \n")
    length = req.headers.get("Content-Length")
    total = int(length) if length else 0

    exp = math.floor(math.log(total, 1024))
    conversion = 1024**exp
    total /= conversion
    steps = 4096 / conversion

    it = req.iter_content(4096)
    desc = f"Downloading {file_desc}"
    chunks = track_progress(it, desc, total, _funits[exp], steps=steps)
    Path(file).parent.mkdir(parents=True, exist_ok=True)
    with open(file, "wb") as out_file:
        for chunk in chunks:
            out_file.write(chunk)
    send_msg(f"{file_desc} downloaded\n")


def fetch_lncrna_files(seq_fasta: StrPath = LNCRNA_FASTA) -> None:
    """Fetch lncRNA FASTA file from the Git repository.

    Parameters
    ----------
    seq_fasta : str or os.PathLike, optional
        Path to store the lncRNA FASTA file. Default:
        `paths.LNCRNA_FASTA`.
    """
    fetch_file(GIT_LNCRNA_FASTA, seq_fasta, file_desc="lncRNA FASTA")


def fetch_mirna_files(
    mat_seq_fasta: StrPath = MAT_MIRNA_FASTA,
    pre_seq_fasta: StrPath = PRE_MIRNA_FASTA,
    precursor_file: StrPath = DATA_MIRNA,
) -> None:
    """Fetch miRNA FASTA and association files from the Git repository.

    Parameters
    ----------
    mat_seq_fasta : str or os.PathLike, optional
        Path to store the mature miRNA FASTA file. Default:
        `paths.MAT_MIRNA_FASTA`.
    pre_seq_fasta : str or os.PathLike, optional
        Path to store the precursor miRNA FASTA file. Default:
        `paths.PRE_MIRNA_FASTA`.
    precursor_file : str or os.PathLike, optional
        Path to store the precursor association CSV file. Default:
        `paths.DATA_MIRNA`.
    """
    fetch_file(
        GIT_DATA_MIRNA, precursor_file, file_desc="Precursor miRNA FASTA"
    )
    fetch_file(
        GIT_MAT_MIRNA_FASTA, mat_seq_fasta, file_desc="Mature miRNA FASTA"
    )
    fetch_file(
        GIT_PRE_MIRNA_FASTA,
        pre_seq_fasta,
        file_desc="Precursor association CSV",
    )


def fetch_interaction_files(
    positive_file: StrPath = PAIRS_FILE_POSITIVE,
) -> None:
    """Fetch the positive pairs CSV file from the GitHub repository.

    Parameters
    ----------
    positive_file : str or os.PathLike, optional
        Path to store the positive pairs CSV file. Default:
        `paths.PAIRS_FILE_POSITIVE`.
    """
    fetch_file(
        GIT_POSITIVE_PAIRS, positive_file, file_desc="Positive pairs CSV"
    )
