import os
import platform
import re
import subprocess
from collections.abc import Iterable
from pathlib import Path

import numpy as np

from compgenlnc.config.paths import BINDING_PREFIX, MIRANDA_PREFIX, TEMP_FOLDER
from compgenlnc.consts.regex import MIRANDA_INFO
from compgenlnc.extractors.characterization.interaction_extractor import load_interactions
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.utils.fasta_manager import load_fasta


def predict_miranda(
    lnc_record: SeqSSRecord, mir_record: SeqSSRecord, folder: os.PathLike
) -> bool:
    folder = Path(folder).resolve()
    folder.mkdir(parents=True, exist_ok=True)

    TEMP_FOLDER.mkdir(parents=True, exist_ok=True)
    lnc_file = TEMP_FOLDER / f"lncrna.fa"
    mir_file = TEMP_FOLDER / f"mirna.fa"
    filename = folder / f"{MIRANDA_PREFIX}{lnc_record.id}_{mir_record.id}.dat"

    def to_wsl(path):
        disk, sub_path = str(path).split(":")
        return "/mnt/" + disk.lower() + sub_path.replace("\\", "/")

    with open(lnc_file, "w") as out_file:
        out_file.write(f">{lnc_record.id}\n")
        out_file.write(f"{lnc_record.seq}\n")
    with open(mir_file, "w") as out_file:
        out_file.write(f">{mir_record.id}\n")
        out_file.write(f"{mir_record.seq}\n")

    command = lambda mir, lnc, out: [
        "conda",
        "run",
        "-n",
        "base",
        "miranda",
        mir,
        lnc,
        "-out",
        out,
    ]

    if platform.system() in ["Linux", "Darwin"]:
        args = command(mir_file, lnc_file, filename)
    elif platform.system() == "Windows":
        cmd = command(to_wsl(mir_file), to_wsl(lnc_file), to_wsl(filename))
        args = ["wsl", "bash", "-i", "-c", " ".join(cmd)]
    else:
        args = []

    success = True
    try:
        subprocess.run(args, check=True)
    except subprocess.CalledProcessError as e:
        print(e)
        success = False

    os.remove(lnc_file)
    os.remove(mir_file)
    return success


def predict_miranda_from_list(
    lnc_list: Iterable[SeqSSRecord],
    mir_list: Iterable[SeqSSRecord],
    folder: os.PathLike,
) -> None:
    records_list = zip(lnc_list, mir_list)
    for records in records_list:
        if None in records:
            continue
        predict_miranda(*records, folder)


def predict_miranda_from_files(
    interaction_file: str | os.PathLike,
    lnc_fasta: str | os.PathLike,
    pre_mir_fasta: str | os.PathLike,
    mat_mir_fasta: str | os.PathLike,
    folder: str | os.PathLike,
) -> None:
    interactions = load_interactions(interaction_file)
    lnc_records = list(load_fasta(lnc_fasta))
    mir_records = list(load_fasta(pre_mir_fasta)) + list(
        load_fasta(mat_mir_fasta)
    )
    lnc_map = {rec.id: rec for rec in lnc_records}
    mir_map = {rec.id: rec for rec in mir_records}
    lnc_list = np.select([interactions["lncRNA"] == key for key in lnc_map.keys()], lnc_map.values(), None)
    mir_list = np.select([interactions["miRNA"] == key for key in mir_map.keys()], mir_map.values(), None)
    predict_miranda_from_list(lnc_list, mir_list, folder)


def extract_binding_zone(
    lnc: str,
    mir: str,
    miranda_folder: os.PathLike,
    binding_folder: os.PathLike,
) -> None:
    miranda_folder = Path(miranda_folder).resolve()
    binding_folder = Path(binding_folder).resolve()

    binding_folder.mkdir(parents=True, exist_ok=True)

    miranda_file = miranda_folder / f"{MIRANDA_PREFIX}{lnc}_{mir}.dat"
    text = miranda_file.read_text()
    hits = list(re.finditer(MIRANDA_INFO, text))
    if hits:
        min_energy_hit = min(hits, key=lambda hit: hit["energy"])
    else:
        min_energy_hit = "No Hits Found above Threshold"
    filename = binding_folder / f"{BINDING_PREFIX}{lnc}_{mir}.dat"
    with open(filename, "w") as out_file:
        out_file.write(min_energy_hit[0])


def extract_binding_zone_from_folder(
    miranda_folder: os.PathLike,
    binding_folder: os.PathLike,
    /,
    lnc_list: Iterable[str] | None = None,
    mir_list: Iterable[str] | None = None,
) -> None:
    miranda_files = os.listdir(miranda_folder)
    for file in miranda_files:
        prefix = file.find(MIRANDA_PREFIX)
        ext = file.find(".dat")
        if prefix < 0 or ext < 0:
            continue

        [lnc, mir] = file[len(MIRANDA_PREFIX) : ext].split("_")
        if (
            lnc_list is not None
            and lnc not in lnc_list
            or mir_list is not None
            and mir not in mir_list
        ):
            continue
        extract_binding_zone(lnc, mir, miranda_folder, binding_folder)
