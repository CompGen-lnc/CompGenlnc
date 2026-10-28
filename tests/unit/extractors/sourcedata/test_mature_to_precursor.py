import pytest

from compgenlnc.extractors.sourcedata import (
    extract_precursor,
    gen_precursors_csv,
    gen_precursors_csv_from_files,
)
from compgenlnc.fileman import load_precursors

from constants import MAT_MIRNA_FASTA, PRE_MIRNA_FASTA


def test_extract_precursor(
    mat_mir_records, pre_mir_records, precursors_consverter, subtests
):
    for mature in mat_mir_records:
        precursor = extract_precursor(mature, pre_mir_records)
        expected = precursors_consverter[mature.id]
        with subtests.test(i=(mature.id, precursor)):
            assert precursor == expected


def test_gen_precursors_csv(
    mat_mir_records, pre_mir_records, precursors_consverter, tmp_path
):
    csv_path = tmp_path / "precursors.csv"
    gen_precursors_csv(csv_path, mat_mir_records, pre_mir_records)
    precursors = load_precursors(csv_path)
    assert precursors == precursors_consverter


def test_gen_precursors_csv_from_files(precursors_consverter, tmp_path):
    csv_path = tmp_path / "precursors.csv"
    gen_precursors_csv_from_files(csv_path, MAT_MIRNA_FASTA, PRE_MIRNA_FASTA)
    precursors = load_precursors(csv_path)
    assert precursors == precursors_consverter
