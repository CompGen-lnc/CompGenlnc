import pytest

from compgenlnc.extractors.molecules import (
    extract_precursor,
    gen_precursors_csv,
    gen_precursors_csv_from_files,
)
from compgenlnc.fileman import load_precursors

from constants import MAT_MIRNA_FASTA, PRE_MIRNA_FASTA


def test_extract_precursor(mat_mirna_id_list, pre_mirna_id_list, precursors_consverter, subtests):
    for mature in mat_mirna_id_list:
        precursor = extract_precursor(mature, pre_mirna_id_list)
        expected = precursors_consverter[mature]
        with subtests.test(i=(mature, precursor)):
            assert precursor == expected


def test_gen_precursors_csv(mat_mirna_id_list, pre_mirna_id_list, precursors_consverter, tmp_path):
    csv_path = tmp_path / "precursors.csv"
    gen_precursors_csv(csv_path, mat_mirna_id_list, pre_mirna_id_list)
    precursors = load_precursors(csv_path)
    assert precursors == precursors_consverter


def test_gen_precursors_csv_from_files(precursors_consverter, tmp_path):
    csv_path = tmp_path / "precursors.csv"
    gen_precursors_csv_from_files(csv_path, MAT_MIRNA_FASTA, PRE_MIRNA_FASTA)
    precursors = load_precursors(csv_path)
    assert precursors == precursors_consverter
