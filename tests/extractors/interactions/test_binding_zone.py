import os

from compgenlnc.config.paths import (
    BINDING_PREFIX,
    MIRANDA_PREFIX,
)
from compgenlnc.extractors.interactions import (
    extract_binding_zone,
    extract_binding_zone_from_folder,
    predict_miranda,
    predict_miranda_from_files,
    predict_miranda_from_list,
)
from compgenlnc.fileman import load_interactions

from constants import (
    INTERACTIONS_CSV,
    INTERACTIONS_FOLDER,
    LNCRNA_FASTA,
    MAT_MIRNA_FASTA,
    PRE_MIRNA_FASTA,
)


def test_predict_miranda(
    interaction_set, lnc_records, mir_records, tmp_path, subtests
):
    folder = tmp_path
    for lnc, mir, _ in interaction_set:
        with subtests.test(i=(lnc, mir)):
            lnc_record = next(filter(lambda rec: rec.id == lnc, lnc_records))
            mir_record = next(filter(lambda rec: rec.id == mir, mir_records))
            predict_miranda(lnc_record, mir_record, folder)
            assert (
                folder / f"{MIRANDA_PREFIX}{lnc_record.id}_{mir_record.id}.dat"
            ).exists()


def test_predict_miranda_from_list(
    interaction_set, lnc_records, mir_records, tmp_path
):
    folder = tmp_path
    lnc_list = map(
        lambda lnc: [rec for rec in lnc_records if rec.id == lnc][0],
        interaction_set["lncRNA"],
    )
    mir_list = map(
        lambda mir: [rec for rec in mir_records if rec.id == mir][0],
        interaction_set["miRNA"],
    )
    predict_miranda_from_list(lnc_list, mir_list, folder)
    expected_folder = tmp_path / "expected"
    for lnc, mir, _ in interaction_set:
        lnc_record = next(filter(lambda rec: rec.id == lnc, lnc_records))
        mir_record = next(filter(lambda rec: rec.id == mir, mir_records))
        predict_miranda(lnc_record, mir_record, expected_folder)
        filename = folder / (
            f"{MIRANDA_PREFIX}{lnc_record.id}_{mir_record.id}.dat"
        )
        expected_filename = expected_folder / (
            f"{MIRANDA_PREFIX}{lnc_record.id}_{mir_record.id}.dat"
        )
        assert filename.exists()
        assert filename.read_text() == expected_filename.read_text()


def test_predict_miranda_from_files(interaction_list, tmp_path):
    folder = tmp_path
    predict_miranda_from_files(
        INTERACTIONS_CSV,
        LNCRNA_FASTA,
        PRE_MIRNA_FASTA,
        MAT_MIRNA_FASTA,
        folder,
    )
    for lnc, mir, _ in interaction_list:
        filename = folder / f"{MIRANDA_PREFIX}{lnc}_{mir}.dat"
        assert filename.exists()


def test_extract_binding_zone(interaction_list, tmp_path, subtests):
    miranda_folder = INTERACTIONS_FOLDER / "miranda"
    expected_folder = INTERACTIONS_FOLDER / "binding"
    folder = tmp_path
    for interaction in interaction_list:
        lnc, mir = interaction[["lncRNA", "miRNA"]]
        extract_binding_zone(lnc, mir, miranda_folder, folder)
        with subtests.test(i=(lnc, mir)):
            filename = folder / f"{BINDING_PREFIX}{lnc}_{mir}.dat"
            assert filename.exists()
            expected = expected_folder / f"{BINDING_PREFIX}{lnc}_{mir}.dat"
            assert filename.read_text() == expected.read_text()


def test_extract_binding_zone_from_folder(
    interaction_list, lncrna_id_filter, mirna_id_filter, tmp_path
):
    miranda_folder = INTERACTIONS_FOLDER / "miranda"
    expected_folder = INTERACTIONS_FOLDER / "binding"
    folder = tmp_path
    extract_binding_zone_from_folder(
        miranda_folder, folder, lncrna_id_filter, mirna_id_filter
    )
    for interaction in interaction_list:
        lnc, mir = interaction[["lncRNA", "miRNA"]]
        filename = folder / f"{BINDING_PREFIX}{lnc}_{mir}.dat"
        expected = expected_folder / f"{BINDING_PREFIX}{lnc}_{mir}.dat"

        lnc_filtered = not lncrna_id_filter or lnc not in lncrna_id_filter
        mir_filtered = not mirna_id_filter or mir not in mirna_id_filter
        file_exists = (
            filename.exists() and filename.read_text() == expected.read_text()
        )
        assert (lnc_filtered or mir_filtered) == (not file_exists)
