import random

import compgenlnc as cgl


def test_predict_miranda(
    interaction_set, lnc_records, mir_records, tmp_path, subtests
):
    folder = tmp_path
    count = len(interaction_set)
    samples = random.sample(range(count), min(10, count))
    for lnc, mir, _ in interaction_set[samples]:
        with subtests.test(i=(lnc, mir)):
            lnc_record = next(filter(lambda rec: rec.id == lnc, lnc_records))
            mir_record = next(filter(lambda rec: rec.id == mir, mir_records))
            cgl.predict_miranda(lnc_record, mir_record, folder)
            filename = (
                f"{cgl.paths.MIRANDA_PREFIX}"
                f"{lnc_record.id}_{mir_record.id}.dat"
            )
            assert (folder / filename).exists()


def test_predict_miranda_from_list(
    interaction_set, lnc_records, mir_records, tmp_path
):
    count = len(interaction_set)
    samples = random.sample(range(count), min(10, count))
    interactions = interaction_set[samples]

    folder = tmp_path
    lnc_list = map(
        lambda lnc: next(rec for rec in lnc_records if rec.id == lnc),
        interactions[cgl.LNCRNA_FIELD],
    )
    mir_list = map(
        lambda mir: next(rec for rec in mir_records if rec.id == mir),
        interactions[cgl.MIRNA_FIELD],
    )
    cgl.predict_miranda_from_list(lnc_list, mir_list, folder)
    expected_folder = tmp_path / "expected"
    for lnc, mir, _ in interactions:
        lnc_record = next(rec for rec in lnc_records if rec.id == lnc)
        mir_record = next(rec for rec in mir_records if rec.id == mir)
        cgl.predict_miranda(lnc_record, mir_record, expected_folder)
        filename = folder / (
            f"{cgl.paths.MIRANDA_PREFIX}{lnc_record.id}_{mir_record.id}.dat"
        )
        expected = expected_folder / (
            f"{cgl.paths.MIRANDA_PREFIX}{lnc_record.id}_{mir_record.id}.dat"
        )
        assert filename.exists()
        file_text, expected_text = filename.read_text(), expected.read_text()
        file_idx = file_text.find("Read Sequence")
        expected_idx = expected_text.find("Read Sequence")
        assert file_text[file_idx:] == expected_text[expected_idx:]


def test_predict_miranda_from_files(interaction_list, tmp_path):
    count = len(interaction_list)
    samples = random.sample(range(count), min(10, count))
    interactions = interaction_list[samples]

    folder = tmp_path
    cgl.predict_miranda_from_files(
        cgl.paths.PAIRS_FILE,
        cgl.paths.LNCRNA_FASTA,
        cgl.paths.PRE_MIRNA_FASTA,
        cgl.paths.MAT_MIRNA_FASTA,
        folder,
        keep=interactions,
    )
    for lnc, mir, _ in interactions:
        filename = folder / f"{cgl.paths.MIRANDA_PREFIX}{lnc}_{mir}.dat"
        assert filename.exists()


def test_extract_binding_zone(interaction_list, tmp_path, subtests):
    miranda_folder = cgl.paths.INTERACTIONS_FOLDER / "miranda"
    expected_folder = cgl.paths.INTERACTIONS_FOLDER / "binding"
    folder = tmp_path
    for interaction in interaction_list:
        lnc, mir = interaction[cgl.PAIR_FIELDS]
        cgl.extract_binding_zone(lnc, mir, miranda_folder, folder)
        with subtests.test(i=(lnc, mir)):
            filename = folder / f"{cgl.paths.BINDING_PREFIX}{lnc}_{mir}.dat"
            assert filename.exists()
            expected = (
                expected_folder / f"{cgl.paths.BINDING_PREFIX}{lnc}_{mir}.dat"
            )
            assert filename.read_text() == expected.read_text()


def test_extract_binding_zone_from_folder(
    interaction_list, lncrna_id_filter, mirna_id_filter, tmp_path
):
    miranda_folder = cgl.paths.MIRANDA_FOLDER
    expected_folder = cgl.paths.BINDING_FOLDER
    folder = tmp_path
    cgl.extract_binding_zone_from_folder(
        miranda_folder, folder, lncrna_id_filter, mirna_id_filter
    )
    for interaction in interaction_list:
        lnc, mir = interaction[cgl.PAIR_FIELDS]
        filename = folder / f"{cgl.paths.BINDING_PREFIX}{lnc}_{mir}.dat"
        expected = (
            expected_folder / f"{cgl.paths.BINDING_PREFIX}{lnc}_{mir}.dat"
        )

        lnc_included = not lncrna_id_filter or lnc in lncrna_id_filter
        mir_included = not mirna_id_filter or mir in mirna_id_filter
        file_exists = (
            filename.exists() and filename.read_text() == expected.read_text()
        )
        assert (lnc_included and mir_included) == file_exists
