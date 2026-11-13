from pathlib import Path

import pytest

from phu.avger.reference import AvgReference, load_avg_reference

FIXTURE_AVG = Path(__file__).parent / "fixtures" / "phu_db" / "prepared" / "avg"

FILTER_HEADER = (
    "id\tdb\tname\tfilter_essential\tfilter_glucan\tfilter_lipid\tfilter_methyl\tfilter_nucleotide\n"
)
CONTENT = {
    "positive_amg.tsv": "id\tdb\tname\tmetabolic_paths\ttotal_paths\tpath_source\tmetabolic_ratio\tVL-score\tAMG_weight\tamg_level\nPF00123.17\tPfam\tenergy enzyme\t1\t2\tPfam\t0.5\t2.0\t0.8\tvery high\nK00001\tKEGG\tko enzyme\t1\t1\tKEGG\t1\t1.0\t0.7\thigh\n",
    "positive_apg.tsv": "id\tdb\tname\nPF00123\tPfam\tenergy enzyme\n",
    "positive_areg.tsv": "id\tdb\tname\nK00002\tKEGG\tregulator\n",
    "filters_amg.tsv": FILTER_HEADER + "PF00123\tPfam\tenergy enzyme\tTrue\tFalse\tFalse\tFalse\tFalse\n",
    "filters_apg.tsv": FILTER_HEADER + "K00001\tKEGG\tko enzyme\tFalse\tTrue\tFalse\tFalse\tFalse\n",
    "filters_areg.tsv": FILTER_HEADER + "K00002\tKEGG\tregulator\tFalse\tFalse\tFalse\tFalse\tTrue\n",
    "v_scores.tsv": "id\tname\tV-score\tVL-score\tdb\nPF00123.17\tenergy enzyme\t9.0\t2.0\tPfam\nK00001\tko enzyme\t8.0\t4.0\tKEGG\n",
}


@pytest.fixture
def root(tmp_path: Path) -> Path:
    for name, content in CONTENT.items():
        (tmp_path / name).write_text(content)
    return tmp_path


def test_load_normalises_and_sorts_every_table(root):
    reference = load_avg_reference(root)

    assert reference.scores == [
        {
            "database": "kofam",
            "accession": "K00001",
            "original_accession": "K00001",
            "name": "ko enzyme",
            "v_score": "8.0",
            "vl_score": "4.0",
        },
        {
            "database": "pfam",
            "accession": "PF00123",
            "original_accession": "PF00123.17",
            "name": "energy enzyme",
            "v_score": "9.0",
            "vl_score": "2.0",
        },
    ]
    assert [
        (row["database"], row["accession"], row["proposed_class"]) for row in reference.positive
    ] == [
        ("kofam", "K00001", "putative_amg"),
        ("kofam", "K00002", "putative_areg"),
        ("pfam", "PF00123", "putative_amg"),
        ("pfam", "PF00123", "putative_apg"),
    ]
    amg = reference.positive[2]
    assert amg["amg_weight"] == "0.8"
    assert amg["amg_level"] == "very high"
    assert [
        (row["database"], row["accession"], row["proposed_class"], row["filter_category"])
        for row in reference.filters
    ] == [
        ("kofam", "K00001", "putative_apg", "filter_glucan"),
        ("kofam", "K00002", "putative_areg", "filter_nucleotide"),
        ("pfam", "PF00123", "putative_amg", "filter_essential"),
    ]


def test_filter_class_comes_from_the_file_it_is_in(root):
    # The three filter tables share a schema and carry no class column.
    (root / "filters_apg.tsv").write_text(
        FILTER_HEADER + "PF00123\tPfam\tenergy enzyme\tTrue\tFalse\tFalse\tFalse\tFalse\n"
    )

    classes = {row["proposed_class"] for row in load_avg_reference(root).filters}

    assert classes == {"putative_amg", "putative_apg", "putative_areg"}


def test_rows_from_other_databases_are_skipped(root):
    (root / "v_scores.tsv").write_text(CONTENT["v_scores.tsv"] + "X1\tother\t1\t1\tCAMPER\n")

    assert len(load_avg_reference(root).scores) == 2


def test_conflicting_duplicate_scores_fail(root):
    (root / "v_scores.tsv").write_text(
        "id\tname\tV-score\tVL-score\tdb\n"
        "PF00123\tfirst\t9\t2\tPfam\n"
        "PF00123\tsecond\t9\t2\tPfam\n"
    )

    with pytest.raises(ValueError, match="Conflicting AVG score rows"):
        load_avg_reference(root)


def test_missing_required_column_fails(root):
    (root / "positive_apg.tsv").write_text("id\tname\nPF00123\tenergy enzyme\n")

    with pytest.raises(ValueError, match="missing columns: db"):
        load_avg_reference(root)


def test_fixture_bundle_avg_database_loads():
    reference = load_avg_reference(FIXTURE_AVG)

    assert isinstance(reference, AvgReference)
    assert reference.scores
    assert {row["proposed_class"] for row in reference.positive} == {
        "putative_amg",
        "putative_apg",
        "putative_areg",
    }


def test_phrog_rows_are_kept_with_their_own_identifiers():
    reference = load_avg_reference(FIXTURE_AVG)

    scores = {(row["database"], row["accession"]): row for row in reference.scores}
    assert scores[("phrog", "phrog_1")]["v_score"] == "10.0"
    positive = {
        (row["database"], row["accession"], row["proposed_class"]) for row in reference.positive
    }
    assert ("phrog", "phrog_2", "putative_apg") in positive
    assert ("phrog", "phrog_3", "putative_areg") in positive
    assert not any(key[0] == "phrog" and key[2] == "putative_amg" for key in positive)
    assert {
        (row["accession"], row["proposed_class"], row["filter_category"])
        for row in reference.filters
        if row["database"] == "phrog"
    } == {
        ("phrog_1", "putative_amg", "filter_essential"),
        ("phrog_1", "putative_apg", "filter_essential"),
        ("phrog_1", "putative_areg", "filter_essential"),
    }
