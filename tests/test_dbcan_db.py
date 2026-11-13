import json
from pathlib import Path

import pytest

from phu.db import dbcan as dbcan_db
from phu.db.dbcan import (
    build_family_offsets,
    is_dbcan_id,
    normalize_dbcan_id,
    parse_pul_rules,
)

FIXTURE_DBCAN = Path(__file__).parent / "fixtures" / "phu_db" / "prepared" / "dbcan"
HEADER = "pul_id\tsubstrate\tcazymes_predicted_dbcan\n"


def _rules(tmp_path: Path, rows: str, family_ids: set[str]):
    path = tmp_path / "pul_rules.tsv"
    path.write_text(HEADER + rows)
    return parse_pul_rules(path, family_ids)


@pytest.fixture
def dbcan_root(tmp_path, monkeypatch) -> Path:
    """A dbCAN database directory laid out like the bundle's."""
    root = tmp_path / "dbcan"
    root.mkdir()
    (root / "profiles.hmm").write_text(
        "HMMER3/f\nNAME  GH5.hmm\n//\n"
        "HMMER3/f\nNAME  CBM48.hmm\n//\n"
        "HMMER3/f\nNAME  SLH.hmm\n//\n"
    )
    (root / "pul_rules.tsv").write_text(
        HEADER + "PUL0001\txylan\tGH5,CBM48\nPUL0002\tstarch\tGH5|SLH\n"
    )
    monkeypatch.setattr(dbcan_db, "_dbcan_root", lambda: root)
    return root


def test_dbcan_family_names_normalize_hmm_suffix():
    assert is_dbcan_id("gh128.hmm")
    assert normalize_dbcan_id("gh128.hmm") == "GH128"


def test_pul0621_rule_keeps_token_order(tmp_path: Path):
    hmm = tmp_path / "dbcan.txt"
    hmm.write_text(
        "HMMER3/f\nNAME GH133.hmm\n//\n"
        "HMMER3/f\nNAME GH3.hmm\n//\n"
        "HMMER3/f\nNAME GH57.hmm\n//\n"
        "HMMER3/f\nNAME GT4.hmm\n//\n"
    )

    rules = _rules(
        tmp_path, "PUL0621\tstarch\tGH133|GH3,GH57|GT4\n", set(build_family_offsets(hmm))
    )

    assert rules["PUL0621"].families == ("GH133", "GH3", "GH57", "GT4")
    assert rules["PUL0621"].raw_rule == "GH133|GH3,GH57|GT4"
    assert rules["PUL0621"].substrate == "starch"


def test_pul_rule_preserves_unsupported_tokens(tmp_path: Path):
    rule = _rules(tmp_path, "PUL0001\t\tPL6_3,CBM35inCE17\n", set())["PUL0001"]

    assert rule.families == ()
    assert rule.unresolved_tokens == ("PL6_3", "CBM35INCE17")
    assert not rule.resolved


def test_pul_rule_marks_ancillary_profiles_unresolved(tmp_path: Path):
    rule = _rules(tmp_path, "PUL0187\t\tGH16_21,SLH|CBM54\n", {"GH16_21", "CBM54"})["PUL0187"]

    assert rule.families == ("GH16_21", "CBM54")
    assert rule.unresolved_tokens == ("SLH",)
    assert not rule.resolved


def test_pul_listed_twice_keeps_first_position_and_last_rule(tmp_path: Path):
    rules = _rules(
        tmp_path, "PUL0001\ta\tGH5\nPUL0002\tb\tGH5\nPUL0001\tc\tCBM48\n", {"GH5", "CBM48"}
    )

    assert list(rules) == ["PUL0001", "PUL0002"]
    assert rules["PUL0001"].substrate == "c"
    assert rules["PUL0001"].families == ("CBM48",)


def test_pul_rule_with_empty_rule_is_unresolved(tmp_path: Path):
    rule = _rules(tmp_path, "PUL0009\tpectin\t\n", {"GH5"})["PUL0009"]

    assert rule.raw_rule == ""
    assert not rule.resolved


def test_unexpected_pul_rules_header_is_refused(tmp_path: Path):
    path = tmp_path / "pul_rules.tsv"
    path.write_text("id\tsubstrate\trule\n")

    with pytest.raises(ValueError, match="unexpected PUL rules header"):
        parse_pul_rules(path, set())


def test_fixture_bundle_pul_rules_parse():
    rules = parse_pul_rules(FIXTURE_DBCAN / "pul_rules.tsv", {"GH5", "CBM48"})

    assert list(rules) == ["PUL0001", "PUL0002", "PUL0003"]
    assert rules["PUL0002"].substrate == "starch glycogen"
    assert rules["PUL0003"].unresolved_tokens == ("GH999",)


def test_inventory_and_rules_use_canonical_models_only(dbcan_root):
    canonical, excluded = dbcan_db.get_dbcan_model_inventory()

    assert canonical == ("GH5", "CBM48")
    assert excluded == ("SLH",)
    assert dbcan_db.get_dbcan_pul("pul0001").resolved
    assert dbcan_db.get_dbcan_pul("PUL0002").unresolved_tokens == ("SLH",)
    assert dbcan_db.get_dbcan_pul("PUL9999") is None
    assert list(dbcan_db.get_dbcan_pul_rules()) == ["PUL0001", "PUL0002"]


def test_extract_writes_models_and_reports_missing(dbcan_root, tmp_path):
    paths, missing = dbcan_db.extract_dbcan_models(["gh5", "GH999"], tmp_path / "out")

    assert [path.name for path in paths] == ["GH5.hmm"]
    assert paths[0].read_text() == "HMMER3/f\nNAME  GH5.hmm\n//\n"
    assert missing == ["GH999"]


def test_offset_index_is_rebuilt_when_profiles_change(dbcan_root, tmp_path):
    dbcan_db.extract_dbcan_models(["GH5"], tmp_path / "out")
    index = json.loads((dbcan_root / "family_offsets.json").read_text())
    assert set(index["offsets"]) == {"GH5", "CBM48", "SLH"}

    # A different profiles.hmm: the index and the cached GH5 model are stale.
    (dbcan_root / "profiles.hmm").write_text(
        "HMMER3/f\nNAME  GH5.hmm\nDESC  changed\n//\nHMMER3/f\nNAME  GH9.hmm\n//\n"
    )
    paths, missing = dbcan_db.extract_dbcan_models(["GH5", "GH9", "CBM48"], tmp_path / "out2")

    assert "changed" in paths[0].read_text()
    assert [path.name for path in paths] == ["GH5.hmm", "GH9.hmm"]
    assert missing == ["CBM48"]
