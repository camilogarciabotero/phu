from pathlib import Path

import pytest

from phu import screen
from phu.dbcan_db import PULRule


def test_resolve_hmm_inputs_resolves_dbcan_families(monkeypatch, tmp_path: Path):
    def fake_ensure(*, force_refresh=False):
        return {"hmm_path": str(tmp_path / "dbCAN-HMMdb-V15.txt")}

    def fake_extract(requested_ids, output_dir):
        assert requested_ids == ["GH128", "CBM89"]
        output_dir.mkdir(parents=True, exist_ok=True)
        paths = []
        for family in requested_ids:
            path = output_dir / f"{family}.hmm"
            path.write_text("HMMER3/f\n//\n")
            paths.append(path)
        return paths, []

    monkeypatch.setattr(screen, "ensure_dbcan_database", fake_ensure)
    monkeypatch.setattr(screen, "extract_dbcan_models", fake_extract)

    resolved, metadata = screen._resolve_hmm_inputs(
        [Path("gh128"), Path("CBM89")], tmp_path / "out"
    )

    assert [path.name for path in resolved] == ["GH128.hmm", "CBM89.hmm"]
    assert metadata == {}


def test_resolve_hmm_inputs_expands_pul(monkeypatch, tmp_path: Path):
    pul = PULRule(
        pul_id="PUL0621",
        substrate="starch",
        raw_rule="GH133,GH3,GH57,GT4",
        families=("GH133", "GH3", "GH57", "GT4"),
        unresolved_tokens=(),
    )

    monkeypatch.setattr(screen, "ensure_dbcan_database", lambda **_: {})
    monkeypatch.setattr(screen, "get_dbcan_pul", lambda _: pul)

    def fake_extract(requested_ids, output_dir):
        assert requested_ids == list(pul.families)
        output_dir.mkdir(parents=True, exist_ok=True)
        paths = []
        for family in requested_ids:
            path = output_dir / f"{family}.hmm"
            path.write_text("HMMER3/f\n//\n")
            paths.append(path)
        return paths, []

    monkeypatch.setattr(screen, "extract_dbcan_models", fake_extract)

    resolved, metadata = screen._resolve_hmm_inputs([Path("pul0621")], tmp_path / "out")

    assert [path.stem for path in resolved] == list(pul.families)
    assert metadata == {}


def test_resolve_hmm_inputs_rejects_pul_mixed_with_other_inputs(tmp_path: Path):
    with pytest.raises(ValueError, match="only database input"):
        screen._resolve_hmm_inputs([Path("PUL0621"), Path("GH128")], tmp_path / "out")


def test_resolve_hmm_inputs_rejects_unresolved_pul(monkeypatch, tmp_path: Path):
    pul = PULRule("PUL0001", "", "PL6_3", (), ("PL6_3",))
    monkeypatch.setattr(screen, "ensure_dbcan_database", lambda **_: {})
    monkeypatch.setattr(screen, "get_dbcan_pul", lambda _: pul)

    with pytest.raises(ValueError, match="unresolved dbCAN rule token"):
        screen._resolve_hmm_inputs([Path("PUL0001")], tmp_path / "out")


@pytest.mark.parametrize(
    ("evalue", "coverage", "kept"),
    [
        (1e-15, 0.9, False),
        (1e-16, 0.35, False),
        (1e-16, 0.3501, True),
    ],
)
def test_dbcan_thresholds_use_strict_boundaries(evalue, coverage, kept):
    hit = screen.Hit(
        contig="c1",
        prot_id="c1|gene1",
        model="GH128",
        bitscore=1.0,
        evalue=evalue,
        hmm_coverage=coverage,
    )

    kept_hits, contigs = screen._choose_best_contigs(
        [hit],
        min_bitscore=None,
        max_evalue=1e-5,
        queried_model_ids=["GH128"],
        dbcan_model_ids=["GH128"],
    )

    assert bool(kept_hits) is kept
    assert bool(contigs) is kept