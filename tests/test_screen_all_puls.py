import csv
import gzip
from pathlib import Path

from typer.testing import CliRunner

from phu import screen
from phu.cli import app
from phu.dbcan_db import PULRule

runner = CliRunner()


def _rule(pul_id: str, families: tuple[str, ...], substrate: str = "") -> PULRule:
    return PULRule(pul_id, substrate, ",".join(families), families, ())


def test_evaluate_pul_signatures_uses_independent_or_and_semantics():
    rules = (
        _rule("PUL0001", ("GH1", "GH2")),
        _rule("PUL0002", ("GH2", "GH3")),
        _rule("PUL0003", ("GH1", "GH2")),
    )
    hits = [
        screen.Hit("contig_A", "a|gene1", "GH1", 1, 1e-20),
        screen.Hit("contig_A", "a|gene2", "GH2", 1, 1e-20),
        screen.Hit("contig_B", "b|gene1", "GH2", 1, 1e-20),
        screen.Hit("contig_B", "b|gene2", "GH3", 1, 1e-20),
        screen.Hit("contig_C", "c|gene1", "GH1", 1, 1e-20),
        screen.Hit("contig_C", "c|gene2", "GH3", 1, 1e-20),
    ]

    matches, contigs = screen.evaluate_pul_signatures(hits, rules)

    assert [(match.contig_id, match.pul_id) for match in matches] == [
        ("contig_A", "PUL0001"),
        ("contig_A", "PUL0003"),
        ("contig_B", "PUL0002"),
    ]
    assert contigs == ["contig_A", "contig_B"]


def test_evaluate_pul_signatures_does_not_match_missing_required_family():
    matches, contigs = screen.evaluate_pul_signatures(
        [screen.Hit("c1", "c1|gene1", "GH1", 1, 1e-20)],
        [_rule("PUL0001", ("GH1", "GH2"))],
    )
    assert matches == []
    assert contigs == []


def test_evaluate_pul_signatures_preserves_identical_signatures():
    rules = (_rule("PUL0001", ("GH1",)), _rule("PUL0002", ("GH1",)))
    matches, _ = screen.evaluate_pul_signatures(
        [screen.Hit("c1", "c1|gene1", "GH1", 1, 1e-20)], rules
    )
    assert [match.pul_id for match in matches] == ["PUL0001", "PUL0002"]


def test_clean_bulk_output_artifacts_preserves_unrelated_files(tmp_path: Path):
    output = tmp_path / "out"
    output.mkdir()
    (output / "resolved_dbcan_hmms").mkdir()
    (output / "resolved_dbcan_hmms" / "GH1.hmm").write_text("internal")
    (output / "hits_GH1.domtblout").write_text("internal")
    (output / "pul_matches.tsv").write_text("old")
    (output / "user-notes.txt").write_text("keep")

    screen._clean_bulk_output_artifacts(output)

    assert not (output / "resolved_dbcan_hmms").exists()
    assert not (output / "hits_GH1.domtblout").exists()
    assert not (output / "pul_matches.tsv").exists()
    assert (output / "user-notes.txt").read_text() == "keep"


def test_write_pul_outputs_has_compressed_tables(tmp_path: Path):
    matches = [
        screen.PULMatch("contig_A", "PUL0001", "starch", ("GH1", "GH2"), ("GH1", "GH2")),
    ]
    hits = [
        screen.Hit("contig_A", "protein_1", "GH1", 80, 1e-20, hmm_coverage=0.8),
        screen.Hit("contig_A", "protein_2", "GH2", 70, 1e-21, hmm_coverage=0.7),
    ]

    screen.write_pul_outputs(matches, hits, ["contig_A"], tmp_path)

    with gzip.open(tmp_path / "pul_matches.tsv.gz", "rt") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["reference_substrate"] == "starch"
    assert rows[0]["required_families"] == "GH1;GH2"
    assert (tmp_path / "pul_family_support.tsv.gz").exists()
    assert (tmp_path / "pul_summary.tsv").exists()
    assert (tmp_path / "substrate_summary.tsv").exists()


def test_resolve_all_puls_deduplicates_family_union(monkeypatch, tmp_path: Path):
    rules = {
        "PUL0001": _rule("PUL0001", ("GH1", "GH2")),
        "PUL0002": _rule("PUL0002", ("GH2", "GH3")),
    }
    calls = []
    monkeypatch.setattr(screen, "ensure_dbcan_database", lambda **_: {})
    monkeypatch.setattr(screen, "get_dbcan_pul_rules", lambda: rules)

    def extract(requested_ids, output_dir):
        calls.append(list(requested_ids))
        return [], []

    monkeypatch.setattr(screen, "extract_dbcan_models", extract)
    _, query, resolved = screen._resolve_all_puls_inputs(tmp_path / "out")
    assert calls == [["GH1", "GH2", "GH3"]]
    assert query.all_puls
    assert len(resolved) == 2


def test_all_puls_help_and_query_validation(tmp_path: Path):
    contigs = tmp_path / "contigs.fa"
    contigs.write_text(">c1\nATG\n")
    help_result = runner.invoke(app, ["screen", "--help"])
    assert help_result.exit_code == 0
    assert "--all-puls" in help_result.stdout

    no_query = runner.invoke(app, ["screen", "-i", str(contigs)])
    assert no_query.exit_code == 1
    assert "No HMM files specified" in no_query.stderr

    mixed = runner.invoke(app, ["screen", "-i", str(contigs), "--all-puls", "GH1"])
    assert mixed.exit_code == 1
    assert "cannot be combined" in mixed.stderr


def test_all_puls_works_without_positional_queries(monkeypatch, tmp_path: Path):
    contigs = tmp_path / "contigs.fa"
    contigs.write_text(">c1\nATG\n")
    captured = {}

    def fake_screen(config):
        captured["all_puls"] = config.all_puls
        captured["hmms"] = config.hmms

    monkeypatch.setattr("phu.cli._screen", fake_screen)

    result = runner.invoke(app, ["screen", "-i", str(contigs), "--all-puls"])

    assert result.exit_code == 0
    assert captured == {"all_puls": True, "hmms": []}