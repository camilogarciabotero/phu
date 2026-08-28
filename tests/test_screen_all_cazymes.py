import csv
import gzip
from pathlib import Path

from typer.testing import CliRunner

from phu import screen
from phu.cli import app

runner = CliRunner()


def test_all_cazyme_inventory_selects_canonical_profiles(monkeypatch, tmp_path: Path):
    inventory = (("AA1", "CBM89", "GH3", "GH5_2", "GT4", "PL6"), ("SLH", "cohesin"))
    calls = []
    monkeypatch.setattr(screen, "ensure_dbcan_database", lambda **_: {})
    monkeypatch.setattr(screen, "get_dbcan_model_inventory", lambda: inventory)

    def extract(requested_ids, output_dir):
        calls.append(list(requested_ids))
        return [tmp_path / f"{family}.hmm" for family in requested_ids], []

    monkeypatch.setattr(screen, "extract_dbcan_models", extract)
    paths, query = screen._resolve_all_cazymes_inputs(tmp_path / "out")

    assert [path.stem for path in paths] == list(inventory[0])
    assert calls == [list(inventory[0])]
    assert query.normalized_models == inventory[0]
    assert query.excluded_ancillary_profiles == inventory[1]
    assert query.total_hmm_profiles == 8


def test_write_cazyme_outputs_preserves_subfamilies_and_order(tmp_path: Path):
    hits = [
        screen.Hit("contig_B", "b|gene1", "GT4", 30, 1e-20, hmm_coverage=0.8),
        screen.Hit("contig_A", "a|gene1", "GH5_2", 40, 1e-25, hmm_coverage=0.9),
        screen.Hit("contig_B", "b|gene2", "CBM89", 50, 1e-30, hmm_coverage=0.7),
    ]
    output = tmp_path / "cazyme_matches.tsv.gz"

    screen.write_cazyme_outputs(hits, hits, ["contig_A", "contig_B"], tmp_path)

    with gzip.open(output, "rt") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [(row["contig_id"], row["cazyme_family"]) for row in rows] == [
        ("contig_A", "GH5_2"),
        ("contig_B", "CBM89"),
        ("contig_B", "GT4"),
    ]


def test_write_cazyme_outputs_aggregates_evidence_and_classes(tmp_path: Path):
    hits = [
        screen.Hit("contig_1", "protein_1", "GH128", 80, 1e-30, domain_i_evalue=1e-31, hmm_coverage=0.8),
        screen.Hit("contig_1", "protein_1", "GH128", 70, 1e-20, domain_i_evalue=1e-21, hmm_coverage=0.7),
        screen.Hit("contig_1", "protein_2", "GH128", 60, 1e-18, domain_i_evalue=1e-19, hmm_coverage=0.6),
        screen.Hit("contig_1", "protein_3", "CBM89", 55, 1e-17, domain_i_evalue=1e-18, hmm_coverage=0.7),
        screen.Hit("contig_2", "protein_4", "GH128", 50, 1e-16, domain_i_evalue=1e-17, hmm_coverage=0.5),
        screen.Hit("contig_2", "protein_5", "GT4", 45, 1e-16, domain_i_evalue=1e-17, hmm_coverage=0.5),
    ]

    counts = screen.write_cazyme_outputs(hits, hits, ["contig_1", "contig_2"], tmp_path)

    with gzip.open(tmp_path / "evidence/cazyme_hits.tsv.gz", "rt") as handle:
        evidence = list(csv.DictReader(handle, delimiter="\t"))
    with gzip.open(tmp_path / "cazyme_matches.tsv.gz", "rt") as handle:
        matches = list(csv.DictReader(handle, delimiter="\t"))
    with (tmp_path / "cazyme_class_summary.tsv").open() as handle:
        classes = list(csv.DictReader(handle, delimiter="\t"))

    assert len(evidence) == 6
    gh128 = next(row for row in matches if row["cazyme_family"] == "GH128")
    assert gh128["contig_id"] == "contig_1"
    assert gh128["protein_count"] == "2"
    assert gh128["hit_count"] == "3"
    assert gh128["best_protein_id"] == "protein_1"
    assert gh128["best_bitscore"] == "80.000000"
    with (tmp_path / "cazyme_summary.tsv").open() as handle:
        families = list(csv.DictReader(handle, delimiter="\t"))
    gh128_summary = next(row for row in families if row["cazyme_family"] == "GH128")
    assert gh128_summary["contig_count"] == "2"
    assert gh128_summary["protein_count"] == "3"
    assert gh128_summary["hit_count"] == "4"
    gh_class = next(row for row in classes if row["cazyme_class"] == "GH")
    assert gh_class["contig_count"] == "2"
    assert gh_class["protein_count"] == "3"
    assert counts == {
        "qualifying_hits": 6,
        "matched_proteins": 5,
        "matched_families": 3,
        "matched_classes": 3,
    }


def test_write_cazyme_outputs_keeps_subfamilies_distinct_and_empty_headers(tmp_path: Path):
    screen.write_cazyme_outputs([], [], [], tmp_path)

    assert (tmp_path / "cazyme_summary.tsv").read_text().splitlines() == [
        "cazyme_family\tcazyme_class\tcontig_count\tprotein_count\thit_count\tbest_bitscore\tbest_evalue"
    ]
    with gzip.open(tmp_path / "cazyme_matches.tsv.gz", "rt") as handle:
        assert handle.readline().strip() == (
            "contig_id\tcazyme_family\tcazyme_class\tprotein_count\thit_count\t"
            "protein_ids\tbest_protein_id\tbest_bitscore\tbest_evalue\t"
            "best_domain_i_evalue\tbest_hmm_coverage"
        )


def test_all_cazymes_cli_contract(tmp_path: Path, monkeypatch):
    contigs = tmp_path / "contigs.fa"
    contigs.write_text(">c1\nATG\n")
    help_result = runner.invoke(app, ["screen", "--help"])
    assert help_result.exit_code == 0
    assert "--all-cazymes" in help_result.stdout
    assert "All canonical CAZymes" in help_result.stdout
    captured = {}

    monkeypatch.setattr(
        "phu.cli._screen",
        lambda config: captured.update(
            {"all_cazymes": config.all_cazymes, "hmms": config.hmms}
        ),
    )
    result = runner.invoke(app, ["screen", "-i", str(contigs), "--all-cazymes"])
    assert result.exit_code == 0
    assert captured == {"all_cazymes": True, "hmms": []}

    for args, message in (
        (["--all-cazymes", "GH128"], "positional"),
        (["--all-cazymes", "PUL0621"], "positional"),
        (["--all-cazymes", "--all-puls"], "--all-puls"),
            (["--all-cazymes", "--combine-mode", "all"], "not applicable"),
    ):
        result = runner.invoke(app, ["screen", "-i", str(contigs), *args])
        assert result.exit_code == 1
        assert message in result.stderr


def test_normal_screen_without_query_still_fails(tmp_path: Path):
    contigs = tmp_path / "contigs.fa"
    contigs.write_text(">c1\nATG\n")
    result = runner.invoke(app, ["screen", "-i", str(contigs)])
    assert result.exit_code == 1
    assert "No HMM files specified" in result.stderr