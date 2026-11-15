"""`phu phrogs`: aliases, forgiving category matching and resolution to ids."""

from __future__ import annotations

import re

import pytest
from typer.testing import CliRunner

import phu.commands.phrogs as phrogs_cmd
from phu.cli import app
from phu.db import phrog_lookup as lookup
from phu.db.phrog import PhrogAnnotation

runner = CliRunner()

CATEGORIES = [
    "tail",
    "head and packaging",
    "lysis",
    "connector",
    "DNA, RNA and nucleotide metabolism",
    "moron, auxiliary metabolic gene and host takeover",
    "integration and excision",
    "unknown function",
]


def _phrogs() -> dict[str, PhrogAnnotation]:
    rows = [
        (1, "integrase", "integration and excision"),
        (2, "terminase large subunit", "head and packaging"),
        (3, "major head protein", "head and packaging"),
        (4, "Major capsid protein", "head and packaging"),
        (5, "minor head protein", "head and packaging"),
        (6, "tail fiber protein", "tail"),
        (7, "tail fiber assembly", "tail"),
        (8, "holin", "lysis"),
        (9, "neck protein", "connector"),
        (10, "DNA polymerase", "DNA, RNA and nucleotide metabolism"),
        (
            11,
            "phosphate starvation-inducible protein",
            "moron, auxiliary metabolic gene and host takeover",
        ),
        (12, "", "unknown function"),
    ]
    return {
        f"phrog_{n}": PhrogAnnotation(f"phrog_{n}", annot, cat)
        for n, annot, cat in rows
    }


@pytest.fixture
def table() -> lookup.AliasTable:
    return lookup.load_aliases()


def _ids(selected):
    return [r.phrog.phrog_id for r in selected]


# ---------- alias file ----------


def test_alias_file_is_well_formed(table):
    assert table.aliases, "no aliases"
    for alias in table.aliases.values():
        assert alias.name == alias.name.lower()
        assert alias.annotations and all(a.strip() for a in alias.annotations)
        assert len(set(alias.annotations)) == len(alias.annotations)
    assert set(table.categories.values()) <= set(CATEGORIES)


def test_check_aliases_reports_unmatched_strings(table):
    problems = lookup.check_aliases(table, _phrogs())
    # Only the strings missing from the small fake table are reported.
    assert problems["mcp"] == ["major head and protease protein"]
    assert problems["portal"] == ["portal protein"]
    assert "terl" not in problems


# ---------- categories ----------


@pytest.mark.parametrize(
    ("query", "expected"),
    [
        ("tail", ["tail"]),
        ("TAIL", ["tail"]),
        ("tails", ["tail"]),
        ("lysys", ["lysis"]),
        ("head", ["head and packaging"]),
        ("pack", ["head and packaging"]),
        ("dna", ["DNA, RNA and nucleotide metabolism"]),
        ("dna rna", ["DNA, RNA and nucleotide metabolism"]),
        ("moron", ["moron, auxiliary metabolic gene and host takeover"]),
        ("host takeover", ["moron, auxiliary metabolic gene and host takeover"]),
        ("connectr", ["connector"]),
        ("xyz", []),
    ],
)
def test_match_category(query, expected):
    assert lookup.match_category(query, CATEGORIES) == expected


def test_match_category_shortcut_and_ambiguity(table):
    assert lookup.match_category("amg", CATEGORIES, table.categories) == [
        "moron, auxiliary metabolic gene and host takeover"
    ]
    assert len(lookup.match_category("a", CATEGORIES)) > 1


# ---------- resolve ----------


def test_resolve_alias_uses_exact_annotations(table):
    selected, notes = lookup.resolve(["alias:mcp"], _phrogs(), table)
    # Case-insensitive exact match; "minor head protein" is not a capsid.
    assert _ids(selected) == ["phrog_3", "phrog_4"]
    assert notes == []


def test_resolve_category_notes_how_it_was_read(table):
    selected, notes = lookup.resolve(["category:tails"], _phrogs(), table)
    assert _ids(selected) == ["phrog_6", "phrog_7"]
    assert notes == ["category:tails read as 'tail'"]


def test_resolve_mixes_queries_without_repeats(table):
    selected, _ = lookup.resolve(
        ["phrog_6", "category:tail", "annot:Holin", "alias:terl"], _phrogs(), table
    )
    assert _ids(selected) == ["phrog_6", "phrog_7", "phrog_8", "phrog_2"]
    assert selected[0].query == "phrog_6" and selected[1].query == "category:tail"


@pytest.mark.parametrize(
    ("query", "message"),
    [
        ("foo", "Unrecognised query"),
        ("alias:mpc", "did you mean"),
        ("alias:portal", "selects no PHROG"),
        ("category:a", "ambiguous"),
        ("category:xyz", "No PHROG category matches"),
        ("phrog_999", "has no model"),
        ("annot:", "Unrecognised query"),
    ],
)
def test_resolve_errors(table, query, message):
    with pytest.raises(lookup.QueryError, match=message):
        lookup.resolve([query], _phrogs(), table)


def _found(text, table=None, **kw):
    return [
        (h.phrog.phrog_id, h.how)
        for h in lookup.search(text, _phrogs(), aliases=table, **kw)
    ]


def test_search_matches_all_words_any_order():
    assert _found("fiber tail") == [("phrog_7", "exact"), ("phrog_6", "exact")]
    assert _found("fiber", category="lysis") == []


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("tail fibre", [("phrog_7", "variant"), ("phrog_6", "variant")]),  # spelling
        ("holins", [("phrog_8", "variant")]),  # plural
        ("tailfiber protein", [("phrog_6", "variant")]),  # written together
        ("intergrase", [("phrog_1", "typo")]),
        ("polymerse", [("phrog_10", "typo")]),
        ("xyzzy", []),
    ],
)
def test_search_forgives_spelling(text, expected):
    assert _found(text) == expected


def test_search_alias_and_synonyms(table):
    assert _found("mcp", table) == [("phrog_4", "alias:mcp"), ("phrog_3", "alias:mcp")]
    assert _found("TerL", table) == [("phrog_2", "alias:terl")]
    # "capsid" also means "head protein"; the literal match ranks first.
    assert _found("capsid", table) == [
        ("phrog_4", "exact"),
        ("phrog_3", "synonym:capsid"),
        ("phrog_5", "synonym:capsid"),
    ]


def test_search_typos_need_five_letters():
    # "holn" is too short to be read as a typo of "holin".
    assert _found("holn") == []


def test_synonyms_are_search_only(table):
    with pytest.raises(lookup.QueryError):
        lookup.resolve(["annot:capsid"], _phrogs(), table)


# ---------- CLI ----------


@pytest.fixture
def cli(monkeypatch):
    monkeypatch.setattr(phrogs_cmd, "_phrogs", _phrogs)


def _plain(text: str) -> str:
    return re.sub(r"\x1b\[[0-?]*[ -/]*[@-~]", "", text)


def test_root_help_lists_phrogs():
    assert "phrogs" in _plain(runner.invoke(app, ["--help"]).stdout)


def test_cli_resolve_ids_only(cli):
    result = runner.invoke(
        app, ["phrogs", "resolve", "alias:mcp", "category:lysis", "--ids-only"]
    )
    assert result.exit_code == 0, result.output
    assert result.stdout.strip() == "phrog_3 phrog_4 phrog_8"
    assert "3 PHROG(s) selected" in result.stderr


def test_cli_resolve_writes_table(cli, tmp_path):
    out = tmp_path / "sel.tsv"
    result = runner.invoke(app, ["phrogs", "resolve", "category:tails", "-o", str(out)])
    assert result.exit_code == 0
    lines = out.read_text().splitlines()
    assert lines[0] == "phrog_id\tannotation\tcategory\tquery"
    assert lines[1] == "phrog_6\ttail fiber protein\ttail\tcategory:tails"
    assert "read as 'tail'" in result.stderr


def test_cli_resolve_error_exits_1(cli):
    result = runner.invoke(app, ["phrogs", "resolve", "category:a"])
    assert result.exit_code == 1
    assert "ambiguous" in result.stderr


def test_cli_categories_and_search(cli):
    result = runner.invoke(app, ["phrogs", "categories"])
    assert result.exit_code == 0
    assert "head and packaging\t4" in result.stdout
    result = runner.invoke(app, ["phrogs", "search", "head", "-c", "head"])
    assert result.exit_code == 0
    assert "phrog_3\tmajor head protein\thead and packaging\texact" in result.stdout
    assert "phrog_5\tminor head protein" in result.stdout
    assert result.stdout.splitlines()[0].endswith("\tmatch")


def test_cli_aliases_check_fails_on_missing(cli):
    result = runner.invoke(app, ["phrogs", "aliases", "--check"])
    assert result.exit_code == 1
    assert "portal" in result.stdout
