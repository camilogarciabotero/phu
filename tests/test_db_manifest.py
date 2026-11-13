import copy
import hashlib
import json
import runpy
import subprocess
import sys
from pathlib import Path

import pyhmmer
import pytest

from phu.db.manifest import (
    MANIFEST_MEMBER,
    ManifestError,
    Step,
    load_manifest,
    loads_manifest,
    parse_step,
)

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "phu_db"
FIXTURE_MANIFEST = FIXTURE_DIR / "phu-db.json"
GENERATOR = FIXTURE_DIR / "make_fixtures.py"

# Contract section 1: the full-bundle member set.
CONTRACT_MEMBERS = frozenset(
    {
        "phu-db/phu-db.json",
        "phu-db/pfam/profiles.hmm",
        "phu-db/kofam/profiles.hmm",
        "phu-db/kofam/ko_list.tsv",
        "phu-db/dbcan/profiles.hmm",
        "phu-db/dbcan/family_substrate.tsv",
        "phu-db/dbcan/pul_rules.tsv",
        "phu-db/avg/v_scores.tsv",
        "phu-db/avg/positive_amg.tsv",
        "phu-db/avg/positive_apg.tsv",
        "phu-db/avg/positive_areg.tsv",
        "phu-db/avg/filters_amg.tsv",
        "phu-db/avg/filters_apg.tsv",
        "phu-db/avg/filters_areg.tsv",
        "phu-db/phrog/profiles.hmm",
        "phu-db/phrog/annotations.tsv",
    }
)

# Contract section 4: every step of the 0.2.0 bundle.
CONTRACT_STEPS = (
    "gunzip:Pfam-A.hmm.gz -> profiles.hmm",
    "extract:profiles.tar.gz -> profiles/",
    "concatenate:profiles/*.hmm -> profiles.hmm",
    "gunzip:ko_list.gz -> ko_list.tsv",
    "rename:dbCAN-HMMdb-V15.txt -> profiles.hmm",
    "rename:dbCAN-PUL-substrate-metadata.txt -> family_substrate.tsv",
    "pul-table:dbCAN-PUL.xlsx -> pul_rules.tsv",
    "rename:vscores.tsv -> v_scores.tsv",
    "rename:AMGs.tsv -> positive_amg.tsv",
    "rename:APGs.tsv -> positive_apg.tsv",
    "rename:AReGs.tsv -> positive_areg.tsv",
    "rename:AMG_filters.tsv -> filters_amg.tsv",
    "rename:APG_filters.tsv -> filters_apg.tsv",
    "rename:AReG_filters.tsv -> filters_areg.tsv",
    "extract:MSA_phrogs.tar.gz -> MSA_Phrogs_M50_FASTA/",
    "build-hmm:MSA_Phrogs_M50_FASTA/*.fma -> profiles.hmm",
    "rename:phrog_annot_v4.tsv -> annotations.tsv",
)

# First lines of the published 0.1.0 bundle (Zenodo 22928399). The three
# filters_* files carry the header 0.1.0 kept when it concatenated them.
REAL_HEADERS = {
    "kofam/ko_list.tsv": (
        "knum\tthreshold\tscore_type\tprofile_type\tF-measure\tnseq\tnseq_used\talen"
        "\tmlen\teff_nseq\tre/pos\tdefinition"
    ),
    "avg/v_scores.tsv": "id\tname\tV-score\tVL-score\tdb",
    "avg/positive_amg.tsv": (
        "id\tdb\tname\tmetabolic_paths\ttotal_paths\tpath_source\tmetabolic_ratio"
        "\tVL-score\tAMG_weight\tamg_level"
    ),
    "avg/positive_apg.tsv": "id\tV-score\tVL-score\tdb\tname",
    "avg/positive_areg.tsv": "id\tV-score\tVL-score\tdb\tname",
    "avg/filters_amg.tsv": (
        "id\tdb\tname\tfilter_essential\tfilter_glucan\tfilter_lipid\tfilter_methyl"
        "\tfilter_nucleotide"
    ),
    "avg/filters_apg.tsv": (
        "id\tdb\tname\tfilter_essential\tfilter_glucan\tfilter_lipid\tfilter_methyl"
        "\tfilter_nucleotide"
    ),
    "avg/filters_areg.tsv": (
        "id\tdb\tname\tfilter_essential\tfilter_glucan\tfilter_lipid\tfilter_methyl"
        "\tfilter_nucleotide"
    ),
}


def _raw() -> dict:
    return json.loads(FIXTURE_MANIFEST.read_text())


def _loads(raw: dict):
    return loads_manifest(json.dumps(raw))


def test_fixture_manifest_loads():
    manifest = load_manifest(FIXTURE_MANIFEST)

    assert manifest.bundle_version == "phu-db-0.2.0-fixture"
    assert manifest.layout_version == 1
    assert list(manifest.databases) == ["pfam", "kofam", "dbcan", "avg", "phrog"]
    assert manifest.expected_members() == CONTRACT_MEMBERS
    assert MANIFEST_MEMBER in manifest.expected_members()
    assert manifest.member_path("kofam", "ko_list.tsv") == "phu-db/kofam/ko_list.tsv"

    kofam = manifest.databases["kofam"]
    assert [source.name for source in kofam.sources] == ["profiles.tar.gz", "ko_list.gz"]
    assert kofam.steps == CONTRACT_STEPS[1:4]
    assert kofam.parsed_steps[0] == Step("extract", "profiles.tar.gz", "profiles")
    assert kofam.parsed_steps[1].is_glob
    assert [item.path for item in manifest.databases["dbcan"].prepared_files] == [
        "profiles.hmm",
        "family_substrate.tsv",
        "pul_rules.tsv",
    ]
    all_steps = tuple(step for entry in manifest.databases.values() for step in entry.steps)
    assert all_steps == CONTRACT_STEPS


def test_fixture_manifest_key_order_matches_bundle_schema():
    raw = _raw()

    assert list(raw) == ["bundle_version", "layout_version", "generator", "bundle", "zenodo", "databases"]
    for entry in raw["databases"].values():
        assert list(entry) == [
            "provenance",
            "redistribution",
            "upstream_version",
            "license",
            "citation",
            "layout_dir",
            "sources",
            "steps",
            "prepared_files",
            "validation",
            "build",
        ]
        assert all(list(source) == ["name", "url", "sha256", "bytes"] for source in entry["sources"])
        assert all(
            list(item) == ["path", "format", "bytes", "sha256", "count"]
            for item in entry["prepared_files"]
        )


def _set(*keys, value):
    def mutate(raw):
        target = raw
        for key in keys[:-1]:
            target = target[key]
        target[keys[-1]] = value

    return mutate


def _delete(*keys):
    def mutate(raw):
        target = raw
        for key in keys[:-1]:
            target = target[key]
        del target[keys[-1]]

    return mutate


def _add_database(raw):
    raw["databases"]["vscore"] = copy.deepcopy(raw["databases"]["avg"])
    raw["databases"]["vscore"]["layout_dir"] = "vscore"


def _duplicate_source(raw):
    sources = raw["databases"]["kofam"]["sources"]
    sources[1]["name"] = sources[0]["name"]


def _swap_kofam_steps(raw):
    steps = raw["databases"]["kofam"]["steps"]
    steps[0], steps[1] = steps[1], steps[0]


def _duplicate_prepared(raw):
    files = raw["databases"]["kofam"]["prepared_files"]
    files.append(copy.deepcopy(files[0]))


def _extend_kofam(step, path, fmt="tsv"):
    """Append a step and a prepared file it produces to kofam."""

    def mutate(raw):
        entry = raw["databases"]["kofam"]
        entry["steps"].append(step)
        prepared = copy.deepcopy(entry["prepared_files"][-1])
        prepared.update(path=path, format=fmt)
        entry["prepared_files"].append(prepared)

    return mutate


def _extend_pfam_with_split(raw):
    """Give pfam two gunzipped file targets and a glob that joins them."""
    entry = raw["databases"]["pfam"]
    entry["steps"] = [
        "gunzip:Pfam-A.hmm.gz -> a.hmm",
        "gunzip:Pfam-A.hmm.gz -> b.hmm",
        "concatenate:*.hmm -> profiles.hmm",
    ]


PFAM_STEP = ("databases", "pfam", "steps", 0)
PFAM_SOURCE = ("databases", "pfam", "sources", 0)
PFAM_PREPARED = ("databases", "pfam", "prepared_files", 0)


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        pytest.param(_set("layout_version", value=2), "layout_version", id="layout-version-2"),
        pytest.param(_set("layout_version", value="1"), "layout_version", id="layout-version-string"),
        pytest.param(_delete("databases", "avg"), "databases must be", id="database-missing"),
        pytest.param(_add_database, "databases must be", id="database-extra"),
        pytest.param(_set(*PFAM_SOURCE, "sha256", value="AB" * 32), "sha256", id="sha256-uppercase"),
        pytest.param(_set(*PFAM_SOURCE, "sha256", value="ab" * 31), "sha256", id="sha256-short"),
        pytest.param(_set(*PFAM_SOURCE, "sha256", value="g" * 64), "sha256", id="sha256-non-hex"),
        pytest.param(_set(*PFAM_PREPARED, "sha256", value="0" * 65), "sha256", id="sha256-prepared-long"),
        pytest.param(_delete(*PFAM_SOURCE, "name"), "missing 'name'", id="source-name-missing"),
        pytest.param(_delete("databases", "pfam", "steps"), "missing 'steps'", id="steps-missing"),
        pytest.param(
            _set("databases", "pfam", "steps", value="gunzip:Pfam-A.hmm.gz -> profiles.hmm"),
            "wrong type",
            id="steps-string",
        ),
        pytest.param(_duplicate_source, "duplicate names", id="source-name-duplicate"),
        pytest.param(_set(*PFAM_SOURCE, "name", value="sub/Pfam-A.hmm.gz"), "unsafe", id="source-name-nested"),
        pytest.param(_set(*PFAM_STEP, value="gunzip:Pfam-A.hmm.gz->profiles.hmm"), "malformed", id="step-grammar"),
        pytest.param(_set(*PFAM_STEP, value="unzip:Pfam-A.hmm.gz -> profiles.hmm"), "malformed", id="step-verb"),
        pytest.param(_set(*PFAM_STEP, value="gunzip:Pfam-B.hmm.gz -> profiles.hmm"), "neither", id="step-unknown-source"),
        pytest.param(_swap_kofam_steps, "neither", id="step-uses-later-target"),
        pytest.param(_set(*PFAM_STEP, value="gunzip:Pfam-A.hmm.gz -> ../profiles.hmm"), "unsafe", id="step-target-dotdot"),
        pytest.param(_set(*PFAM_STEP, value="gunzip:Pfam-A.hmm.gz -> /profiles.hmm"), "unsafe", id="step-target-absolute"),
        pytest.param(_set(*PFAM_STEP, value="gunzip:Pfam-A.hmm.gz -> a\\profiles.hmm"), "unsafe", id="step-target-backslash"),
        pytest.param(_set(*PFAM_PREPARED, "path", value="/profiles.hmm"), "unsafe", id="prepared-absolute"),
        pytest.param(_set(*PFAM_PREPARED, "path", value="../profiles.hmm"), "unsafe", id="prepared-dotdot"),
        pytest.param(_set(*PFAM_PREPARED, "path", value="pfam\\profiles.hmm"), "unsafe", id="prepared-backslash"),
        pytest.param(_set(*PFAM_PREPARED, "path", value="other.hmm"), "not the target", id="prepared-not-produced"),
        pytest.param(
            _set("databases", "kofam", "prepared_files", 0, "path", value="profiles"),
            "not the target",
            id="prepared-is-extract-directory",
        ),
        pytest.param(
            _extend_kofam("concatenate:*.txt -> all.txt", "all.txt"),
            "neither",
            id="step-glob-matches-no-earlier-file",
        ),
        pytest.param(
            _extend_kofam("rename:ko_list.gz -> profiles", "profiles"),
            "conflicts",
            id="step-file-target-is-earlier-directory",
        ),
        pytest.param(
            _extend_kofam("rename:ko_list.gz -> ko_list.tsv/nested.tsv", "ko_list.tsv/nested.tsv"),
            "conflicts",
            id="step-target-below-earlier-file",
        ),
        pytest.param(_set("databases", "pfam", "layout_dir", value="../pfam"), "unsafe", id="layout-dir-dotdot"),
        pytest.param(_set("databases", "pfam", "layout_dir", value="/pfam"), "unsafe", id="layout-dir-absolute"),
        pytest.param(_set("databases", "pfam", "layout_dir", value="pfam\\x"), "unsafe", id="layout-dir-backslash"),
        pytest.param(_set("databases", "pfam", "layout_dir", value="Pfam"), "layout_dir must be", id="layout-dir-mismatch"),
        pytest.param(_set(*PFAM_SOURCE, "bytes", value="12"), "wrong type", id="source-bytes-string"),
        pytest.param(_set(*PFAM_SOURCE, "bytes", value=None), "wrong type", id="source-bytes-null"),
        pytest.param(_set(*PFAM_SOURCE, "bytes", value=-1), "negative", id="source-bytes-negative"),
        pytest.param(_set(*PFAM_PREPARED, "bytes", value=1.5), "wrong type", id="prepared-bytes-float"),
        pytest.param(_set(*PFAM_PREPARED, "bytes", value=True), "wrong type", id="prepared-bytes-bool"),
        pytest.param(_set(*PFAM_PREPARED, "bytes", value=-5), "negative", id="prepared-bytes-negative"),
        pytest.param(_set(*PFAM_PREPARED, "count", value="2"), "wrong type", id="prepared-count-string"),
        pytest.param(_set(*PFAM_PREPARED, "count", value=-1), "negative", id="prepared-count-negative"),
        pytest.param(_set(*PFAM_PREPARED, "format", value="fasta"), "format", id="prepared-format-unknown"),
        pytest.param(_duplicate_prepared, "duplicate paths", id="prepared-duplicate"),
        pytest.param(_set("bundle_version", value="../phu-db"), "unsafe", id="bundle-version-unsafe"),
    ],
)
def test_load_manifest_rejects(mutate, message):
    raw = _raw()
    mutate(raw)

    with pytest.raises(ManifestError, match=message):
        _loads(raw)


def test_step_may_consume_an_earlier_literal_target():
    raw = _raw()
    _extend_kofam("rename:ko_list.tsv -> copy.tsv", "copy.tsv")(raw)

    entry = _loads(raw).databases["kofam"]

    assert entry.steps[-1] == "rename:ko_list.tsv -> copy.tsv"
    assert entry.prepared_files[-1].path == "copy.tsv"


def test_concatenate_glob_may_range_over_earlier_file_targets():
    raw = _raw()
    _extend_pfam_with_split(raw)

    entry = _loads(raw).databases["pfam"]

    assert entry.parsed_steps[-1] == Step("concatenate", "*.hmm", "profiles.hmm")
    assert entry.parsed_steps[-1].is_glob


def test_step_target_may_sit_inside_an_earlier_directory():
    raw = _raw()
    _extend_kofam("rename:ko_list.gz -> profiles/ko_list.gz", "profiles/ko_list.gz")(raw)

    assert _loads(raw).databases["kofam"].prepared_files[-1].path == "profiles/ko_list.gz"


def test_manifest_is_hashable():
    manifest = load_manifest(FIXTURE_MANIFEST)

    assert hash(manifest) == hash(load_manifest(FIXTURE_MANIFEST))


def test_load_manifest_rejects_duplicate_json_keys():
    text = FIXTURE_MANIFEST.read_text().replace(
        '"layout_version": 1,', '"layout_version": 1,\n  "layout_version": 1,', 1
    )

    with pytest.raises(ManifestError, match="duplicate key"):
        loads_manifest(text)


def test_load_manifest_rejects_invalid_json():
    with pytest.raises(ManifestError, match="not valid JSON"):
        loads_manifest(b'{"layout_version": 1')


def test_load_manifest_rejects_deeply_nested_json():
    with pytest.raises(ManifestError, match="not valid JSON"):
        loads_manifest(b"[" * 200_000)


def test_unknown_keys_are_ignored():
    raw = _raw()
    raw["future_top_level"] = {"anything": [1, 2]}
    raw["databases"]["avg"]["validation"] = {"passed": False, "reason": "ignored"}
    raw["databases"]["avg"]["build"] = "not an object"
    raw["databases"]["pfam"]["future_db_key"] = 3
    raw["databases"]["pfam"]["sources"][0]["mirror"] = "https://example.invalid"
    raw["databases"]["pfam"]["prepared_files"][0]["note"] = None

    assert _loads(raw) == load_manifest(FIXTURE_MANIFEST)


@pytest.mark.parametrize("text", CONTRACT_STEPS)
def test_parse_step_accepts_contract_steps(text):
    step = parse_step(text)

    verb, rest = text.split(":", 1)
    source, target = rest.split(" -> ")
    assert step == Step(verb, source, target.rstrip("/"))
    assert step.is_directory == (verb == "extract")


@pytest.mark.parametrize(
    "text",
    [
        pytest.param("gunzip:Pfam-A.hmm.gz profiles.hmm", id="missing-arrow"),
        pytest.param("gunzip:Pfam-A.hmm.gz->profiles.hmm", id="arrow-without-spaces"),
        pytest.param("unzip:Pfam-A.hmm.gz -> profiles.hmm", id="unknown-verb"),
        pytest.param("Gunzip:Pfam-A.hmm.gz -> profiles.hmm", id="verb-case"),
        pytest.param("rename:AMGs.tsv -> positive/", id="trailing-slash-rename"),
        pytest.param("gunzip:ko_list.gz -> ko_list/", id="trailing-slash-gunzip"),
        pytest.param("concatenate:profiles/*.hmm -> profiles/", id="trailing-slash-concatenate"),
        pytest.param("gunzip:Pfam-A.hmm.gz -> ../profiles.hmm", id="dotdot-target"),
        pytest.param("gunzip:Pfam-A.hmm.gz -> pfam/../../profiles.hmm", id="nested-dotdot-target"),
        pytest.param("extract:profiles.tar.gz -> ../profiles/", id="dotdot-extract-target"),
        pytest.param("gunzip:Pfam-A.hmm.gz -> /profiles.hmm", id="absolute-target"),
        pytest.param("gunzip:/Pfam-A.hmm.gz -> profiles.hmm", id="absolute-source"),
        pytest.param("concatenate:../*.hmm -> profiles.hmm", id="dotdot-source"),
        pytest.param("gunzip:Pfam-A.hmm.gz -> a\\b.hmm", id="backslash-target"),
        pytest.param("gunzip:Pfam-A.hmm.gz -> a//b.hmm", id="empty-segment-target"),
        pytest.param("extract:profiles.tar.gz -> profiles//", id="double-trailing-slash"),
        pytest.param(":Pfam-A.hmm.gz -> profiles.hmm", id="empty-verb"),
        pytest.param("gunzip: -> profiles.hmm", id="empty-source"),
        pytest.param("gunzip:Pfam-A.hmm.gz -> ", id="empty-target"),
        pytest.param("gunzip:Pfam-A.hmm.gz -> profiles.hmm\n", id="trailing-newline"),
        pytest.param(" gunzip:Pfam-A.hmm.gz -> profiles.hmm", id="leading-space"),
        pytest.param("gunzip:Pfam-A.hmm.gz -> profiles.hmm -> x", id="two-arrows"),
        pytest.param("", id="empty"),
    ],
)
def test_parse_step_rejects_malformed(text):
    with pytest.raises(ManifestError):
        parse_step(text)


# What make_fixtures.py owns. Later packages keep their own fixtures beside it.
GENERATED = ("upstream", "prepared", "phu-db.json")


def _generated_tree(root: Path) -> dict[str, bytes]:
    files = {}
    for name in GENERATED:
        top = root / name
        for path in sorted(top.rglob("*")) if top.is_dir() else [top]:
            if path.is_file():
                files[path.relative_to(root).as_posix()] = path.read_bytes()
    return files


def test_fixture_generator_is_byte_reproducible(tmp_path: Path):
    generator = runpy.run_path(str(GENERATOR), run_name="make_fixtures")
    if pyhmmer.__version__ != generator["PYHMMER_VERSION"]:
        pytest.skip(
            f"fixture HMMs were written by pyhmmer {generator['PYHMMER_VERSION']}, "
            f"this is {pyhmmer.__version__}"
        )
    subprocess.run([sys.executable, str(GENERATOR), "--out", str(tmp_path)], check=True)

    generated = _generated_tree(tmp_path)
    committed = _generated_tree(FIXTURE_DIR)
    assert sorted(generated) == sorted(committed)
    for name, data in generated.items():
        assert data == committed[name], f"{name} differs from the committed fixture"


def test_fixture_files_match_manifest():
    raw = _raw()

    for name, entry in raw["databases"].items():
        upstream = FIXTURE_DIR / "upstream" / name
        assert sorted(path.name for path in upstream.iterdir()) == sorted(
            source["name"] for source in entry["sources"]
        )
        for source in entry["sources"]:
            data = (upstream / source["name"]).read_bytes()
            assert source["bytes"] == len(data)
            assert source["sha256"] == hashlib.sha256(data).hexdigest()

        prepared = FIXTURE_DIR / "prepared" / name
        assert sorted(path.name for path in prepared.iterdir()) == sorted(
            item["path"] for item in entry["prepared_files"]
        )
        for item in entry["prepared_files"]:
            data = (prepared / item["path"]).read_bytes()
            lines = data.split(b"\n")
            if item["format"] == "hmm":
                count = sum(1 for line in lines if line.startswith(b"NAME "))
            else:
                count = sum(1 for line in lines if line.strip())
            assert item["bytes"] == len(data), item["path"]
            assert item["sha256"] == hashlib.sha256(data).hexdigest(), item["path"]
            assert item["count"] == count, item["path"]


def test_fixture_counts_follow_contract_semantics():
    manifest = load_manifest(FIXTURE_MANIFEST)
    counts = {
        f"{name}/{item.path}": item.count
        for name, entry in manifest.databases.items()
        for item in entry.prepared_files
    }

    assert counts["pfam/profiles.hmm"] == 2
    assert counts["kofam/profiles.hmm"] == 3
    assert counts["dbcan/profiles.hmm"] == 2
    # Header included where there is one; family_substrate.tsv has none.
    assert counts["kofam/ko_list.tsv"] == 1 + 3
    assert counts["dbcan/family_substrate.tsv"] == 3
    assert counts["dbcan/pul_rules.tsv"] == 1 + 3
    for path in REAL_HEADERS:
        if path.startswith("avg/"):
            rows = {"avg/v_scores.tsv": 8, "avg/positive_amg.tsv": 5}.get(path, 6)
            assert counts[path] == 1 + rows
    assert counts["phrog/profiles.hmm"] == 3
    # One annotation row has no model, like the singleton PHROGs upstream.
    assert counts["phrog/annotations.tsv"] == 1 + 4


@pytest.mark.parametrize("path", sorted(REAL_HEADERS))
def test_fixture_headers_match_real_bundle(path):
    first_line = (FIXTURE_DIR / "prepared" / path).read_text().split("\n", 1)[0]

    assert first_line == REAL_HEADERS[path]


def test_fixture_family_substrate_is_headerless_like_real_bundle():
    first_line = (FIXTURE_DIR / "prepared" / "dbcan" / "family_substrate.tsv").read_text().split("\n", 1)[0]

    assert first_line.split("\t")[0] == "PUL0001"


def test_fixture_ko_list_has_a_dash_threshold_row():
    rows = (FIXTURE_DIR / "prepared" / "kofam" / "ko_list.tsv").read_text().splitlines()[1:]

    assert any(row.split("\t")[1] == "-" for row in rows)


def test_fixture_pul_rules_follow_pul_table_semantics():
    text = (FIXTURE_DIR / "prepared" / "dbcan" / "pul_rules.tsv").read_bytes()

    assert text == (
        b"pul_id\tsubstrate\tcazymes_predicted_dbcan\n"
        b"PUL0001\txylan\tGH5,CBM48\n"
        b"PUL0002\tstarch glycogen\tGH5|CBM48\n"
        b"PUL0003\tbeta-glucan\tGH5, GH999\n"
    )


def test_fixture_pul_rules_agree_with_phu_rule_reader():
    from phu.db.dbcan import parse_pul_rules

    rules = parse_pul_rules(FIXTURE_DIR / "prepared" / "dbcan" / "pul_rules.tsv", {"GH5", "CBM48"})
    rows = (FIXTURE_DIR / "prepared" / "dbcan" / "pul_rules.tsv").read_text().splitlines()[1:]

    assert [row.split("\t") for row in rows] == [
        [rule.pul_id, rule.substrate, rule.raw_rule] for rule in rules.values()
    ]
    assert rules["PUL0003"].unresolved_tokens == ("GH999",)
