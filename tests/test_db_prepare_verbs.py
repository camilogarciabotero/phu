import gzip
import hashlib
import io
import tarfile
import zipfile
from pathlib import Path

import pytest

from phu.db_fetch import DatabaseError
from phu.db_manifest import Step
from phu.db_prepare_verbs import (
    PUL_TABLE_HEADER,
    PrepareError,
    UnsafeArchive,
    concatenate,
    extract,
    gunzip,
    pul_table,
    rename,
    run_steps,
)

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "phu_db"


def _files_below(path: Path) -> set[str]:
    return {item.relative_to(path).as_posix() for item in path.rglob("*")}


# The fixture upstreams prepare into exactly the fixture manifest's files.


@pytest.mark.parametrize("database", ["pfam", "kofam", "dbcan", "avg"])
def test_run_steps_reproduces_fixture_prepared_files(tmp_path, fixture_manifest, database):
    entry = fixture_manifest.databases[database]
    workdir = tmp_path / "work"

    produced = run_steps(
        entry.steps, workdir=workdir, downloads=FIXTURE_DIR / "upstream" / database
    )

    assert produced == [workdir / step.target for step in entry.parsed_steps]
    for item in entry.prepared_files:
        data = (workdir / item.path).read_bytes()
        assert (len(data), hashlib.sha256(data).hexdigest()) == (item.bytes, item.sha256), item.path
        assert data == (FIXTURE_DIR / "prepared" / database / item.path).read_bytes()


def test_run_steps_leaves_kofam_profiles_beside_the_concatenation(tmp_path, fixture_manifest):
    workdir = tmp_path / "work"
    run_steps(
        fixture_manifest.databases["kofam"].steps,
        workdir=workdir,
        downloads=FIXTURE_DIR / "upstream" / "kofam",
    )
    assert _files_below(workdir) == {
        "ko_list.tsv",
        "profiles",
        "profiles.hmm",
        "profiles/K00001.hmm",
        "profiles/K00002.hmm",
        "profiles/K00003.hmm",
    }


# run_steps source resolution.


def test_earlier_target_shadows_a_download_of_the_same_name(tmp_path):
    downloads = tmp_path / "downloads"
    downloads.mkdir()
    (downloads / "a.txt").write_bytes(b"from a\n")
    (downloads / "b.txt").write_bytes(b"downloaded b\n")
    workdir = tmp_path / "work"

    run_steps(["rename:a.txt -> b.txt", "rename:b.txt -> c.txt"], workdir=workdir, downloads=downloads)

    assert (workdir / "c.txt").read_bytes() == b"from a\n"


def test_concatenate_prefers_workdir_matches_over_downloads(tmp_path):
    downloads = tmp_path / "downloads"
    (downloads / "parts").mkdir(parents=True)
    (downloads / "parts" / "x.hmm").write_bytes(b"download\n")
    workdir = tmp_path / "work"
    (workdir / "parts").mkdir(parents=True)
    (workdir / "parts" / "y.hmm").write_bytes(b"workdir\n")

    concatenate("parts/*.hmm", workdir, downloads, workdir / "all.hmm")

    assert (workdir / "all.hmm").read_bytes() == b"workdir\n"


def test_concatenate_falls_back_to_downloads(tmp_path):
    downloads = tmp_path / "downloads"
    downloads.mkdir()
    (downloads / "x.hmm").write_bytes(b"download\n")
    workdir = tmp_path / "work"
    workdir.mkdir()

    concatenate("*.hmm", workdir, downloads, workdir / "all.hmm")

    assert (workdir / "all.hmm").read_bytes() == b"download\n"


def test_run_steps_missing_source_raises(tmp_path):
    with pytest.raises(PrepareError, match="neither an earlier target nor a download"):
        run_steps(["rename:absent.txt -> b.txt"], workdir=tmp_path / "w", downloads=tmp_path / "d")


@pytest.mark.parametrize(
    "step", ["rename:a.txt->b.txt", "copy:a.txt -> b.txt", "rename:a.txt -> ../b.txt", "rename:a -> b/"]
)
def test_run_steps_malformed_step_raises(tmp_path, step):
    with pytest.raises(PrepareError):
        run_steps([step], workdir=tmp_path / "w", downloads=tmp_path)


def test_run_steps_accepts_parsed_steps(tmp_path, fixture_manifest):
    entry = fixture_manifest.databases["kofam"]
    workdir = tmp_path / "work"
    run_steps(entry.parsed_steps, workdir=workdir, downloads=FIXTURE_DIR / "upstream" / "kofam")
    for item in entry.prepared_files:
        assert hashlib.sha256((workdir / item.path).read_bytes()).hexdigest() == item.sha256


@pytest.mark.parametrize(
    "step",
    [
        Step("copy", "a.txt", "b.txt"),
        Step("rename", "../secret.txt", "b.txt"),
        Step("rename", "a.txt", "../escaped.txt"),
        Step("gunzip", "/abs.gz", "b.txt"),
        Step("rename", "a.txt", "b/"),
    ],
    ids=["unknown-verb", "climbing-source", "climbing-target", "absolute-source", "slash-target"],
)
def test_run_steps_checks_hand_built_steps(tmp_path, step):
    downloads = tmp_path / "d"
    downloads.mkdir()
    (downloads / "a.txt").write_bytes(b"a\n")
    (tmp_path / "secret.txt").write_bytes(b"secret\n")
    with pytest.raises(PrepareError):
        run_steps([step], workdir=tmp_path / "w" / "work", downloads=downloads)
    assert _files_below(tmp_path) == {"d", "d/a.txt", "secret.txt", "w", "w/work"}


def test_prepare_error_is_a_database_error():
    assert issubclass(PrepareError, DatabaseError)


# gunzip


def test_gunzip_decompresses(tmp_path):
    src = tmp_path / "x.gz"
    src.write_bytes(gzip.compress(b"HMMER3/f\n" * 1000))
    gunzip(src, tmp_path / "out" / "x")
    assert (tmp_path / "out" / "x").read_bytes() == b"HMMER3/f\n" * 1000


@pytest.mark.parametrize(
    "data",
    [
        b"not gzip at all",
        gzip.compress(b"HMMER3/f\n" * 1000)[:-30],
        gzip.compress(b"HMMER3/f\n" * 1000)[:12] + b"\xff" * 40,
    ],
    ids=["not-gzip", "truncated", "corrupt-deflate"],
)
def test_gunzip_bad_input_raises_and_leaves_nothing(tmp_path, data):
    src = tmp_path / "in" / "x.gz"
    src.parent.mkdir()
    src.write_bytes(data)
    out = tmp_path / "out"
    out.mkdir()

    with pytest.raises(PrepareError, match="cannot gunzip x.gz"):
        gunzip(src, out / "x")
    assert _files_below(out) == set()


def test_failed_verb_keeps_an_existing_target(tmp_path):
    src = tmp_path / "x.gz"
    src.write_bytes(b"not gzip")
    dst = tmp_path / "x"
    dst.write_bytes(b"previous")

    with pytest.raises(PrepareError):
        gunzip(src, dst)
    assert dst.read_bytes() == b"previous"
    assert _files_below(tmp_path) == {"x.gz", "x"}


# extract


def _info(name, type=tarfile.REGTYPE, linkname="", size=0):
    info = tarfile.TarInfo(name)
    info.type = type
    info.linkname = linkname
    info.size = size
    return info


def _tar_gz(path: Path, entries) -> Path:
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as tar:
        for info, data in entries:
            if data is None:
                tar.addfile(info)
            else:
                info.size = len(data)
                tar.addfile(info, io.BytesIO(data))
    path.write_bytes(gzip.compress(buffer.getvalue()))
    return path


def test_extract_unpacks_into_workdir_root(tmp_path):
    archive = _tar_gz(
        tmp_path / "a.tar.gz",
        [(_info("profiles", tarfile.DIRTYPE), None), (_info("profiles/K1.hmm"), b"HMM\n")],
    )
    workdir = tmp_path / "work"
    workdir.mkdir()

    assert extract(archive, workdir, Path("profiles")) == workdir / "profiles"
    assert _files_below(workdir) == {"profiles", "profiles/K1.hmm"}
    assert (workdir / "profiles" / "K1.hmm").read_bytes() == b"HMM\n"


def test_extract_without_directory_members_still_produces_target(tmp_path):
    archive = _tar_gz(tmp_path / "a.tar.gz", [(_info("profiles/K1.hmm"), b"HMM\n")])
    workdir = tmp_path / "work"
    workdir.mkdir()
    extract(archive, workdir, Path("profiles"))
    assert _files_below(workdir) == {"profiles", "profiles/K1.hmm"}


@pytest.mark.parametrize(
    "entry",
    [
        (_info("profiles/K2.hmm", tarfile.SYMTYPE, "/etc/passwd"), None),
        (_info("profiles/K2.hmm", tarfile.LNKTYPE, "profiles/K1.hmm"), None),
        (_info("../K2.hmm"), b"x"),
        (_info("profiles/../../K2.hmm"), b"x"),
        (_info("/tmp/K2.hmm"), b"x"),
    ],
    ids=["symlink", "hardlink", "leading-dotdot", "inner-dotdot", "absolute"],
)
def test_extract_rejects_unsafe_members_and_writes_nothing(tmp_path, entry):
    archive = _tar_gz(tmp_path / "a.tar.gz", [(_info("profiles/K1.hmm"), b"HMM\n"), entry])
    workdir = tmp_path / "work"
    workdir.mkdir()

    with pytest.raises(PrepareError) as raised:
        extract(archive, workdir, Path("profiles"))
    assert isinstance(raised.value.__cause__, UnsafeArchive)
    assert _files_below(workdir) == set()


def test_extract_raises_when_target_dir_is_absent(tmp_path):
    archive = _tar_gz(tmp_path / "a.tar.gz", [(_info("models/K1.hmm"), b"HMM\n")])
    workdir = tmp_path / "work"
    workdir.mkdir()

    with pytest.raises(PrepareError, match="did not produce profiles/"):
        extract(archive, workdir, Path("profiles"))
    assert _files_below(workdir) == set()


def test_extract_corrupt_archive_raises(tmp_path):
    archive = tmp_path / "a.tar.gz"
    archive.write_bytes(b"\x1f\x8b not really")
    workdir = tmp_path / "work"
    workdir.mkdir()

    with pytest.raises(PrepareError, match="cannot extract a.tar.gz"):
        extract(archive, workdir, Path("profiles"))
    assert _files_below(workdir) == set()


def test_extract_into_missing_workdir_raises(tmp_path):
    archive = _tar_gz(tmp_path / "a.tar.gz", [(_info("profiles/K1.hmm"), b"x")])
    with pytest.raises(PrepareError, match="cannot extract"):
        extract(archive, tmp_path / "absent", Path("profiles"))


def test_extract_target_must_stay_in_workdir(tmp_path):
    archive = _tar_gz(tmp_path / "a.tar.gz", [(_info("profiles/K1.hmm"), b"HMM\n")])
    workdir = tmp_path / "work"
    workdir.mkdir()
    for target in (Path("../profiles"), tmp_path / "elsewhere", Path(".")):
        with pytest.raises(PrepareError, match="not inside"):
            extract(archive, workdir, target)


# concatenate


def _write(base: Path, files: dict[str, bytes]) -> Path:
    for name, data in files.items():
        path = base / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    return base


def test_concatenate_sorts_by_relative_path_string(tmp_path):
    # Path ordering would put a/x.hmm first; the posix string puts a-b/ first.
    workdir = _write(
        tmp_path / "work",
        {"a/x.hmm": b"a\n", "a-b/x.hmm": b"a-b\n", "B/x.hmm": b"B\n", "c/x.hmm": b"c\n"},
    )
    concatenate("*/*.hmm", workdir, tmp_path / "none", workdir / "all.hmm")
    assert (workdir / "all.hmm").read_bytes() == b"B\na-b\na\nc\n"


def test_concatenate_hmm_is_byte_exact(tmp_path):
    workdir = _write(
        tmp_path / "work",
        {"p/1.hmm": b"HEAD\r\nno newline", "p/2.hmm": b"HEAD\nsecond\n", "p/3.hmm": b""},
    )
    concatenate("p/*.hmm", workdir, tmp_path / "none", workdir / "profiles.hmm")
    assert (workdir / "profiles.hmm").read_bytes() == b"HEAD\r\nno newlineHEAD\nsecond\n"


def test_concatenate_tsv_drops_later_headers_and_inserts_newlines(tmp_path):
    workdir = _write(
        tmp_path / "work",
        {
            "t/1.tsv": b"id\tv\nA\t1",
            "t/2.tsv": b"id\tv\nB\t2\n",
            "t/3.tsv": b"id\tv\n",
            "t/4.tsv": b"id\tv\nC\t3\n",
        },
    )
    concatenate("t/*.tsv", workdir, tmp_path / "none", workdir / "all.tsv")
    assert (workdir / "all.tsv").read_bytes() == b"id\tv\nA\t1\nB\t2\nC\t3\n"


def test_concatenate_tsv_header_without_newline_in_a_later_file(tmp_path):
    workdir = _write(tmp_path / "work", {"t/1.tsv": b"id\nA", "t/2.tsv": b"id", "t/3.tsv": b"id\nB\n"})
    concatenate("t/*.tsv", workdir, tmp_path / "none", workdir / "all.tsv")
    assert (workdir / "all.tsv").read_bytes() == b"id\nA\nB\n"


def test_concatenate_never_reads_its_own_target(tmp_path):
    workdir = _write(tmp_path / "work", {"a.hmm": b"a\n", "all.hmm": b"stale\n"})
    concatenate("*.hmm", workdir, tmp_path / "none", workdir / "all.hmm")
    assert (workdir / "all.hmm").read_bytes() == b"a\n"


def test_concatenate_zero_matches_raises(tmp_path):
    workdir = _write(tmp_path / "work", {"p/1.txt": b"x"})
    downloads = _write(tmp_path / "downloads", {"q/1.hmm": b"x"})
    with pytest.raises(PrepareError, match="matched no files"):
        concatenate("p/*.hmm", workdir, downloads, workdir / "profiles.hmm")
    assert not (workdir / "profiles.hmm").exists()


@pytest.mark.parametrize("pattern", ["", "/abs/*.hmm", "../*.hmm"])
def test_concatenate_bad_pattern_raises(tmp_path, pattern):
    workdir = _write(tmp_path / "work", {"p/1.hmm": b"x"})
    with pytest.raises(PrepareError):
        concatenate(pattern, workdir, workdir, workdir / "profiles.hmm")


def test_concatenate_read_failure_leaves_no_partial_target(tmp_path):
    # A directory matching the glob fails mid-write, after p/1.hmm was copied.
    workdir = _write(tmp_path / "work", {"p/1.hmm": b"x\n"})
    (workdir / "p" / "2.hmm").mkdir()
    with pytest.raises(PrepareError, match="cannot concatenate"):
        concatenate("p/*.hmm", workdir, workdir, workdir / "profiles.hmm")
    assert _files_below(workdir) == {"p", "p/1.hmm", "p/2.hmm"}


# rename


def test_rename_is_byte_identical(tmp_path):
    data = bytes(range(256)) + b"\r\n\x00tail without newline"
    (tmp_path / "src.bin").write_bytes(data)
    rename(tmp_path / "src.bin", tmp_path / "sub" / "dst.bin")
    assert (tmp_path / "sub" / "dst.bin").read_bytes() == data
    assert (tmp_path / "src.bin").read_bytes() == data


def test_rename_missing_source_raises(tmp_path):
    with pytest.raises(PrepareError):
        rename(tmp_path / "absent", tmp_path / "dst")
    assert not (tmp_path / "dst").exists()


# pul-table

_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"


def _escape(text: str) -> str:
    # XML parsers turn a literal CR into LF; a character reference survives.
    return text.replace("&", "&amp;").replace("<", "&lt;").replace("\r", "&#13;")


def _s(ref, index):
    return f'<c r="{ref}" t="s"><v>{index}</v></c>'


def _inline(ref, text):
    return f'<c r="{ref}" t="inlineStr"><is><t>{_escape(text)}</t></is></c>'


def _v(ref, text):
    return f'<c r="{ref}"><v>{_escape(text)}</v></c>'


def _workbook(
    path: Path,
    rows: list[list[str]],
    *,
    shared=(),
    target: str = "worksheets/sheet1.xml",
    first_sheet: str = "rId1",
) -> Path:
    """A minimal workbook; ``rows`` hold already-encoded <c> elements."""
    sheet_data = "".join(
        f'<row r="{number}">{"".join(cells)}</row>' for number, cells in enumerate(rows, start=1)
    )
    strings = "".join(f"<si><t>{_escape(text)}</t></si>" for text in shared)
    members = {
        "xl/workbook.xml": (
            f'<workbook xmlns="{_NS}" xmlns:r="{_R}"><sheets>'
            f'<sheet name="first" sheetId="1" r:id="{first_sheet}"/>'
            '<sheet name="other" sheetId="2" r:id="rId9"/></sheets></workbook>'
        ),
        "xl/_rels/workbook.xml.rels": (
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'<Relationship Id="rId1" Target="{target}"/>'
            '<Relationship Id="rId9" Target="worksheets/other.xml"/></Relationships>'
        ),
        "xl/worksheets/sheet1.xml": f'<worksheet xmlns="{_NS}"><sheetData>{sheet_data}</sheetData></worksheet>',
        "xl/worksheets/other.xml": (
            f'<worksheet xmlns="{_NS}"><sheetData><row r="1">{_inline("A1", "id")}'
            f'{_inline("B1", "cazymes_predicted_dbcan")}</row><row r="2">{_inline("A2", "PUL9")}'
            f'{_inline("B2", "GH1")}</row></sheetData></worksheet>'
        ),
    }
    if shared:
        members["xl/sharedStrings.xml"] = f'<sst xmlns="{_NS}">{strings}</sst>'
    with zipfile.ZipFile(path, "w") as archive:
        for name, text in members.items():
            archive.writestr(name, text)
    return path


HEADER = [_inline("A1", " ID "), _inline("B1", "Substrate"), _inline("C1", "cazymes_predicted_dbcan")]


def _table(tmp_path, rows, **kwargs) -> list[str]:
    xlsx = _workbook(tmp_path / "pul.xlsx", rows, **kwargs)
    pul_table(xlsx, tmp_path / "pul_rules.tsv")
    data = (tmp_path / "pul_rules.tsv").read_bytes()
    assert data.endswith(b"\n") and b"\r" not in data
    lines = data.decode("utf-8").split("\n")
    assert lines[0] == PUL_TABLE_HEADER and lines[-1] == ""
    return lines[1:-1]


def test_pul_table_decodes_shared_inline_and_plain_cells(tmp_path):
    shared = ["PUL0001", "xylan", "GH5,CBM48"]
    rows = [
        HEADER,
        [_s("A2", 0), _s("B2", 1), _s("C2", 2)],
        [_inline("A3", "PUL0002"), _inline("B3", "starch"), _inline("C3", "GH13|CBM48")],
        [_v("A4", "PUL0003"), _v("B4", "42"), _v("C4", "GH5")],
    ]
    assert _table(tmp_path, rows, shared=shared) == [
        "PUL0001\txylan\tGH5,CBM48",
        "PUL0002\tstarch\tGH13|CBM48",
        "PUL0003\t42\tGH5",
    ]


def test_pul_table_missing_cells_are_empty(tmp_path):
    rows = [
        HEADER,
        # No B cell; then a row ending before the rule column.
        [_inline("A2", "PUL1"), _inline("C2", "GH5")],
        [_inline("A3", "PUL2"), _inline("B3", "pectin")],
        # A cell with neither <v> nor <is>, and a shared-string cell with empty <v>.
        [_inline("A4", "PUL3"), '<c r="B4"/>', '<c r="C4" t="s"><v></v></c>'],
    ]
    assert _table(tmp_path, rows) == ["PUL1\t\tGH5", "PUL2\tpectin\t", "PUL3\t\t"]


def test_pul_table_columns_come_from_cell_refs_not_position(tmp_path):
    rows = [
        [_inline("C1", "cazymes_predicted_dbcan"), _inline("A1", "pul_id"), _inline("AA1", "substrate")],
        [_inline("AA2", "chitin"), _inline("A2", "PUL7"), _inline("C2", "GH18")],
    ]
    assert _table(tmp_path, rows) == ["PUL7\tchitin\tGH18"]


def test_pul_table_keeps_only_pul_rows_and_normalises_ids(tmp_path):
    rows = [
        HEADER,
        [_inline("A2", "see notes"), _inline("C2", "GH5")],
        [_inline("A3", " pul0003 "), _inline("B3", "  beta-glucan "), _inline("C3", " GH5, GH999 ")],
        [_inline("A4", "PUL"), _inline("C4", "GH5")],
        [_inline("A5", "PUL12a"), _inline("C5", "GH5")],
        [_inline("A6", "XPUL1"), _inline("C6", "GH5")],
        [_inline("B7", "orphan substrate")],
        [],
        [_inline("A9", "Pul0010"), _inline("C9", "GH1")],
    ]
    assert _table(tmp_path, rows) == ["PUL0003\tbeta-glucan\tGH5, GH999", "PUL0010\t\tGH1"]


def test_pul_table_replaces_embedded_tabs_and_line_breaks(tmp_path):
    rows = [
        HEADER,
        [_inline("A2", "PUL1"), _inline("B2", "starch\nglycogen"), _inline("C2", "GH5\tCBM48")],
        [_inline("A3", "PUL2"), _inline("B3", "a\r\nb"), _inline("C3", "\nGH5\r\n")],
    ]
    assert _table(tmp_path, rows) == ["PUL1\tstarch glycogen\tGH5 CBM48", "PUL2\ta  b\tGH5"]


def test_pul_table_keeps_duplicate_ids_in_order(tmp_path):
    rows = [
        HEADER,
        [_inline("A2", "PUL1"), _inline("B2", "first"), _inline("C2", "GH5")],
        [_inline("A3", "PUL2"), _inline("B3", "other"), _inline("C3", "GH1")],
        [_inline("A4", "pul1"), _inline("B4", "second"), _inline("C4", "GH6")],
    ]
    assert _table(tmp_path, rows) == ["PUL1\tfirst\tGH5", "PUL2\tother\tGH1", "PUL1\tsecond\tGH6"]


@pytest.mark.parametrize("name", ["id", "PUL", "pul_id", "Pul ID"])
def test_pul_table_accepts_every_pul_column_name(tmp_path, name):
    rows = [
        [_inline("A1", name), _inline("B1", "cazymes_predicted_dbcan")],
        [_inline("A2", "PUL1"), _inline("B2", "GH5")],
    ]
    assert _table(tmp_path, rows) == ["PUL1\t\tGH5"]


def test_pul_table_prefers_id_over_other_pul_columns(tmp_path):
    rows = [
        [_inline("A1", "pul_id"), _inline("B1", "id"), _inline("C1", "cazymes_predicted_dbcan")],
        [_inline("A2", "PUL1"), _inline("B2", "PUL2"), _inline("C2", "GH5")],
    ]
    assert _table(tmp_path, rows) == ["PUL2\t\tGH5"]


def test_pul_table_substrate_is_first_header_containing_substrate(tmp_path):
    rows = [
        [
            _inline("A1", "id"),
            _inline("B1", "Predicted Substrate"),
            _inline("C1", "substrate_final"),
            _inline("D1", "cazymes_predicted_dbcan"),
        ],
        [_inline("A2", "PUL1"), _inline("B2", "xylan"), _inline("C2", "other"), _inline("D2", "GH5")],
    ]
    assert _table(tmp_path, rows) == ["PUL1\txylan\tGH5"]


def test_pul_table_repeated_header_resolves_to_its_last_column(tmp_path):
    # As phu main's parse_pul_workbook and phu-db-prep read it: names collapse
    # into a dict, so a header repeated after normalising takes its last column.
    rows = [
        [
            _inline("A1", "id"),
            _inline("B1", "substrate"),
            _inline("C1", "cazymes_predicted_dbcan"),
            _inline("D1", "Substrate"),
            _inline("E1", "ID"),
        ],
        [
            _inline("A2", "PUL1"),
            _inline("B2", "first"),
            _inline("C2", "GH5"),
            _inline("D2", "last"),
            _inline("E2", "PUL2"),
        ],
    ]
    assert _table(tmp_path, rows) == ["PUL2\tlast\tGH5"]


def test_pul_table_without_substrate_column_writes_empty_substrate(tmp_path):
    rows = [
        [_inline("A1", "id"), _inline("B1", "cazymes_predicted_dbcan")],
        [_inline("A2", "PUL1"), _inline("B2", "GH5")],
    ]
    assert _table(tmp_path, rows) == ["PUL1\t\tGH5"]


@pytest.mark.parametrize(
    "header, message",
    [
        ([_inline("A1", "name"), _inline("B1", "cazymes_predicted_dbcan")], "missing columns: pul"),
        ([_inline("A1", "id"), _inline("B1", "substrate")], "missing columns: cazymes_predicted_dbcan"),
    ],
    ids=["no-pul", "no-rule"],
)
def test_pul_table_missing_required_column_raises(tmp_path, header, message):
    xlsx = _workbook(tmp_path / "pul.xlsx", [header, [_inline("A2", "PUL1"), _inline("B2", "GH5")]])
    with pytest.raises(PrepareError, match=message):
        pul_table(xlsx, tmp_path / "pul_rules.tsv")
    assert _files_below(tmp_path) == {"pul.xlsx"}


def test_pul_table_reads_only_the_first_listed_sheet(tmp_path):
    rows = [HEADER, [_inline("A2", "PUL1"), _inline("C2", "GH5")]]
    assert _table(tmp_path, rows) == ["PUL1\t\tGH5"]
    # Listing the other sheet first switches to it.
    assert _table(tmp_path, rows, first_sheet="rId9") == ["PUL9\t\tGH1"]


@pytest.mark.parametrize("target", ["/xl/worksheets/sheet1.xml", "xl/worksheets/sheet1.xml"])
def test_pul_table_resolves_rels_targets_below_xl(tmp_path, target):
    rows = [HEADER, [_inline("A2", "PUL1"), _inline("C2", "GH5")]]
    assert _table(tmp_path, rows, target=target) == ["PUL1\t\tGH5"]


def test_pul_table_empty_sheet_raises(tmp_path):
    xlsx = _workbook(tmp_path / "pul.xlsx", [[]])
    with pytest.raises(PrepareError, match="is empty"):
        pul_table(xlsx, tmp_path / "pul_rules.tsv")


def _patch_zip_headers(data: bytes, damage: str) -> bytes:
    """Set the encryption flag or an unknown compression method on every member."""
    patched = bytearray(data)
    # Field offsets in the local (PK\3\4) and central (PK\1\2) headers.
    offsets = {"encrypted": (6, 8), "unknown-compression": (8, 10)}[damage]
    for signature, offset in zip((b"PK\x03\x04", b"PK\x01\x02"), offsets):
        start = patched.find(signature)
        while start != -1:
            field = slice(start + offset, start + offset + 2)
            value = int.from_bytes(patched[field], "little")
            value = value | 1 if damage == "encrypted" else 99
            patched[field] = value.to_bytes(2, "little")
            start = patched.find(signature, start + 4)
    return bytes(patched)


@pytest.mark.parametrize(
    "damage",
    ["not-a-zip", "missing-sheet", "bad-shared-index", "bad-xml", "encrypted", "unknown-compression"],
)
def test_pul_table_damaged_workbook_raises(tmp_path, damage):
    xlsx = tmp_path / "pul.xlsx"
    if damage == "not-a-zip":
        xlsx.write_bytes(b"PK not really")
    elif damage == "missing-sheet":
        _workbook(xlsx, [HEADER], target="worksheets/absent.xml")
    elif damage == "bad-shared-index":
        _workbook(xlsx, [HEADER, [_s("A2", 5)]], shared=["only"])
    elif damage in ("encrypted", "unknown-compression"):
        _workbook(xlsx, [HEADER, [_s("A2", 0)]], shared=["PUL1"])
        xlsx.write_bytes(_patch_zip_headers(xlsx.read_bytes(), damage))
    else:
        with zipfile.ZipFile(xlsx, "w") as archive:
            archive.writestr("xl/workbook.xml", "<workbook")
    dst = tmp_path / "pul_rules.tsv"
    dst.write_bytes(b"previous")

    with pytest.raises(PrepareError, match="not a readable xlsx workbook"):
        pul_table(xlsx, dst)
    assert dst.read_bytes() == b"previous"
    assert _files_below(tmp_path) == {"pul.xlsx", "pul_rules.tsv"}


def test_pul_table_matches_fixture_workbook(tmp_path):
    pul_table(FIXTURE_DIR / "upstream" / "dbcan" / "dbCAN-PUL.xlsx", tmp_path / "pul_rules.tsv")
    assert (tmp_path / "pul_rules.tsv").read_bytes() == (
        FIXTURE_DIR / "prepared" / "dbcan" / "pul_rules.tsv"
    ).read_bytes()


def test_extract_without_tarfile_data_filter(tmp_path, monkeypatch, fixture_manifest):
    # Python 3.10 and 3.11 patch releases before the filter use the manual write path.
    monkeypatch.delattr(tarfile, "data_filter", raising=False)
    entry = fixture_manifest.databases["kofam"]
    workdir = tmp_path / "work"
    run_steps(entry.steps, workdir=workdir, downloads=FIXTURE_DIR / "upstream" / "kofam")
    for item in entry.prepared_files:
        assert hashlib.sha256((workdir / item.path).read_bytes()).hexdigest() == item.sha256
