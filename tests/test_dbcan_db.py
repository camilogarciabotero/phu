import io
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from typer.testing import CliRunner

import phu.cli as cli_module
from phu import dbcan_db
from phu.cli import app
from phu.dbcan_db import (
    build_family_offsets,
    is_dbcan_id,
    normalize_dbcan_id,
    parse_pul_workbook,
)

runner = CliRunner()


def test_download_streams_to_atomic_destination(monkeypatch, tmp_path: Path):
    class StreamingResponse(io.BytesIO):
        def read(self, size=-1):
            assert size != -1
            return super().read(size)

    payload = b"dbCAN" * 1024
    monkeypatch.setattr(dbcan_db, "urlopen", lambda _: StreamingResponse(payload))
    destination = tmp_path / "database" / "dbcan.txt"

    dbcan_db._download("https://example.test/dbcan.txt", destination)

    assert destination.read_bytes() == payload


def _write_workbook(path: Path, rows: list[list[str]]) -> None:
    cells = []
    for row_number, row in enumerate(rows, start=1):
        row_cells = []
        for column_number, value in enumerate(row, start=1):
            column = chr(64 + column_number)
            row_cells.append(
                f'<c r="{column}{row_number}" t="inlineStr"><is><t>{value}</t></is></c>'
            )
        cells.append(f'<row r="{row_number}">' + "".join(row_cells) + "</row>")
    sheet = (
        '<?xml version="1.0"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        "<sheetData>" + "".join(cells) + "</sheetData></worksheet>"
    )
    workbook = (
        '<?xml version="1.0"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
        '<sheet name="PUL" sheetId="1" r:id="rId1"/></sheets></workbook>'
    )
    rels = (
        '<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Target="worksheets/sheet1.xml" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet"/></Relationships>'
    )
    with ZipFile(path, "w", ZIP_DEFLATED) as archive:
        archive.writestr("xl/workbook.xml", workbook)
        archive.writestr("xl/_rels/workbook.xml.rels", rels)
        archive.writestr("xl/worksheets/sheet1.xml", sheet)


def test_dbcan_family_names_normalize_hmm_suffix():
    assert is_dbcan_id("gh128.hmm")
    assert normalize_dbcan_id("gh128.hmm") == "GH128"


def test_pul0621_rule_is_header_selected_and_ordered(tmp_path: Path):
    hmm = tmp_path / "dbcan.txt"
    hmm.write_text(
        "HMMER3/f\nNAME GH133.hmm\n//\n"
        "HMMER3/f\nNAME GH3.hmm\n//\n"
        "HMMER3/f\nNAME GH57.hmm\n//\n"
        "HMMER3/f\nNAME GT4.hmm\n//\n"
    )
    workbook = tmp_path / "pul.xlsx"
    _write_workbook(
        workbook,
        [
            ["substrate_final", "cazymes_predicted_dbcan", "ID"],
            ["starch", "GH133|GH3,GH57|GT4", "PUL0621"],
        ],
    )

    rules = parse_pul_workbook(workbook, set(build_family_offsets(hmm)))

    assert rules["PUL0621"].families == ("GH133", "GH3", "GH57", "GT4")
    assert rules["PUL0621"].raw_rule == "GH133|GH3,GH57|GT4"
    assert rules["PUL0621"].substrate == "starch"


def test_pul_rule_preserves_unsupported_tokens(tmp_path: Path):
    workbook = tmp_path / "pul.xlsx"
    _write_workbook(
        workbook,
        [["PUL", "cazymes_predicted_dbcan"], ["PUL0001", "PL6_3,CBM35inCE17"]],
    )

    rule = parse_pul_workbook(workbook, set())["PUL0001"]

    assert rule.families == ()
    assert rule.unresolved_tokens == ("PL6_3", "CBM35INCE17")
    assert not rule.resolved


def test_pul_rule_marks_ancillary_profiles_unresolved(tmp_path: Path):
    workbook = tmp_path / "pul.xlsx"
    _write_workbook(
        workbook,
        [["ID", "cazymes_predicted_dbcan"], ["PUL0187", "GH16_21,SLH|CBM54"]],
    )

    rule = parse_pul_workbook(workbook, {"GH16_21", "CBM54"})["PUL0187"]

    assert rule.families == ("GH16_21", "CBM54")
    assert rule.unresolved_tokens == ("SLH",)
    assert not rule.resolved


def test_dbcan_cli_lifecycle_dispatch(monkeypatch):
    manifest = {
        "files": {
            "hmm": {"path": "/tmp/hmm"},
            "family_offsets": {"path": "/tmp/family_offsets"},
            "pul_rules": {"path": "/tmp/pul_rules"},
        }
    }
    monkeypatch.setattr(cli_module, "get_dbcan_database_status", lambda: {
        "downloaded": True,
        "indexed": True,
        "manifest_exists": True,
        "model_count": 4,
        "sparse_cached_count": 0,
        "root": "/tmp/dbcan",
    })
    monkeypatch.setattr(cli_module, "prepare_dbcan_database", lambda **_: manifest)
    monkeypatch.setattr(cli_module, "refresh_dbcan_database", lambda: manifest)
    monkeypatch.setattr(cli_module, "remove_dbcan_database", lambda: True)

    assert runner.invoke(app, ["dbs", "list"]).exit_code == 0
    assert "dbcan\tready" in runner.invoke(app, ["dbs", "list"]).stdout
    assert runner.invoke(app, ["dbs", "status", "dbcan"]).exit_code == 0
    assert runner.invoke(app, ["dbs", "prepare", "dbcan"]).exit_code == 0
    assert runner.invoke(app, ["dbs", "refresh", "dbcan"]).exit_code == 0
    assert runner.invoke(app, ["dbs", "remove", "dbcan", "--yes"]).exit_code == 0