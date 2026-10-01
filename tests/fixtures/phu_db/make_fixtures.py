"""Generate the miniature phu-db fixture: fake upstreams, prepared tree and manifest.

Run from the repository root to rewrite the committed fixture:

    uv run python tests/fixtures/phu_db/make_fixtures.py

``--out DIR`` writes the same files somewhere else, which is how the tests check
that the committed copy is reproducible. Output is byte-deterministic: HMMs are
built by pyhmmer with a fixed seed and no creation time, gzip members use stored
deflate blocks with mtime 0 (so no zlib version can change them), and tar and
zip members carry fixed metadata in sorted order.

Only the stdlib and pyhmmer are used. The prepared files are produced by running
each database's ``steps`` over the upstream bytes with the verb semantics of the
phu-db 0.1.1 contract, so they double as the expected output of phu's verbs.
"""

from __future__ import annotations

import argparse
import fnmatch
import gzip
import hashlib
import io
import json
import re
import shutil
import struct
import tarfile
import zipfile
import zlib
from pathlib import Path, PurePosixPath
from xml.etree import ElementTree

from pyhmmer.easel import Alphabet, TextMSA, TextSequence
from pyhmmer.plan7 import Background, Builder

HERE = Path(__file__).resolve().parent
BUNDLE_VERSION = "phu-db-0.1.1-fixture"
# The committed HMMs embed this pyhmmer/HMMER build's version line and numbers;
# another pyhmmer writes different bytes, so the reproducibility test skips.
PYHMMER_VERSION = "0.12.3"
URL_ROOT = "https://phu-db-fixture.invalid"

# Real header lines of the CheckAMG v1.1.1 tables (and of the 0.1.0 bundle).
VSCORES_HEADER = "id\tname\tV-score\tVL-score\tdb"
AMG_HEADER = (
    "id\tdb\tname\tmetabolic_paths\ttotal_paths\tpath_source\tmetabolic_ratio"
    "\tVL-score\tAMG_weight\tamg_level"
)
APG_HEADER = "id\tV-score\tVL-score\tdb\tname"
FILTERS_HEADER = (
    "id\tdb\tname\tfilter_essential\tfilter_glucan\tfilter_lipid\tfilter_methyl"
    "\tfilter_nucleotide"
)
KO_LIST_HEADER = (
    "knum\tthreshold\tscore_type\tprofile_type\tF-measure\tnseq\tnseq_used\talen"
    "\tmlen\teff_nseq\tre/pos\tdefinition"
)

# (name, accession, description, aligned sequences); every model is synthetic.
PFAM_MODELS = (
    ("FixtureA", "PF00001.1", "phu fixture domain A", ("MKVLAAGIVGLLLAQ", "MKVLSAGIVGLLIAQ", "MRVLAAGLVGLLLAE")),
    ("FixtureB", "PF00002.1", "phu fixture domain B", ("GSWTKHEPLRCDNYF", "GSWTRHEPLKCDNYF", "GAWTKHDPLRCENYF")),
)
KOFAM_MODELS = (
    ("K00001", None, None, ("MSTAKLIEHGRWQAV", "MSTAKVIEHGRWQSV", "MATSKLIDHGRWQAV")),
    ("K00002", None, None, ("PELNDVYCRTAGHKM", "PELNEVYCRTSGHKM", "PDLNDVFCRTAGHRM")),
    ("K00003", None, None, ("WQHFRGTIKPMSNLE", "WQHYRGTVKPMSNLE", "WEHFRGSIKPLSNLD")),
)
DBCAN_MODELS = (
    ("GH5.hmm", None, None, ("FTDGWNYRSLEKMPA", "FTDGWNYKSLEKLPA", "FSDGWNFRSIEKMPA")),
    ("CBM48.hmm", None, None, ("KVYLTGSFNNWEPIA", "KVYLAGSFNDWQPIA", "RVHLTGSFNNWEPLS")),
)

KO_LIST_ROWS = (
    "K00001\t120.50\tdomain\tall\t0.883405\t97\t92\t560\t15\t1.26\t0.589\tfixture dehydrogenase A [EC:1.1.1.1]",
    "K00002\t-\t-\t-\t-\t1\t1\t15\t15\t0.39\t0.593\tfixture monooxygenase without threshold",
    "K00003\t80.25\tfull\tall\t0.512000\t42\t40\t300\t15\t2.10\t0.601\tfixture transferase C",
)
SUBSTRATE_ROWS = ("PUL0001\txylan", "PUL0002\tstarch", "PUL0003\tbeta-glucan")
# Header plus rows for the PUL workbook. The workbook encodes these values with
# a mix of shared strings, inline strings and plain values (see _xlsx).
PUL_HEADER = ("ID", "organism_name", "substrate", "num_cazymes", "cazymes_predicted_dbcan")
PUL_ROWS = (
    ("PUL0001", "Fixture bacterium alpha", "xylan", "2", "GH5,CBM48"),
    ("PUL0002", "Fixture bacterium beta", "starch\nglycogen", "2", "GH5|CBM48"),
    ("see notes", "", "", "", ""),
    (" pul0003 ", "Fixture bacterium gamma", "beta-glucan", "2", "GH5, GH999"),
)

VSCORES_ROWS = (
    "K00001\tfixture dehydrogenase A\t10.0\t4.890951685\tKEGG",
    "K00002\tfixture monooxygenase\t3.5\t2.104233101\tKEGG",
    "PF00001.1\tphu fixture domain A\t1.25\t1.502991400\tPfam",
    "PF00002.1\tphu fixture domain B\t0.5\t0.811402377\tPfam",
    "K00003\tfixture transferase C\t\t\tKEGG",
)
AMG_ROWS = (
    "K00001\tKEGG\tfixture dehydrogenase A [EC:1.1.1.1]\t5.0\t5.0\tKEGG\t1.0\t\t1.0\tvery high",
    "K00002\tKEGG\tfixture monooxygenase\t2.0\t4.0\tKEGG\t0.5\t2.104233101\t0.5\tmedium",
    "K00003\tKEGG\tfixture transferase C\t1.0\t5.0\tKEGG\t0.2\t\t0.2\tlow",
    "PF00001.1\tPfam\tphu fixture domain A\t4.0\t4.0\tGO\t1.0\t0.0\t1.0\tvery high",
    "PF00002.1\tPfam\tphu fixture domain B\t1.0\t2.0\tGO\t0.5\t0.811402377\t0.5\tmedium",
)
APG_ROWS = (
    "K00001\t1.33\t2.123851641\tKEGG\tfixture dehydrogenase A",
    "K00002\t1.49\t2.173186268\tKEGG\tfixture monooxygenase",
    "K00003\t\t\tKEGG\tfixture transferase C",
    "PF00001.1\t\t\tPfam\tphu fixture domain A",
    "PF00002.1\t0.02\t0.410000000\tPfam\tphu fixture domain B",
)
AREG_ROWS = (
    "K00001\t10.0\t3.148910993\tKEGG\tfixture dehydrogenase A",
    "K00002\t\t\tKEGG\tfixture monooxygenase",
    "K00003\t0.01\t0.120000000\tKEGG\tfixture transferase C",
    "PF00001.1\t1.0\t1.000000000\tPfam\tphu fixture domain A",
    "PF00002.1\t\t\tPfam\tphu fixture domain B",
)
FILTER_ROWS = {
    "AMG_filters.tsv": (
        "K00001\tKEGG\tfixture dehydrogenase A\tFalse\tFalse\tTrue\tFalse\tFalse",
        "K00002\tKEGG\tfixture monooxygenase\tFalse\tTrue\tFalse\tFalse\tFalse",
        "K00003\tKEGG\tfixture transferase C\tTrue\tFalse\tFalse\tFalse\tFalse",
        "PF00001.1\tPfam\tphu fixture domain A\tTrue\tFalse\tFalse\tFalse\tFalse",
        "PF00002.1\tPfam\tphu fixture domain B\tFalse\tFalse\tFalse\tTrue\tFalse",
    ),
    "APG_filters.tsv": (
        "K00001\tKEGG\tfixture dehydrogenase A\tFalse\tFalse\tFalse\tFalse\tTrue",
        "K00002\tKEGG\tfixture monooxygenase\tTrue\tFalse\tFalse\tFalse\tFalse",
        "K00003\tKEGG\tfixture transferase C\tFalse\tFalse\tFalse\tTrue\tFalse",
        "PF00001.1\tPfam\tphu fixture domain A\tFalse\tTrue\tFalse\tFalse\tFalse",
        "PF00002.1\tPfam\tphu fixture domain B\tFalse\tFalse\tTrue\tFalse\tFalse",
    ),
    "AReG_filters.tsv": (
        "K00001\tKEGG\tfixture dehydrogenase A\tTrue\tFalse\tFalse\tFalse\tFalse",
        "K00002\tKEGG\tfixture monooxygenase\tFalse\tFalse\tFalse\tFalse\tFalse",
        "K00003\tKEGG\tfixture transferase C\tFalse\tTrue\tFalse\tFalse\tFalse",
        "PF00001.1\tPfam\tphu fixture domain A\tFalse\tFalse\tFalse\tFalse\tTrue",
        "PF00002.1\tPfam\tphu fixture domain B\tFalse\tFalse\tTrue\tFalse\tFalse",
    ),
}

# Per-database manifest metadata, steps and prepared files, in manifest order.
DATABASES = {
    "pfam": {
        "provenance": "fetch",
        "redistribution": "pin-only",
        "upstream_version": "37.0",
        "license": "CC0-1.0",
        "citation": "Mistry J, et al. Nucleic Acids Res. 2021;49(D1):D412-D419.",
        "sources": ("Pfam-A.hmm.gz",),
        "steps": ("gunzip:Pfam-A.hmm.gz -> profiles.hmm",),
        "prepared_files": (("profiles.hmm", "hmm"),),
        "build": None,
    },
    "kofam": {
        "provenance": "fetch",
        "redistribution": "pin-only",
        "upstream_version": "genome.jp archives/2026-07-02",
        "license": "UNKNOWN",
        "citation": "KOfamKOALA",
        "sources": ("profiles.tar.gz", "ko_list.gz"),
        "steps": (
            "extract:profiles.tar.gz -> profiles/",
            "concatenate:profiles/*.hmm -> profiles.hmm",
            "gunzip:ko_list.gz -> ko_list.tsv",
        ),
        "prepared_files": (("profiles.hmm", "hmm"), ("ko_list.tsv", "tsv")),
        "build": None,
    },
    "dbcan": {
        "provenance": "fetch",
        "redistribution": "pin-only",
        "upstream_version": "V15",
        "license": "UNKNOWN",
        "citation": "dbCAN",
        "sources": ("dbCAN-HMMdb-V15.txt", "dbCAN-PUL-substrate-metadata.txt", "dbCAN-PUL.xlsx"),
        "steps": (
            "rename:dbCAN-HMMdb-V15.txt -> profiles.hmm",
            "rename:dbCAN-PUL-substrate-metadata.txt -> family_substrate.tsv",
            "pul-table:dbCAN-PUL.xlsx -> pul_rules.tsv",
        ),
        "prepared_files": (
            ("profiles.hmm", "hmm"),
            ("family_substrate.tsv", "tsv"),
            ("pul_rules.tsv", "tsv"),
        ),
        "build": None,
    },
    "avg": {
        "provenance": "build",
        "redistribution": "pin-only",
        "upstream_version": "v1.1.1",
        "license": "UNKNOWN",
        "citation": "CheckAMG",
        "sources": (
            "AMGs.tsv",
            "APGs.tsv",
            "AReGs.tsv",
            "AMG_filters.tsv",
            "APG_filters.tsv",
            "AReG_filters.tsv",
            "vscores.tsv",
        ),
        "steps": (
            "rename:vscores.tsv -> v_scores.tsv",
            "rename:AMGs.tsv -> positive_amg.tsv",
            "rename:APGs.tsv -> positive_apg.tsv",
            "rename:AReGs.tsv -> positive_areg.tsv",
            "rename:AMG_filters.tsv -> filters_amg.tsv",
            "rename:APG_filters.tsv -> filters_apg.tsv",
            "rename:AReG_filters.tsv -> filters_areg.tsv",
        ),
        "prepared_files": (
            ("v_scores.tsv", "tsv"),
            ("positive_amg.tsv", "tsv"),
            ("positive_apg.tsv", "tsv"),
            ("positive_areg.tsv", "tsv"),
            ("filters_amg.tsv", "tsv"),
            ("filters_apg.tsv", "tsv"),
            ("filters_areg.tsv", "tsv"),
        ),
        "build": {
            "source_repo": "https://github.com/AnantharamanLab/CheckAMG",
            "source_tag": "v1.1.1",
            "tool_version_key": "phu.avg_reference_db",
        },
    },
}

_XLSX_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"


def _table(header: str | None, rows: tuple[str, ...]) -> bytes:
    lines = ([header] if header else []) + list(rows)
    return ("\n".join(lines) + "\n").encode()


def _hmms(models: tuple) -> bytes:
    alphabet = Alphabet.amino()
    background = Background(alphabet)
    output = io.BytesIO()
    for name, accession, description, sequences in models:
        msa = TextMSA(
            name=name.encode(),
            sequences=[
                TextSequence(name=f"seq{index}".encode(), sequence=sequence)
                for index, sequence in enumerate(sequences, start=1)
            ],
        ).digitize(alphabet)
        hmm, _, _ = Builder(alphabet, seed=42).build_msa(msa, background)
        if accession:
            hmm.accession = accession.encode()
        if description:
            hmm.description = description.encode()
        # Without these the text carries the build clock and this script's path.
        hmm.creation_time = None
        hmm.command_line = None
        hmm.write(output)
    return output.getvalue()


def _gzip(data: bytes) -> bytes:
    """Gzip with stored deflate blocks: no compressor version can change the bytes."""
    chunks = [data[start : start + 0xFFFF] for start in range(0, len(data), 0xFFFF)] or [b""]
    blocks = b"".join(
        struct.pack("<BHH", int(index == len(chunks) - 1), len(chunk), len(chunk) ^ 0xFFFF) + chunk
        for index, chunk in enumerate(chunks)
    )
    header = b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\xff"
    return header + blocks + struct.pack("<II", zlib.crc32(data), len(data) & 0xFFFFFFFF)


def _tar_gz(directories: tuple[str, ...], files: dict[str, bytes]) -> bytes:
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w", format=tarfile.USTAR_FORMAT) as archive:
        entries = [(name, None) for name in directories] + list(files.items())
        for name, data in sorted(entries):
            info = tarfile.TarInfo(name)
            info.mtime = 0
            info.uid = info.gid = 0
            info.uname = info.gname = "root"
            if data is None:
                info.type = tarfile.DIRTYPE
                info.mode = 0o755
                archive.addfile(info)
            else:
                info.mode = 0o644
                info.size = len(data)
                archive.addfile(info, io.BytesIO(data))
    return _gzip(buffer.getvalue())


def _cell(ref: str, value: str, shared: list[str]) -> str:
    """Encode a cell the ways real workbooks do, varying by column.

    Shared strings by default, inline strings in column B, bare numbers, and
    plain-valued formula strings (t="str") in column E below the header.
    """
    escaped = value.replace("&", "&amp;").replace("<", "&lt;")
    if value.isdigit():
        return f'<c r="{ref}"><v>{value}</v></c>'
    if ref.startswith("B"):
        return f'<c r="{ref}" t="inlineStr"><is><t>{escaped}</t></is></c>'
    if ref.startswith("E") and ref != "E1":
        return f'<c r="{ref}" t="str"><f>CONCAT("")</f><v>{escaped}</v></c>'
    if value not in shared:
        shared.append(value)
    return f'<c r="{ref}" t="s"><v>{shared.index(value)}</v></c>'


def _xlsx() -> bytes:
    shared: list[str] = []
    rows = []
    for number, values in enumerate((PUL_HEADER, *PUL_ROWS), start=1):
        # Empty cells are left out, as spreadsheet writers do.
        cells = "".join(
            _cell(f"{chr(65 + column)}{number}", value, shared)
            for column, value in enumerate(values)
            if value
        )
        rows.append(f'<row r="{number}">{cells}</row>')
    sheet = (
        f'<?xml version="1.0" encoding="UTF-8"?><worksheet xmlns="{_XLSX_NS}">'
        f"<sheetData>{''.join(rows)}</sheetData></worksheet>"
    )
    notes = (
        f'<?xml version="1.0" encoding="UTF-8"?><worksheet xmlns="{_XLSX_NS}"><sheetData>'
        '<row r="1"><c r="A1" t="inlineStr"><is><t>not the PUL table</t></is></c></row>'
        "</sheetData></worksheet>"
    )
    strings = "".join(
        f'<si><t xml:space="preserve">{value.replace("&", "&amp;").replace("<", "&lt;")}</t></si>'
        for value in shared
    )
    members = {
        "[Content_Types].xml": (
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            '<Override PartName="/xl/worksheets/sheet1.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            '<Override PartName="/xl/worksheets/sheet2.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            '<Override PartName="/xl/sharedStrings.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml"/>'
            "</Types>"
        ),
        "_rels/.rels": (
            f'<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="{_PKG_REL_NS}">'
            f'<Relationship Id="rId1" Type="{_REL_NS}/officeDocument" Target="xl/workbook.xml"/>'
            "</Relationships>"
        ),
        "xl/workbook.xml": (
            f'<?xml version="1.0" encoding="UTF-8"?><workbook xmlns="{_XLSX_NS}" xmlns:r="{_REL_NS}">'
            '<sheets><sheet name="PUL" sheetId="1" r:id="rId1"/>'
            '<sheet name="Notes" sheetId="2" r:id="rId2"/></sheets></workbook>'
        ),
        "xl/_rels/workbook.xml.rels": (
            f'<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="{_PKG_REL_NS}">'
            f'<Relationship Id="rId1" Type="{_REL_NS}/worksheet" Target="worksheets/sheet1.xml"/>'
            f'<Relationship Id="rId2" Type="{_REL_NS}/worksheet" Target="worksheets/sheet2.xml"/>'
            f'<Relationship Id="rId3" Type="{_REL_NS}/sharedStrings" Target="sharedStrings.xml"/>'
            "</Relationships>"
        ),
        "xl/sharedStrings.xml": (
            f'<?xml version="1.0" encoding="UTF-8"?><sst xmlns="{_XLSX_NS}" '
            f'count="{len(shared)}" uniqueCount="{len(shared)}">{strings}</sst>'
        ),
        "xl/worksheets/sheet1.xml": sheet,
        "xl/worksheets/sheet2.xml": notes,
    }
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_STORED) as archive:
        for name in sorted(members):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.create_system = 3
            info.external_attr = 0o644 << 16
            archive.writestr(info, members[name].encode())
    return buffer.getvalue()


def upstream_files() -> dict[str, dict[str, bytes]]:
    filters = {name: _table(FILTERS_HEADER, rows) for name, rows in FILTER_ROWS.items()}
    kofam_models = {
        f"profiles/{model[0]}.hmm": _hmms((model,)) for model in KOFAM_MODELS
    }
    return {
        "pfam": {"Pfam-A.hmm.gz": _gzip(_hmms(PFAM_MODELS))},
        "kofam": {
            "profiles.tar.gz": _tar_gz(("profiles",), kofam_models),
            "ko_list.gz": _gzip(_table(KO_LIST_HEADER, KO_LIST_ROWS)),
        },
        "dbcan": {
            "dbCAN-HMMdb-V15.txt": _hmms(DBCAN_MODELS),
            "dbCAN-PUL-substrate-metadata.txt": _table(None, SUBSTRATE_ROWS),
            "dbCAN-PUL.xlsx": _xlsx(),
        },
        "avg": {
            "AMGs.tsv": _table(AMG_HEADER, AMG_ROWS),
            "APGs.tsv": _table(APG_HEADER, APG_ROWS),
            "AReGs.tsv": _table(APG_HEADER, AREG_ROWS),
            **filters,
            "vscores.tsv": _table(VSCORES_HEADER, VSCORES_ROWS),
        },
    }


def _xlsx_rows(data: bytes) -> list[list[str]]:
    """Decode the first worksheet exactly like phu's dbcan_db._xlsx_values (phu main)."""
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
            shared = ["".join(node.itertext()) for node in root.findall(f"{{{_XLSX_NS}}}si")]
        workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
        rels = ElementTree.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {item.attrib["Id"]: item.attrib["Target"] for item in rels}
        sheet = workbook.find(f".//{{{_XLSX_NS}}}sheet")
        if sheet is None:
            raise ValueError("PUL workbook has no worksheets")
        target = targets[sheet.attrib[f"{{{_REL_NS}}}id"]]
        sheet_path = "xl/" + target.removeprefix("/") if not target.startswith("xl/") else target
        root = ElementTree.fromstring(archive.read(sheet_path))
    rows = []
    for row in root.findall(f".//{{{_XLSX_NS}}}row"):
        values: dict[int, str] = {}
        for cell in row.findall(f"{{{_XLSX_NS}}}c"):
            match = re.search(r"[A-Z]+", cell.attrib.get("r", ""))
            if not match:
                continue
            column = 0
            for character in match.group():
                column = column * 26 + ord(character) - 64
            value = cell.find(f"{{{_XLSX_NS}}}v")
            text = "" if value is None else value.text or ""
            if cell.attrib.get("t") == "s" and text:
                text = shared[int(text)]
            elif cell.attrib.get("t") == "inlineStr":
                inline = cell.find(f"{{{_XLSX_NS}}}is")
                text = "" if inline is None else "".join(inline.itertext())
            values[column - 1] = text
        if values:
            rows.append([values.get(index, "") for index in range(max(values) + 1)])
    return rows


def _pul_table(data: bytes) -> bytes:
    rows = _xlsx_rows(data)
    header = [value.strip().lower() for value in rows[0]]
    pul = next(header.index(name) for name in ("id", "pul", "pul_id", "pul id") if name in header)
    rule = header.index("cazymes_predicted_dbcan")
    substrate = next((index for index, name in enumerate(header) if "substrate" in name), None)

    def value(row: list[str], index: int | None) -> str:
        text = row[index] if index is not None and index < len(row) else ""
        return re.sub(r"[\t\r\n]", " ", text.strip())

    lines = ["pul_id\tsubstrate\tcazymes_predicted_dbcan"]
    for row in rows[1:]:
        pul_id = value(row, pul).upper()
        if re.fullmatch(r"PUL\d+", pul_id):
            lines.append(f"{pul_id}\t{value(row, substrate)}\t{value(row, rule)}")
    return ("\n".join(lines) + "\n").encode()


def run_steps(steps: tuple[str, ...], downloads: dict[str, bytes]) -> dict[str, bytes]:
    """Apply the contract verbs to in-memory files; returns the workdir contents."""
    workdir: dict[str, bytes] = {}
    for step in steps:
        verb, source, target = re.fullmatch(
            r"(gunzip|extract|concatenate|rename|pul-table):(\S+) -> (\S+)", step
        ).groups()
        data = workdir.get(source, downloads.get(source))
        if verb == "gunzip":
            workdir[target] = gzip.decompress(data)
        elif verb == "extract":
            with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as archive:
                for member in archive.getmembers():
                    if member.isfile():
                        workdir[member.name] = archive.extractfile(member).read()
        elif verb == "concatenate":
            depth = len(PurePosixPath(source).parts)
            matches = sorted(
                name
                for name in workdir
                if len(PurePosixPath(name).parts) == depth and fnmatch.fnmatchcase(name, source)
            )
            workdir[target] = b"".join(workdir[name] for name in matches)
        elif verb == "rename":
            workdir[target] = data
        else:
            workdir[target] = _pul_table(data)
    return workdir


def _count(data: bytes, kind: str) -> int:
    if kind == "hmm":
        return sum(1 for line in data.splitlines() if line.startswith(b"NAME "))
    return sum(1 for line in data.splitlines() if line.strip())


def build(out: Path) -> None:
    upstream = upstream_files()
    databases = {}
    for name, spec in DATABASES.items():
        downloads = upstream[name]
        workdir = run_steps(spec["steps"], downloads)
        prepared_files = []
        for path, kind in spec["prepared_files"]:
            data = workdir[path]
            prepared_files.append(
                {
                    "path": path,
                    "format": kind,
                    "bytes": len(data),
                    "sha256": hashlib.sha256(data).hexdigest(),
                    "count": _count(data, kind),
                }
            )
            target = out / "prepared" / name / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
        for source in spec["sources"]:
            target = out / "upstream" / name / source
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(downloads[source])
        databases[name] = {
            "provenance": spec["provenance"],
            "redistribution": spec["redistribution"],
            "upstream_version": spec["upstream_version"],
            "license": spec["license"],
            "citation": spec["citation"],
            "layout_dir": name,
            "sources": [
                {
                    "name": source,
                    "url": f"{URL_ROOT}/{name}/{source}",
                    "sha256": hashlib.sha256(downloads[source]).hexdigest(),
                    "bytes": len(downloads[source]),
                }
                for source in spec["sources"]
            ],
            "steps": list(spec["steps"]),
            "prepared_files": prepared_files,
            "validation": {"passed": True, "reason": None},
            "build": spec["build"],
        }
    payload = {
        "bundle_version": BUNDLE_VERSION,
        "layout_version": 1,
        "generator": {"name": "phu-db-prep", "version": BUNDLE_VERSION.removeprefix("phu-db-")},
        "bundle": {"name": "phu-db.tar.gz", "layout_root": "phu-db", "url": None},
        "zenodo": {"concept_doi": None, "version_doi": None, "record_url": None},
        "databases": databases,
    }
    # Bytes, not text mode, so no platform newline translation can change the digest.
    (out / "phu-db.json").write_bytes((json.dumps(payload, indent=2) + "\n").encode())


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=HERE, help="fixture directory to (re)write")
    args = parser.parse_args(argv)
    # Only the generated subtrees are replaced; scripts beside them are kept.
    for subtree in ("upstream", "prepared"):
        shutil.rmtree(args.out / subtree, ignore_errors=True)
    args.out.mkdir(parents=True, exist_ok=True)
    build(args.out)


if __name__ == "__main__":
    main()
