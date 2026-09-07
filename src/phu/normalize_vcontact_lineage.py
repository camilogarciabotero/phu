from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

import pandas as pd

CANDIDATE_DELIM = "||"

RANKS: tuple[str, ...] = (
    "realm",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "subfamily",
    "genus",
)
CODES: dict[str, str] = {
    "realm": "NR",
    "kingdom": "NK",
    "phylum": "NP",
    "class": "NC",
    "order": "NO",
    "family": "NF",
    "subfamily": "NSF",
    "genus": "NG",
}
RANK_INDEX = {rank: index for index, rank in enumerate(RANKS)}

_NODE_RE = re.compile(r"^novel_(" + "|".join(RANKS) + r")_([0-9]+)$")
_COMPACT_RE = re.compile(r"^[^:]+(?::N(?:SF|R|K|P|C|O|F|G)[0-9]+)+$")
_COMPACT_CODE_RE = re.compile(r":(N(?:SF|R|K|P|C|O|F|G))[0-9]+")
_CODE_TO_RANK = {code: rank for rank, code in CODES.items()}
MISSING_SENTINELS = frozenset({"", "-", "<NA>"})


class Node(NamedTuple):
    rank: str
    identifier: str


@dataclass(frozen=True)
class ParseResult:
    status: str
    value: str | None
    rank_mismatch: bool = False
    deepest_rank: str | None = None
    reason: str = ""


@dataclass
class QualitySummary:
    cells_examined: int = 0
    normalized: int = 0
    unchanged: int = 0
    missing: int = 0
    unparsed: int = 0
    rank_mismatch: int = 0
    multi_candidate: int = 0
    lineage_skipped_unparsed: int = 0


@dataclass
class NormalizeLineageConfig:
    input_file: Path
    output_file: Path
    add_lineage: bool = False
    lineage_col: str = "compact_lineage"
    sep: str | None = None
    strict: bool = False


def _column_key(name: object) -> str:
    return str(name).strip().lower().replace("-", "_").replace(" ", "_")


def prediction_columns(columns: list[object]) -> dict[str, str]:
    found: dict[str, str] = {}
    for original in columns:
        key = _column_key(original)
        if (
            key.endswith("_prediction")
            and key.removesuffix("_prediction") in RANK_INDEX
        ):
            if key in found:
                raise ValueError(f"Duplicate normalized column name: {key}")
            found[key] = str(original)
    return found


def normalize_cell(
    value: object, column_rank: str | None, summary: QualitySummary
) -> object:
    summary.cells_examined += 1
    if pd.isna(value):
        summary.missing += 1
        return pd.NA

    text = str(value).strip()
    parts = text.split(CANDIDATE_DELIM)
    if len(parts) > 1:
        summary.multi_candidate += 1

    results = [parse_candidate(part, column_rank) for part in parts]
    if any(result.status == "unparsed" for result in results):
        summary.unparsed += 1
    elif all(result.status == "missing" for result in results):
        summary.missing += 1
        return pd.NA
    elif any(result.status == "normalized" for result in results):
        summary.normalized += 1
    else:
        summary.unchanged += 1

    if any(result.rank_mismatch for result in results):
        summary.rank_mismatch += 1

    values = [result.value for result in results]
    if all(value is None for value in values):
        return pd.NA
    return CANDIDATE_DELIM.join(value or "" for value in values)


def normalize_dataframe(
    dataframe: pd.DataFrame,
    *,
    add_lineage: bool = False,
    lineage_col: str = "compact_lineage",
) -> tuple[pd.DataFrame, QualitySummary]:
    columns = prediction_columns(list(dataframe.columns))
    summary = QualitySummary()
    result = dataframe.copy()
    for key, original in columns.items():
        rank = key.removesuffix("_prediction")
        result[original] = result[original].map(
            lambda value, rank=rank: normalize_cell(value, rank, summary)
        )

    if add_lineage:
        ordered = [
            columns[f"{rank}_prediction"]
            for rank in reversed(RANKS)
            if f"{rank}_prediction" in columns
        ]

        def choose(row: pd.Series) -> object:
            for column in ordered:
                value = row[column]
                if pd.isna(value):
                    continue
                text = str(value).strip()
                if not text:
                    continue
                if any(
                    parse_candidate(part, None).status == "unparsed"
                    for part in text.split(CANDIDATE_DELIM)
                ):
                    summary.lineage_skipped_unparsed += 1
                    continue
                return value
            return pd.NA

        result[lineage_col] = result.apply(choose, axis=1)
    return result, summary


def read_table(path: Path, sep: str | None = None) -> pd.DataFrame:
    delimiter = sep or ("\t" if path.suffix.lower() == ".tsv" else ",")
    return pd.read_csv(
        path, sep=delimiter, dtype=str, keep_default_na=False, na_filter=False
    )


def write_table(dataframe: pd.DataFrame, path: Path, sep: str | None = None) -> None:
    delimiter = sep or ("\t" if path.suffix.lower() == ".tsv" else ",")
    dataframe.to_csv(path, sep=delimiter, index=False, na_rep="")


def normalize_lineage_file(config: NormalizeLineageConfig) -> QualitySummary:
    dataframe = read_table(config.input_file, config.sep)
    normalized, summary = normalize_dataframe(
        dataframe, add_lineage=config.add_lineage, lineage_col=config.lineage_col
    )
    if config.strict and (summary.unparsed or summary.rank_mismatch):
        raise ValueError(
            "Strict mode rejected unparsed or rank-mismatched lineage values."
        )
    config.output_file.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        import tempfile

        with tempfile.NamedTemporaryFile(
            mode="w",
            suffix=config.output_file.suffix,
            dir=config.output_file.parent,
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
        write_table(normalized, temporary_path, config.sep)
        temporary_path.replace(config.output_file)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return summary


def parse_candidate(raw: str, column_rank: str | None) -> ParseResult:
    text = raw.strip()

    if text in MISSING_SENTINELS:
        return ParseResult("missing", None)

    if _COMPACT_RE.fullmatch(text):
        codes = _COMPACT_CODE_RE.findall(text)
        ranks = [_CODE_TO_RANK[code] for code in codes]
        indices = [RANK_INDEX[rank] for rank in ranks]
        if indices != sorted(set(indices)):
            return ParseResult("unparsed", text, reason="compact codes out of order")
        deepest_rank = ranks[-1]
        mismatch = column_rank is not None and deepest_rank != column_rank
        return ParseResult("unchanged", text, mismatch, deepest_rank=deepest_rank)

    segments = text.split("_of_")
    nodes: list[Node] = []
    index = 0
    while index < len(segments):
        match = _NODE_RE.fullmatch(segments[index])
        if match is None:
            break
        nodes.append(Node(match.group(1), match.group(2)))
        index += 1

    if not nodes:
        if "novel_" in text or ":" in text:
            return ParseResult(
                "unparsed", text, reason="novel_ present but unparseable"
            )
        return ParseResult("unchanged", text)

    anchor = "_of_".join(segments[index:])
    if not anchor:
        return ParseResult("unparsed", text, reason="anchorless chain")
    if ":" in anchor:
        return ParseResult("unparsed", text, reason="anchor contains ':'")

    indices = [RANK_INDEX[node.rank] for node in nodes]
    if any(left <= right for left, right in zip(indices, indices[1:])):
        return ParseResult("unparsed", text, reason="rank order invalid")

    rendered = anchor + "".join(
        f":{CODES[node.rank]}{node.identifier}" for node in reversed(nodes)
    )
    deepest_rank = nodes[0].rank
    mismatch = column_rank is not None and deepest_rank != column_rank
    return ParseResult("normalized", rendered, mismatch, deepest_rank)
