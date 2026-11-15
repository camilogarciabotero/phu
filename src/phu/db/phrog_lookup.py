"""Find PHROGs by name: aliases, functional categories and annotation text.

``resolve`` turns queries into PHROG identifiers that ``phu screen`` accepts:

- ``phrog_216``          one PHROG
- ``alias:mcp``          a curated group (src/phu/data/phrog_aliases.json)
- ``category:tail``      a PHROG functional category, matched forgivingly
- ``annot:portal protein`` an exact annotation, ignoring letter case

``search`` is for exploring: annotation text that contains the words given.
"""

from __future__ import annotations

import difflib
import json
import re
from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from importlib import resources

from phu.db.phrog import PhrogAnnotation, is_phrog_id, normalize_phrog_id

UNKNOWN_CATEGORY = "unknown function"
QUERY_KINDS = ("alias", "category", "annot")


@dataclass(frozen=True)
class Alias:
    name: str
    description: str
    annotations: tuple[str, ...]


@dataclass(frozen=True)
class AliasTable:
    aliases: dict[str, Alias]
    categories: dict[str, str]  # shortcut -> category


@dataclass(frozen=True)
class Resolved:
    """One PHROG selected by a query."""

    phrog: PhrogAnnotation
    query: str  # the query that selected it, e.g. "alias:mcp"


def load_aliases() -> AliasTable:
    data = json.loads(
        resources.files("phu").joinpath("data", "phrog_aliases.json").read_text("utf-8")
    )
    aliases = {
        name.lower(): Alias(
            name.lower(), entry["description"], tuple(entry["annotations"])
        )
        for name, entry in data["aliases"].items()
    }
    return AliasTable(aliases, {k.lower(): v for k, v in data["categories"].items()})


def _normal(text: str) -> str:
    return " ".join(re.sub(r"[^0-9a-z]+", " ", text.lower()).split())


def category_counts(phrogs: Iterable[PhrogAnnotation]) -> Counter[str]:
    return Counter(p.category for p in phrogs)


def match_category(
    query: str, categories: Iterable[str], shortcuts: Mapping[str, str] | None = None
) -> list[str]:
    """Categories meant by ``query``, most specific rule first.

    1. A shortcut (``amg``) or the exact category, ignoring case and punctuation.
    2. Every query word starts a word of the category (``head`` or ``pack`` ->
       "head and packaging", ``dna`` -> "DNA, RNA and nucleotide metabolism").
    3. Close spelling (``tails``, ``lysys``), by difflib.

    Returns the matches of the first rule that finds any: one for a clear
    query, several when the query is ambiguous, none when nothing is close.
    """
    names = sorted(set(categories))
    wanted = _normal(query)
    if not wanted:
        return []
    shortcut = (shortcuts or {}).get(wanted) or (shortcuts or {}).get(
        query.strip().lower()
    )
    if shortcut in names:
        return [shortcut]
    exact = [name for name in names if _normal(name) == wanted]
    if exact:
        return exact

    words = wanted.split()
    prefix = [
        name
        for name in names
        if all(any(w.startswith(q) for w in _normal(name).split()) for q in words)
    ]
    if prefix:
        return prefix

    by_normal = {_normal(name): name for name in names}
    close = difflib.get_close_matches(wanted, list(by_normal), n=3, cutoff=0.6)
    if not close:
        # A misspelt word against each category's words ("lysys" -> "lysis").
        scored = []
        for normal, name in by_normal.items():
            score = max(
                difflib.SequenceMatcher(None, wanted, word).ratio()
                for word in normal.split()
            )
            if score >= 0.75:
                scored.append((score, name))
        best = max((score for score, _ in scored), default=0)
        return [name for score, name in scored if score == best]
    return [by_normal[c] for c in close[:1]]


class QueryError(ValueError):
    """A query that selects nothing, or more than intended."""


def _split(query: str) -> tuple[str, str]:
    if is_phrog_id(query):
        return "phrog", normalize_phrog_id(query)
    kind, sep, value = query.partition(":")
    kind = kind.strip().lower()
    if not sep or kind not in QUERY_KINDS or not value.strip():
        raise QueryError(
            f"Unrecognised query {query!r}: use phrog_<n>, alias:<name>, "
            "category:<name> or annot:<annotation>"
        )
    return kind, value.strip()


def resolve(
    queries: Iterable[str],
    phrogs: Mapping[str, PhrogAnnotation],
    aliases: AliasTable,
) -> tuple[list[Resolved], list[str]]:
    """PHROGs selected by ``queries``, in query then identifier order.

    ``phrogs`` holds the PHROGs that can be searched (those with a model).
    Returns the selection, without repeats, and notes for the user (how a
    fuzzy category was read). Raises QueryError for a query that selects
    nothing or names an ambiguous category.
    """
    selected: dict[str, Resolved] = {}
    notes: list[str] = []
    categories = {p.category for p in phrogs.values()}

    def number(phrog_id: str) -> int:
        return int(phrog_id.split("_", 1)[1])

    for query in queries:
        kind, value = _split(query)
        if kind == "phrog":
            if value not in phrogs:
                raise QueryError(f"{value} has no model in the PHROG database")
            matches = [phrogs[value]]
        elif kind == "alias":
            alias = aliases.aliases.get(value.lower())
            if alias is None:
                close = difflib.get_close_matches(
                    value.lower(), list(aliases.aliases), n=3
                )
                hint = f"; did you mean {', '.join(close)}?" if close else ""
                raise QueryError(
                    f"Unknown alias {value!r}{hint} (see `phu phrogs aliases`)"
                )
            wanted = {a.lower() for a in alias.annotations}
            matches = [p for p in phrogs.values() if p.annotation.lower() in wanted]
        elif kind == "category":
            found = match_category(value, categories, aliases.categories)
            if not found:
                raise QueryError(
                    f"No PHROG category matches {value!r}; categories: "
                    + "; ".join(sorted(categories))
                )
            if len(found) > 1:
                raise QueryError(
                    f"category:{value} is ambiguous: "
                    + "; ".join(found)
                    + ". Be more specific."
                )
            if _normal(found[0]) != _normal(value):
                notes.append(f"category:{value} read as {found[0]!r}")
            matches = [p for p in phrogs.values() if p.category == found[0]]
        else:  # annot
            matches = [
                p for p in phrogs.values() if p.annotation.lower() == value.lower()
            ]
        if not matches:
            raise QueryError(f"{query} selects no PHROG with a model")
        for phrog in sorted(matches, key=lambda p: number(p.phrog_id)):
            selected.setdefault(phrog.phrog_id, Resolved(phrog, query))
    return list(selected.values()), notes


def search(
    text: str,
    phrogs: Mapping[str, PhrogAnnotation],
    *,
    category: str | None = None,
) -> list[PhrogAnnotation]:
    """PHROGs whose annotation contains every word of ``text`` (any order)."""
    words = _normal(text).split()
    hits = [
        p
        for p in phrogs.values()
        if p.annotation
        and all(word in _normal(p.annotation) for word in words)
        and (category is None or p.category == category)
    ]
    return sorted(
        hits, key=lambda p: (p.annotation.lower(), int(p.phrog_id.split("_")[1]))
    )


def check_aliases(
    aliases: AliasTable, phrogs: Mapping[str, PhrogAnnotation]
) -> dict[str, list[str]]:
    """Annotation strings of each alias that match no PHROG, by alias."""
    present = {p.annotation.lower() for p in phrogs.values()}
    categories = {p.category for p in phrogs.values()}
    problems = {
        alias.name: [a for a in alias.annotations if a.lower() not in present]
        for alias in aliases.aliases.values()
    }
    problems.update(
        {
            f"category shortcut {name}": [target]
            for name, target in aliases.categories.items()
            if target not in categories
        }
    )
    return {name: missing for name, missing in problems.items() if missing}
