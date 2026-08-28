from __future__ import annotations

import csv
import gzip
import json
import logging
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from multiprocessing.pool import ThreadPool
from pathlib import Path
from typing import Optional

import pyhmmer
import pyhmmer.easel
import pyhmmer.plan7
import typer

from ._click import run_click_task
from ._exec import _executable
from .dbcan_db import (
    ensure_dbcan_database,
    extract_dbcan_models,
    get_dbcan_database_status,
    get_dbcan_model_inventory,
    get_dbcan_pul,
    get_dbcan_pul_rules,
    is_dbcan_id,
    is_dbcan_pul_id,
    normalize_dbcan_id,
    normalize_dbcan_pul_id,
)
from .gene_prediction_core import (
    PredictionInputs,
    get_or_predict_proteins,
    predict_genes_pyrodigal,
    write_predicted_proteins_fasta,
    write_prediction_metadata,
)
from .kofam_db import (
    KOFamMetadata,
    ensure_kofam_database,
    extract_kofam_models,
    get_kofam_metadata_map,
    is_kofam_id,
    normalize_kofam_id,
)
from .pfam_db import (
    ensure_pfam_database,
    extract_pfam_models,
    is_pfam_id,
    normalize_pfam_id,
)

logger = logging.getLogger(__name__)

app = typer.Typer(
    help="Screen contigs for a protein family using pyHMMER on predicted CDS."
)


@dataclass(frozen=True)
class ScreenQuery:
    original_inputs: tuple[str, ...]
    normalized_models: tuple[str, ...]
    pul_id: str | None = None
    pul_substrate: str | None = None
    pul_raw_rule: str | None = None
    matching_rule: str = "any"
    all_puls: bool = False
    total_pul_count: int = 0
    resolvable_pul_count: int = 0
    skipped_puls: tuple[dict[str, object], ...] = ()
    all_cazymes: bool = False
    total_hmm_profiles: int = 0
    excluded_ancillary_profiles: tuple[str, ...] = ()
    database_version: str | None = None
    database_digest: str | None = None


@dataclass(frozen=True)
class PULMatch:
    contig_id: str
    pul_id: str
    substrate: str
    required_families: tuple[str, ...]
    matched_families: tuple[str, ...]


CAZYME_CLASSES = ("AA", "CBM", "CE", "GH", "GT", "PL")


def cazyme_class(family: str) -> str:
    return next((item for item in CAZYME_CLASSES if family.startswith(item)), "")


def _best_hit_key(hit: Hit) -> tuple[float, float, float, str]:
    return (-hit.bitscore, hit.evalue, -(hit.hmm_coverage or 0.0), hit.prot_id)


def _dbcan_hashes() -> dict[str, str]:
    manifest_path = Path(str(get_dbcan_database_status()["manifest_path"]))
    if not manifest_path.exists():
        return {}
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    return {
        name: item["sha256"]
        for name, item in manifest.get("files", {}).items()
        if isinstance(item, dict) and item.get("sha256")
    }


def _write_tsv(
    path: Path,
    fieldnames: list[str],
    rows: Iterable[dict[str, object]],
    *,
    compressed: bool = False,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if compressed else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_cazyme_evidence(
    hits: Iterable[Hit], contig_order: Iterable[str], output_path: Path
) -> None:
    order = {contig: index for index, contig in enumerate(contig_order)}
    rows = sorted(
        hits,
        key=lambda hit: (
            order.get(hit.contig, len(order)),
            hit.model,
            hit.prot_id,
            hit.evalue,
            -hit.bitscore,
        ),
    )
    fields = [
        "contig_id", "protein_id", "cazyme_family", "cazyme_class", "bitscore",
        "evalue", "domain_i_evalue", "hmm_coverage", "hmm_from", "hmm_to",
        "target_from", "target_to",
    ]
    values = []
    for hit in rows:
        values.append({
            "contig_id": hit.contig,
            "protein_id": hit.prot_id,
            "cazyme_family": hit.model,
            "cazyme_class": cazyme_class(hit.model),
            "bitscore": f"{hit.bitscore:.6f}",
            "evalue": f"{hit.evalue:.6g}",
            "domain_i_evalue": "" if hit.domain_i_evalue is None else f"{hit.domain_i_evalue:.6g}",
            "hmm_coverage": "" if hit.hmm_coverage is None else f"{hit.hmm_coverage:.6f}",
            "hmm_from": "" if hit.hmm_from is None else hit.hmm_from,
            "hmm_to": "" if hit.hmm_to is None else hit.hmm_to,
            "target_from": "" if hit.target_from is None else hit.target_from,
            "target_to": "" if hit.target_to is None else hit.target_to,
        })
    _write_tsv(output_path, fields, values, compressed=True)


def write_cazyme_outputs(
    evidence_hits: Iterable[Hit],
    selected_hits: Iterable[Hit],
    contig_order: Iterable[str],
    output_dir: Path,
) -> dict[str, int]:
    evidence = list(evidence_hits)
    selected = list(selected_hits)
    order = {contig: index for index, contig in enumerate(contig_order)}
    evidence.sort(key=lambda hit: (order.get(hit.contig, len(order)), hit.model, hit.prot_id, hit.evalue, -hit.bitscore))
    evidence_fields = [
        "contig_id", "protein_id", "cazyme_family", "cazyme_class", "bitscore",
        "evalue", "domain_i_evalue", "hmm_coverage", "hmm_from", "hmm_to",
        "target_from", "target_to",
    ]
    evidence_rows = []
    for hit in evidence:
        evidence_rows.append({
            "contig_id": hit.contig,
            "protein_id": hit.prot_id,
            "cazyme_family": hit.model,
            "cazyme_class": cazyme_class(hit.model),
            "bitscore": f"{hit.bitscore:.6f}",
            "evalue": f"{hit.evalue:.6g}",
            "domain_i_evalue": "" if hit.domain_i_evalue is None else f"{hit.domain_i_evalue:.6g}",
            "hmm_coverage": "" if hit.hmm_coverage is None else f"{hit.hmm_coverage:.6f}",
            "hmm_from": hit.hmm_from if hit.hmm_from is not None else "",
            "hmm_to": hit.hmm_to if hit.hmm_to is not None else "",
            "target_from": hit.target_from if hit.target_from is not None else "",
            "target_to": hit.target_to if hit.target_to is not None else "",
        })
    _write_tsv(output_dir / "evidence" / "cazyme_hits.tsv.gz", evidence_fields, evidence_rows, compressed=True)

    grouped: dict[tuple[str, str], list[Hit]] = defaultdict(list)
    for hit in selected:
        grouped[(hit.contig, hit.model)].append(hit)
    match_rows = []
    for (contig, family), hits in grouped.items():
        hits.sort(key=_best_hit_key)
        best = hits[0]
        match_rows.append({
            "contig_id": contig,
            "cazyme_family": family,
            "cazyme_class": cazyme_class(family),
            "protein_count": len({hit.prot_id for hit in hits}),
            "hit_count": len(hits),
            "protein_ids": ";".join(sorted({hit.prot_id for hit in hits})),
            "best_protein_id": best.prot_id,
            "best_bitscore": f"{best.bitscore:.6f}",
            "best_evalue": f"{best.evalue:.6g}",
            "best_domain_i_evalue": "" if best.domain_i_evalue is None else f"{best.domain_i_evalue:.6g}",
            "best_hmm_coverage": "" if best.hmm_coverage is None else f"{best.hmm_coverage:.6f}",
        })
    match_fields = [
        "contig_id", "cazyme_family", "cazyme_class", "protein_count", "hit_count",
        "protein_ids", "best_protein_id", "best_bitscore", "best_evalue",
        "best_domain_i_evalue", "best_hmm_coverage",
    ]
    match_rows.sort(key=lambda row: (order.get(str(row["contig_id"]), len(order)), str(row["cazyme_family"])))
    _write_tsv(output_dir / "cazyme_matches.tsv.gz", match_fields, match_rows, compressed=True)

    family_groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in match_rows:
        family_groups[str(row["cazyme_family"])].append(row)
    summary_rows = []
    for family in sorted(family_groups):
        rows = family_groups[family]
        best = min(rows, key=lambda row: (-float(row["best_bitscore"]), float(row["best_evalue"]), str(row["best_protein_id"])))
        summary_rows.append({
            "cazyme_family": family,
            "cazyme_class": cazyme_class(family),
            "contig_count": len({str(row["contig_id"]) for row in rows}),
            "protein_count": len({protein for row in rows for protein in str(row["protein_ids"]).split(";")}),
            "hit_count": sum(int(row["hit_count"]) for row in rows),
            "best_bitscore": best["best_bitscore"],
            "best_evalue": best["best_evalue"],
        })
    _write_tsv(output_dir / "cazyme_summary.tsv", [
        "cazyme_family", "cazyme_class", "contig_count", "protein_count", "hit_count", "best_bitscore", "best_evalue"
    ], summary_rows)

    class_rows = []
    for klass in CAZYME_CLASSES:
        rows = [row for row in match_rows if row["cazyme_class"] == klass]
        if not rows:
            continue
        class_rows.append({
            "cazyme_class": klass,
            "matched_family_count": len({str(row["cazyme_family"]) for row in rows}),
            "contig_count": len({str(row["contig_id"]) for row in rows}),
            "protein_count": len({protein for row in rows for protein in str(row["protein_ids"]).split(";")}),
            "hit_count": sum(int(row["hit_count"]) for row in rows),
        })
    _write_tsv(output_dir / "cazyme_class_summary.tsv", [
        "cazyme_class", "matched_family_count", "contig_count", "protein_count", "hit_count"
    ], class_rows)
    return {
        "qualifying_hits": len(evidence),
        "matched_proteins": len({hit.prot_id for hit in selected}),
        "matched_families": len({hit.model for hit in selected}),
        "matched_classes": len({cazyme_class(hit.model) for hit in selected}),
    }


def evaluate_pul_signatures(
    hits: Iterable[Hit], pul_rules: Iterable[object]
) -> tuple[list[PULMatch], list[str]]:
    """Evaluate independent PUL signatures against threshold-passing model hits."""
    observed: dict[str, set[str]] = defaultdict(set)
    contig_order: list[str] = []
    for hit in hits:
        if hit.contig not in observed:
            contig_order.append(hit.contig)
        observed[hit.contig].add(hit.model)

    matches: list[PULMatch] = []
    for contig in contig_order:
        families = observed[contig]
        for rule in sorted(pul_rules, key=lambda item: item.pul_id):
            required = tuple(rule.families)
            if set(required).issubset(families):
                matches.append(
                    PULMatch(
                        contig_id=contig,
                        pul_id=rule.pul_id,
                        substrate=rule.substrate,
                        required_families=required,
                        matched_families=tuple(
                            family for family in required if family in families
                        ),
                    )
                )
    return matches, list(dict.fromkeys(match.contig_id for match in matches))


def write_pul_outputs(
    matches: Iterable[PULMatch],
    qualifying_hits: Iterable[Hit],
    contig_order: Iterable[str],
    output_dir: Path,
) -> dict[str, int]:
    matches = list(matches)
    hits = list(qualifying_hits)
    order = {contig: index for index, contig in enumerate(contig_order)}
    hits_by_key: dict[tuple[str, str], list[Hit]] = defaultdict(list)
    for hit in hits:
        hits_by_key[(hit.contig, hit.model)].append(hit)

    def best(rows: list[Hit]) -> Hit:
        return sorted(rows, key=_best_hit_key)[0]

    support_rows = []
    match_rows = []
    for match in sorted(matches, key=lambda item: (order.get(item.contig_id, len(order)), item.pul_id)):
        supporting = [
            hit
            for family in match.required_families
            for hit in hits_by_key.get((match.contig_id, family), [])
        ]
        for family in match.required_families:
            family_hits = hits_by_key.get((match.contig_id, family), [])
            if not family_hits:
                continue
            selected = best(family_hits)
            support_rows.append({
                "contig_id": match.contig_id,
                "pul_id": match.pul_id,
                "cazyme_family": family,
                "protein_count": len({hit.prot_id for hit in family_hits}),
                "hit_count": len(family_hits),
                "protein_ids": ";".join(sorted({hit.prot_id for hit in family_hits})),
                "best_protein_id": selected.prot_id,
                "best_bitscore": f"{selected.bitscore:.6f}",
                "best_evalue": f"{selected.evalue:.6g}",
                "best_domain_i_evalue": "" if selected.domain_i_evalue is None else f"{selected.domain_i_evalue:.6g}",
                "best_hmm_coverage": "" if selected.hmm_coverage is None else f"{selected.hmm_coverage:.6f}",
                "evidence_start": "",
                "evidence_end": "",
            })
        match_rows.append({
            "contig_id": match.contig_id,
            "contig_length_bp": "",
            "pul_id": match.pul_id,
            "reference_substrate": match.substrate,
            "required_family_count": len(match.required_families),
            "required_families": ";".join(match.required_families),
            "matched_family_count": len(match.matched_families),
            "matched_families": ";".join(match.matched_families),
            "matched_protein_count": len({hit.prot_id for hit in supporting}),
            "hmm_hit_count": len(supporting),
            "matched_protein_ids": ";".join(sorted({hit.prot_id for hit in supporting})),
            "evidence_start": "",
            "evidence_end": "",
            "evidence_span_bp": "",
        })

    _write_tsv(output_dir / "pul_matches.tsv.gz", [
        "contig_id", "contig_length_bp", "pul_id", "reference_substrate",
        "required_family_count", "required_families", "matched_family_count",
        "matched_families", "matched_protein_count", "hmm_hit_count",
        "matched_protein_ids", "evidence_start", "evidence_end", "evidence_span_bp",
    ], match_rows, compressed=True)
    _write_tsv(output_dir / "pul_family_support.tsv.gz", [
        "contig_id", "pul_id", "cazyme_family", "protein_count", "hit_count",
        "protein_ids", "best_protein_id", "best_bitscore", "best_evalue",
        "best_domain_i_evalue", "best_hmm_coverage", "evidence_start", "evidence_end",
    ], support_rows, compressed=True)

    pul_groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in match_rows:
        pul_groups[str(row["pul_id"])].append(row)
    pul_summary = []
    for pul_id, rows in sorted(pul_groups.items()):
        pul_summary.append({
            "pul_id": pul_id,
            "reference_substrate": rows[0]["reference_substrate"],
            "required_family_count": rows[0]["required_family_count"],
            "required_families": rows[0]["required_families"],
            "matched_contig_count": len({row["contig_id"] for row in rows}),
            "matched_protein_count": len({protein for row in rows for protein in str(row["matched_protein_ids"]).split(";") if protein}),
            "hmm_hit_count": sum(int(row["hmm_hit_count"]) for row in rows),
        })
    _write_tsv(output_dir / "pul_summary.tsv", [
        "pul_id", "reference_substrate", "required_family_count", "required_families",
        "matched_contig_count", "matched_protein_count", "hmm_hit_count",
    ], pul_summary)
    substrate_groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in pul_summary:
        substrate_groups[str(row["reference_substrate"])].append(row)
    substrate_summary = []
    for substrate, rows in sorted(substrate_groups.items()):
        pul_ids = {str(row["pul_id"]) for row in rows}
        contigs = {str(match["contig_id"]) for match in match_rows if str(match["pul_id"]) in pul_ids}
        proteins = {protein for match in match_rows if str(match["pul_id"]) in pul_ids for protein in str(match["matched_protein_ids"]).split(";") if protein}
        substrate_summary.append({
            "reference_substrate": substrate,
            "supporting_pul_count": len(pul_ids),
            "matched_contig_count": len(contigs),
            "matched_protein_count": len(proteins),
        })
    _write_tsv(output_dir / "substrate_summary.tsv", [
        "reference_substrate", "supporting_pul_count", "matched_contig_count", "matched_protein_count",
    ], substrate_summary)
    return {"complete_pul_match_count": len(matches)}


# ---------- Utilities ----------


def _cmd_exists(exe: str) -> bool:
    return shutil.which(exe) is not None


def _resolve_all_puls_inputs(
    outdir: Path,
) -> tuple[list[Path], ScreenQuery, tuple[object, ...]]:
    """Resolve every resolvable PUL into one ordered union of family models."""
    run_click_task("Preparing dbCAN database", ensure_dbcan_database)
    rules = get_dbcan_pul_rules()
    resolvable = tuple(rule for rule in rules.values() if rule.resolved)
    skipped = tuple(
        {
            "pul_id": rule.pul_id,
            "reasons": list(rule.unresolved_tokens) or ["empty rule"],
        }
        for rule in rules.values()
        if not rule.resolved
    )
    if not resolvable:
        raise ValueError("No resolvable dbCAN PUL rules remain; refresh the dbCAN database")

    families = list(
        dict.fromkeys(family for rule in resolvable for family in rule.families)
    )
    dbcan_hmms, missing = run_click_task(
        "Resolving dbCAN PUL families",
        extract_dbcan_models,
        requested_ids=families,
        output_dir=outdir / "resolved_dbcan_hmms",
    )
    if missing:
        raise FileNotFoundError(
            "dbCAN family model(s) missing for --all-puls: " + ", ".join(missing)
        )
    query = ScreenQuery(
        original_inputs=("--all-puls",),
        normalized_models=tuple(families),
        matching_rule="any_pul_all_families",
        all_puls=True,
        total_pul_count=len(rules),
        resolvable_pul_count=len(resolvable),
        skipped_puls=skipped,
    )
    if skipped:
        print(f"Warning: skipped {len(skipped)} unresolved dbCAN PUL rule(s)")
    print(
        f"Resolved {len(dbcan_hmms)} unique dbCAN family model(s) for "
        f"{len(resolvable)} PUL signature(s)"
    )
    return dbcan_hmms, query, resolvable


def _resolve_all_cazymes_inputs(
    outdir: Path,
) -> tuple[list[Path], ScreenQuery]:
    """Resolve every indexed canonical dbCAN family exactly once."""
    database_manifest = run_click_task("Preparing dbCAN database", ensure_dbcan_database)
    canonical, excluded = get_dbcan_model_inventory()
    if not canonical:
        raise ValueError("No canonical dbCAN CAZyme profiles are available")
    dbcan_hmms, missing = run_click_task(
        "Resolving all dbCAN CAZyme families",
        extract_dbcan_models,
        requested_ids=list(canonical),
        output_dir=outdir / "resolved_dbcan_hmms",
    )
    if missing:
        raise FileNotFoundError(
            "dbCAN profiles missing from the installed index: " + ", ".join(missing)
        )
    query = ScreenQuery(
        original_inputs=("--all-cazymes",),
        normalized_models=canonical,
        matching_rule="any",
        all_cazymes=True,
        total_hmm_profiles=len(canonical) + len(excluded),
        excluded_ancillary_profiles=excluded,
        database_version=database_manifest.get("release"),
        database_digest=database_manifest.get("files", {}).get("hmm", {}).get("sha256"),
    )
    print(
        f"Resolved {len(dbcan_hmms)} canonical dbCAN CAZyme profile(s); "
        f"excluded {len(excluded)} ancillary profile(s)"
    )
    return dbcan_hmms, query


def _resolve_hmm_inputs(
    hmm_inputs: list[Path],
    outdir: Path,
    *,
    return_query: bool = False,
) -> tuple[list[Path], dict[str, KOFamMetadata]] | tuple[list[Path], dict[str, KOFamMetadata], ScreenQuery]:
    """
    Resolve positional inputs into concrete HMM file paths.

    Input tokens may be local HMM paths, PFAM/KOfam identifiers, dbCAN families, or PUL IDs.
    """
    local_hmms: list[Path] = []
    pfam_ids: list[str] = []
    ko_ids: list[str] = []
    dbcan_ids: list[str] = []
    pul_ids: list[str] = []

    for token_path in hmm_inputs:
        token = str(token_path)

        if token_path.exists():
            local_hmms.append(token_path)
            continue

        if is_pfam_id(token):
            pfam_ids.append(normalize_pfam_id(token))
            continue

        if is_kofam_id(token):
            ko_ids.append(normalize_kofam_id(token))
            continue

        if is_dbcan_pul_id(token):
            pul_ids.append(normalize_dbcan_pul_id(token))
            continue

        if is_dbcan_id(token):
            dbcan_ids.append(normalize_dbcan_id(token))
            continue

        raise FileNotFoundError(f"HMM file not found: {token_path}")

    dedup_pfam_ids = list(dict.fromkeys(pfam_ids))
    dedup_ko_ids = list(dict.fromkeys(ko_ids))
    dedup_dbcan_ids = list(dict.fromkeys(dbcan_ids))
    dedup_pul_ids = list(dict.fromkeys(pul_ids))

    if dedup_pul_ids and (
        len(dedup_pul_ids) != 1
        or local_hmms
        or dedup_pfam_ids
        or dedup_ko_ids
        or dedup_dbcan_ids
    ):
        raise ValueError(
            "A dbCAN PUL query must be the only database input and exactly one PUL is supported"
        )

    resolved: list[Path] = list(local_hmms)

    if dedup_pfam_ids:
        pfam_meta = run_click_task("Preparing PFAM database", ensure_pfam_database)
        pfam_hmm_db_path = Path(pfam_meta["hmm_path"])

        pfam_dir = outdir / "resolved_pfam_hmms"
        pfam_hmms, missing = run_click_task(
            "Resolving PFAM accessions",
            extract_pfam_models,
            hmm_db_path=pfam_hmm_db_path,
            requested_ids=dedup_pfam_ids,
            output_dir=pfam_dir,
        )

        if missing:
            missing_msg = ", ".join(missing)
            raise FileNotFoundError(
                f"PFAM accession(s) not found in local database: {missing_msg}. "
                "Run again later to refresh the DB source if needed."
            )

        print(
            f"Resolved {len(pfam_hmms)} PFAM model(s) from local database: "
            f"{', '.join(dedup_pfam_ids)}"
        )
        resolved.extend(pfam_hmms)

    ko_metadata: dict[str, KOFamMetadata] = {}
    if dedup_ko_ids:
        run_click_task("Preparing KOFam database", ensure_kofam_database)

        kofam_dir = outdir / "resolved_kofam_hmms"
        kofam_hmms, missing = run_click_task(
            "Resolving KO identifiers",
            extract_kofam_models,
            requested_ids=dedup_ko_ids,
            output_dir=kofam_dir,
        )

        if missing:
            missing_msg = ", ".join(missing)
            raise FileNotFoundError(
                f"KO identifier(s) not found in local KOFam database: {missing_msg}. "
                "Run `phu dbs refresh kofam` and try again."
            )

        ko_metadata = get_kofam_metadata_map(dedup_ko_ids)
        print(
            f"Resolved {len(kofam_hmms)} KOFam model(s) from local database: "
            f"{', '.join(dedup_ko_ids)}"
        )
        resolved.extend(kofam_hmms)

    if dedup_pul_ids:
        run_click_task("Preparing dbCAN database", ensure_dbcan_database)
        pul = get_dbcan_pul(dedup_pul_ids[0])
        if pul is None:
            raise FileNotFoundError(
                f"PUL identifier not found in local dbCAN database: {dedup_pul_ids[0]}"
            )
        if not pul.resolved:
            unresolved = ", ".join(pul.unresolved_tokens) or "empty rule"
            raise ValueError(
                f"PUL {pul.pul_id} has unresolved dbCAN rule token(s): {unresolved}"
            )
        dedup_dbcan_ids = list(pul.families)
        dbcan_dir = outdir / "resolved_dbcan_hmms"
        dbcan_hmms, missing = run_click_task(
            "Resolving dbCAN PUL families",
            extract_dbcan_models,
            requested_ids=dedup_dbcan_ids,
            output_dir=dbcan_dir,
        )
        if missing:
            raise FileNotFoundError(
                "dbCAN family model(s) not found for PUL "
                f"{pul.pul_id}: {', '.join(missing)}"
            )
        query = ScreenQuery(
            original_inputs=tuple(str(value) for value in hmm_inputs),
            normalized_models=tuple(dedup_dbcan_ids),
            pul_id=pul.pul_id,
            pul_substrate=pul.substrate,
            pul_raw_rule=pul.raw_rule,
            matching_rule="all",
        )
        print(
            f"Expanded {pul.pul_id} to dbCAN CAZyme signature: "
            f"{', '.join(dedup_dbcan_ids)}"
        )
        resolved.extend(dbcan_hmms)
    elif dedup_dbcan_ids:
        run_click_task("Preparing dbCAN database", ensure_dbcan_database)
        dbcan_dir = outdir / "resolved_dbcan_hmms"
        dbcan_hmms, missing = run_click_task(
            "Resolving dbCAN families",
            extract_dbcan_models,
            requested_ids=dedup_dbcan_ids,
            output_dir=dbcan_dir,
        )
        if missing:
            raise FileNotFoundError(
                "dbCAN family model(s) not found in local database: "
                f"{', '.join(missing)}. Run `phu dbs refresh dbcan` and try again."
            )
        query = ScreenQuery(
            original_inputs=tuple(str(value) for value in hmm_inputs),
            normalized_models=tuple(dedup_dbcan_ids),
        )
        print(
            f"Resolved {len(dbcan_hmms)} dbCAN model(s) from local database: "
            f"{', '.join(dedup_dbcan_ids)}"
        )
        resolved.extend(dbcan_hmms)

    if "query" not in locals():
        query = ScreenQuery(
            original_inputs=tuple(str(value) for value in hmm_inputs),
            normalized_models=(*dedup_pfam_ids, *dedup_ko_ids, *dedup_dbcan_ids),
        )
    if return_query:
        return resolved, ko_metadata, query
    return resolved, ko_metadata


@dataclass
class ScreenConfig:
    """Configuration for screening contigs for protein families."""

    input_contigs: Path
    hmms: list[Path]  # Changed from hmm: Path to support multiple HMMs
    all_puls: bool = False
    all_cazymes: bool = False
    outdir: Path = Path("phu-screen")
    mode: str = "meta"  # pyrodigal mode: meta|single
    threads: int = 1
    min_bitscore: Optional[float] = None
    max_evalue: Optional[float] = 1e-5
    cut_ga: bool = True
    use_kofam_thresholds: bool = True
    top_per_contig: int = 1
    min_protein_len_aa: int = 30
    translation_table: int = 11
    keep_proteins: bool = False
    keep_domtbl: bool = True
    combine_mode: str = "any"  # New: how to combine hits from multiple HMMs
    min_hmm_hits: int = 1  # New: minimum number of HMMs that must hit a contig
    save_target_proteins: bool = False  # New: save matched proteins per HMM
    save_target_hmms: bool = False  # New: build and save HMMs from target proteins
    hmm_mode: str = (
        "pure"  # New: "pure" for single models, "mixed" for pressed/concatenated HMMs
    )

    def __post_init__(self):
        """Validate configuration parameters."""
        if self.threads < 0:
            raise ValueError("threads must be >= 0")
        if self.threads == 0:
            # HMMER interprets --cpu 0 as "turn off multithreading"
            # This is valid, so we allow it
            pass
        if not self.hmms and not self.all_puls and not self.all_cazymes:
            raise ValueError("At least one HMM file must be provided")
        if self.all_puls and self.all_cazymes:
            raise ValueError("all_puls and all_cazymes cannot both be enabled")
        if self.all_cazymes and self.combine_mode != "any":
            raise ValueError("--combine-mode is not applicable with --all-cazymes")
        if self.combine_mode not in {"any", "all", "threshold"}:
            raise ValueError("combine_mode must be 'any', 'all', or 'threshold'")
        if self.hmm_mode not in {"pure", "mixed"}:
            raise ValueError("hmm_mode must be 'pure' or 'mixed'")
        if self.save_target_hmms and not self.save_target_proteins:
            raise ValueError(
                "save_target_hmms requires save_target_proteins to be True"
            )
        if self.min_protein_len_aa < 1:
            raise ValueError("min_protein_len_aa must be >= 1")

    def plan(self) -> "ScreenPlan":
        """Create execution plan from configuration."""
        if self.mode not in {"meta", "single"}:
            raise ValueError("mode must be 'meta' or 'single'")

        effective_threads = self.threads

        # Create domtbl paths for each HMM
        domtbl_paths = {}
        for hmm in self.hmms:
            hmm_name = hmm.stem  # filename without extension
            domtbl_paths[hmm_name] = self.outdir / f"hits_{hmm_name}.domtblout"

        return ScreenPlan(
            hmmer_bin="",
            seqkit_bin="",
            input_contigs=self.input_contigs,
            hmms=self.hmms,
            outdir=self.outdir,
            mode=self.mode,
            threads=effective_threads,
            min_bitscore=self.min_bitscore,
            max_evalue=self.max_evalue,
            cut_ga=self.cut_ga,
            use_kofam_thresholds=self.use_kofam_thresholds,
            top_per_contig=self.top_per_contig,
            min_protein_len_aa=self.min_protein_len_aa,
            translation_table=self.translation_table,
            keep_proteins=self.keep_proteins,
            keep_domtbl=self.keep_domtbl,
            combine_mode=self.combine_mode,
            min_hmm_hits=self.min_hmm_hits,
            save_target_proteins=self.save_target_proteins,
            save_target_hmms=self.save_target_hmms,
            hmm_mode=self.hmm_mode,
            proteins_fa=self.outdir / "proteins.faa",
            domtbl_paths=domtbl_paths,
            kept_ids=self.outdir / "kept_contigs.txt",
            out_contigs=self.outdir / "screened_contigs.fasta",
        )


@dataclass
class ScreenPlan:
    """Execution plan for screening operation."""

    hmmer_bin: str
    seqkit_bin: str
    input_contigs: Path
    hmms: list[Path]  # Changed from hmm: Path
    outdir: Path
    mode: str
    threads: int
    min_bitscore: Optional[float]
    max_evalue: Optional[float]
    cut_ga: bool
    use_kofam_thresholds: bool
    top_per_contig: int
    min_protein_len_aa: int
    translation_table: int
    keep_proteins: bool
    keep_domtbl: bool
    combine_mode: str
    min_hmm_hits: int
    save_target_proteins: bool  # New
    save_target_hmms: bool  # New
    hmm_mode: str  # New
    proteins_fa: Path
    domtbl_paths: dict[str, Path]  # Changed from domtbl: Path
    kept_ids: Path
    out_contigs: Path


def _binaries() -> str:
    """
    Discover required binaries for screening.
    Only seqkit is needed since we use pyHMMER instead of HMMER binary.
    """
    seqkit = _executable(["seqkit"])
    return seqkit


def _read_fasta(fp: Path) -> Iterable[tuple[str, str]]:
    """
    Read FASTA sequences with two strategies:
    1) Python parser (stdlib; robust for compressed/edge-case inputs)
    2) pyhmmer/easel fallback (fast native parser when available)

    This keeps a dependency-light default while preserving a pyhmmer-backed
    alternative in environments where easel parsing is preferred.
    """
    try:
        yield from _read_fasta_python(fp)
        return
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        # Fallback to easel parser when Python parsing fails unexpectedly.
        logger.debug("Python FASTA parsing failed for %s; trying easel: %s", fp, exc)
        yield from _read_fasta_easel(fp)


def _read_fasta_easel(fp: Path) -> Iterable[tuple[str, str]]:
    """Read FASTA records using pyhmmer.easel.SequenceFile."""
    with pyhmmer.easel.SequenceFile(str(fp)) as seq_file:
        for seq in seq_file:
            seq_id = seq.name.decode() if isinstance(seq.name, bytes) else seq.name
            seq_seq = seq.sequence
            if isinstance(seq_seq, (bytes, bytearray, memoryview)):
                seq_seq = bytes(seq_seq).decode()
            elif not isinstance(seq_seq, str):
                seq_seq = str(seq_seq)
            yield seq_id, seq_seq


def _read_fasta_python(fp: Path) -> Iterable[tuple[str, str]]:
    """Read FASTA records using Python text IO (supports .gz)."""
    opener = gzip.open if fp.suffix == ".gz" else open
    with opener(fp, "rt") as handle:
        seq_id: Optional[str] = None
        seq_chunks: list[str] = []

        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(seq_chunks)
                # FASTA header id is first token after '>'
                seq_id = line[1:].split(None, 1)[0]
                seq_chunks = []
                continue

            if seq_id is None:
                raise ValueError(
                    f"Invalid FASTA format in {fp}: sequence before header"
                )

            seq_chunks.append(line)

        if seq_id is not None:
            yield seq_id, "".join(seq_chunks)


@dataclass
class Hit:
    contig: str
    prot_id: str
    model: str
    bitscore: float
    evalue: float
    domain_bitscore: Optional[float] = None
    hmm_coverage: Optional[float] = None
    domain_i_evalue: Optional[float] = None
    hmm_from: Optional[int] = None
    hmm_to: Optional[int] = None
    target_from: Optional[int] = None
    target_to: Optional[int] = None


# ---------- Core pipeline ----------


def _predict_proteins_pyrodigal(
    contigs_fa: Path,
    output_prot_fa: Path,
    mode: str = "meta",
    min_len: int = 90,
    min_protein_len_aa: int = 30,
    translation_table: int = 11,
    threads: int = 1,
) -> int:
    """Backward-compatible wrapper around shared prediction core."""
    inputs = PredictionInputs(
        input_contigs=contigs_fa,
        mode=mode,
        min_gene_len=min_len,
        min_protein_len_aa=min_protein_len_aa,
        translation_table=translation_table,
    )
    genes = predict_genes_pyrodigal(inputs, threads=threads)
    return write_predicted_proteins_fasta(genes, output_prot_fa)


def _hmmsearch(
    hmm_paths: list[Path],
    proteins_fa: Path,
    domtbl_paths: dict[str, Path],
    threads: int = 1,
    hmm_mode: str = "pure",
    keep_domtbl: bool = True,
    cut_ga: bool = True,
) -> Iterable[Hit]:
    """
    Run pyhmmer.hmmsearch on loaded HMMs and proteins.
    Returns hits directly as Hit objects and optionally writes domtbl files.
    """
    # Load all HMMs into memory
    hmms = []
    hmm_names = []
    for hmm_path in hmm_paths:
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            for hmm in hmm_file:
                hmms.append(hmm)
                hmm_names.append(hmm_path.stem)  # Use filename for pure mode

    # Load proteins into memory
    with pyhmmer.easel.SequenceFile(proteins_fa, digital=True) as seq_file:
        proteins = seq_file.read_block()

    # Run hmmsearch with pyHMMER
    bit_cutoffs = "gathering" if cut_ga else None
    try:
        hits_list = list(
            pyhmmer.hmmsearch(hmms, proteins, cpus=threads, bit_cutoffs=bit_cutoffs)
        )
    except Exception as exc:
        missing_cutoffs_exc = getattr(pyhmmer.plan7, "MissingCutoffs", None)
        if (
            cut_ga
            and missing_cutoffs_exc is not None
            and isinstance(exc, missing_cutoffs_exc)
        ):
            print(
                " Warning: One or more models are missing gathering cutoffs; "
                "retrying without --cut-ga."
            )
            hits_list = list(
                pyhmmer.hmmsearch(hmms, proteins, cpus=threads, bit_cutoffs=None)
            )
        else:
            raise

    # Write domtbl files if requested
    if keep_domtbl:
        for i, top_hits in enumerate(hits_list):
            # Map back to original HMM file for domtbl naming
            hmm_name = hmm_names[i] if i < len(hmm_names) else f"hmm_{i}"
            if hmm_name in domtbl_paths:
                domtbl_path = domtbl_paths[hmm_name]
                with domtbl_path.open("wb") as f:
                    top_hits.write(f, format="domains")

    # Process hits and yield Hit objects
    for i, top_hits in enumerate(hits_list):
        query_name = top_hits.query.name
        model_name = (
            query_name.decode() if isinstance(query_name, bytes) else query_name
        )

        # Determine model identifier based on hmm_mode
        if hmm_mode == "pure":
            # Use filename as model ID for pure mode
            model_id = hmm_names[i] if i < len(hmm_names) else model_name
        else:
            # Use actual HMM name for mixed mode
            model_id = model_name

        for hit in top_hits:
            if hit.included:  # pyHMMER's inclusion check
                prot_id = hit.name if isinstance(hit.name, str) else hit.name.decode()

                # Extract contig from prot_id with robust handling of multiple "|" characters
                # Expected format: "contig_name|gene<idx>" where contig_name may contain "|"
                # Use regex to find the last "|gene" pattern
                import re

                gene_pattern = r"\|gene\d+$"
                match = re.search(gene_pattern, prot_id)
                if match:
                    # Split at the position where "|gene" starts
                    contig = prot_id[: match.start()]
                else:
                    # Fallback: if no "|gene" pattern found, try simple split
                    if "|" in prot_id:
                        # Take everything before the last "|" as contig ID
                        contig = prot_id.rsplit("|", 1)[0]
                    else:
                        # No "|" found, use entire protein ID as contig ID
                        contig = prot_id

                domain_bitscore: Optional[float] = None
                hmm_coverage: Optional[float] = None
                domain_i_evalue: Optional[float] = None
                hmm_from: Optional[int] = None
                hmm_to: Optional[int] = None
                target_from: Optional[int] = None
                target_to: Optional[int] = None
                hit_domains = getattr(hit, "domains", None)
                if hit_domains is not None:
                    included_domains = [
                        domain
                        for domain in hit_domains
                        if getattr(domain, "included", True)
                    ]
                    if included_domains:
                        best_domain = max(included_domains, key=lambda domain: domain.score)
                        domain_bitscore = best_domain.score
                        domain_i_evalue = getattr(best_domain, "i_evalue", None)
                        model_length = getattr(top_hits.query, "M", None)
                        alignment = getattr(best_domain, "alignment", None)
                        if model_length and alignment is not None:
                            hmm_from = alignment.hmm_from
                            hmm_to = alignment.hmm_to
                            target_from = alignment.target_from
                            target_to = alignment.target_to
                            hmm_coverage = (
                                alignment.hmm_to - alignment.hmm_from + 1
                            ) / model_length

                yield Hit(
                    contig=contig,
                    prot_id=prot_id,
                    model=model_id,
                    bitscore=hit.score,
                    domain_bitscore=domain_bitscore,
                    evalue=hit.evalue,
                    hmm_coverage=hmm_coverage,
                    domain_i_evalue=domain_i_evalue,
                    hmm_from=hmm_from,
                    hmm_to=hmm_to,
                    target_from=target_from,
                    target_to=target_to,
                )


def _query_model_ids(hmm_paths: list[Path], hmm_mode: str) -> list[str]:
    """Return the complete, ordered model inventory before searching."""
    model_ids: list[str] = []
    sources: dict[str, Path] = {}

    for hmm_path in hmm_paths:
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            models = list(hmm_file)

        if hmm_mode == "pure" and len(models) != 1:
            raise ValueError(
                f"Pure HMM mode requires exactly one model per file: {hmm_path} "
                f"contains {len(models)}"
            )

        for model in models:
            model_name = model.name
            if isinstance(model_name, bytes):
                model_name = model_name.decode()
            model_id = hmm_path.stem if hmm_mode == "pure" else str(model_name)
            prev = sources.get(model_id)
            if prev is not None:
                raise ValueError(
                    f"Duplicate HMM model ID: {model_id} (seen in {prev} and {hmm_path})"
                )
            sources[model_id] = hmm_path
            model_ids.append(model_id)

    if not model_ids:
        raise ValueError("No HMM models found in the supplied files")
    return model_ids


def _effective_hit_score(hit: Hit, score_type: str) -> float:
    if score_type == "domain" and hit.domain_bitscore is not None:
        return hit.domain_bitscore
    return hit.bitscore


def _hit_sort_key(hit: Hit, score_type: str) -> tuple[float, float]:
    return (-_effective_hit_score(hit, score_type), hit.evalue)


def _filter_qualifying_hits(
    hits: Iterable[Hit],
    min_bitscore: Optional[float],
    max_evalue: Optional[float],
    kofam_metadata_by_model: Optional[dict[str, KOFamMetadata]] = None,
    use_kofam_thresholds: bool = True,
    dbcan_model_ids: Optional[Iterable[str]] = None,
) -> list[Hit]:
    kofam_metadata_by_model = kofam_metadata_by_model or {}
    dbcan_models = set(dbcan_model_ids or ())
    qualifying: list[Hit] = []
    for hit in hits:
        score_type = "full"
        ko_threshold: Optional[float] = None
        ko_meta = kofam_metadata_by_model.get(hit.model)
        if ko_meta is not None:
            score_type = ko_meta.score_type
            ko_threshold = ko_meta.threshold
        effective_score = _effective_hit_score(hit, score_type)
        if hit.model in dbcan_models and (
            hit.evalue >= 1e-15
            or hit.hmm_coverage is None
            or hit.hmm_coverage <= 0.35
        ):
            continue
        effective_min_bitscore = min_bitscore
        if use_kofam_thresholds and ko_threshold is not None:
            effective_min_bitscore = (
                ko_threshold
                if effective_min_bitscore is None
                else max(effective_min_bitscore, ko_threshold)
            )
        if effective_min_bitscore is not None and effective_score < effective_min_bitscore:
            continue
        if max_evalue is not None and hit.evalue > max_evalue:
            continue
        qualifying.append(hit)
    return qualifying


def _choose_best_contigs(
    hits: Iterable[Hit],
    min_bitscore: Optional[float],
    max_evalue: Optional[float],
    top_per_contig: int = 1,
    combine_mode: str = "any",
    min_hmm_hits: int = 1,
    total_hmm_models: int = 1,
    queried_model_ids: Optional[Iterable[str]] = None,
    kofam_metadata_by_model: Optional[dict[str, KOFamMetadata]] = None,
    use_kofam_thresholds: bool = True,
    dbcan_model_ids: Optional[Iterable[str]] = None,
) -> tuple[list[Hit], list[str]]:
    """
    Filter by thresholds, then pick top N hits per contig by bitscore.
    For KO models, ko_list thresholds are used by their score_type (full/domain).

    Returns (kept_hits, list_of_contig_ids).
    """
    kofam_metadata_by_model = kofam_metadata_by_model or {}
    queried_models = set(queried_model_ids or ())
    if queried_models:
        total_hmm_models = len(queried_models)

    per_contig: dict[str, list[Hit]] = defaultdict(list)
    for hit in _filter_qualifying_hits(
        hits,
        min_bitscore,
        max_evalue,
        kofam_metadata_by_model,
        use_kofam_thresholds,
        dbcan_model_ids,
    ):
        per_contig[hit.contig].append(hit)

    kept: list[Hit] = []
    kept_contigs: list[str] = []

    for contig, contig_hits in per_contig.items():
        if combine_mode == "any":
            if contig_hits:
                hits_per_model = defaultdict(list)
                for hit in contig_hits:
                    hits_per_model[hit.model].append(hit)

                for model_hits in hits_per_model.values():
                    score_type = "full"
                    if model_hits[0].model in kofam_metadata_by_model:
                        score_type = kofam_metadata_by_model[
                            model_hits[0].model
                        ].score_type
                    model_hits.sort(key=lambda hit: _hit_sort_key(hit, score_type))
                    kept.extend(model_hits[: max(1, top_per_contig)])
                kept_contigs.append(contig)

        elif combine_mode == "all":
            model_names = set(hit.model for hit in contig_hits)
            if queried_models:
                has_all_models = queried_models.issubset(model_names)
            else:
                has_all_models = len(model_names) == total_hmm_models
            if has_all_models:
                hits_per_model = defaultdict(list)
                for hit in contig_hits:
                    hits_per_model[hit.model].append(hit)

                for model_hits in hits_per_model.values():
                    score_type = "full"
                    if model_hits[0].model in kofam_metadata_by_model:
                        score_type = kofam_metadata_by_model[
                            model_hits[0].model
                        ].score_type
                    model_hits.sort(key=lambda hit: _hit_sort_key(hit, score_type))
                    kept.extend(model_hits[:1])

                kept_contigs.append(contig)

        elif combine_mode == "threshold":
            model_names = set(hit.model for hit in contig_hits)
            if len(model_names) >= min_hmm_hits:
                contig_hits.sort(
                    key=lambda hit: _hit_sort_key(
                        hit,
                        kofam_metadata_by_model[hit.model].score_type
                        if hit.model in kofam_metadata_by_model
                        else "full",
                    )
                )
                kept.extend(contig_hits[: max(1, top_per_contig)])
                kept_contigs.append(contig)

    return kept, kept_contigs


def _seqkit_extract(
    input_fa: Path,
    ids: list[str],
    output_fa: Path,
    seqkit_bin: str = "seqkit",
) -> None:
    if not ids:
        # write empty file but succeed, useful for pipelines
        output_fa.write_text("")
        return
    tmp = output_fa.parent / (output_fa.name + ".ids.txt")
    tmp.write_text("\n".join(ids) + "\n")
    cmd = [seqkit_bin, "grep", "-f", str(tmp), str(input_fa)]
    with output_fa.open("w") as out:
        p = subprocess.run(cmd, stdout=out, text=True, check=False)
    if p.returncode != 0:
        raise RuntimeError(f"seqkit grep failed with code {p.returncode}")
    tmp.unlink()  # cleanup


def _extract_target_proteins(
    kept_hits: list[Hit],
    kept_contig_ids: list[str],
    proteins_fa: Path,
    outdir: Path,
    hmm_mode: str,
    seqkit_bin: str = "seqkit",
    kofam_metadata_by_model: Optional[dict[str, KOFamMetadata]] = None,
) -> None:
    """
    Extract matched proteins per HMM model from the final kept contigs only.
    For KOFam models, headers are enriched as KO|definition.
    """
    # Create a set of kept contig IDs for fast lookup
    kept_contig_set = set(kept_contig_ids)
    kofam_metadata_by_model = kofam_metadata_by_model or {}

    # Group protein IDs by model - only from kept hits AND kept contigs
    proteins_per_model: dict[str, list[str]] = defaultdict(list)

    for hit in kept_hits:
        if hit.contig in kept_contig_set:
            proteins_per_model[hit.model].append(hit.prot_id)

    # Create target_proteins directory
    target_proteins_dir = outdir / "target_proteins"
    target_proteins_dir.mkdir(parents=True, exist_ok=True)

    # Extract proteins for each model
    for model_id, protein_ids in proteins_per_model.items():
        # Create a safe filename from the model identifier
        safe_model_name = re.sub(r"[^\w\-_.]", "_", model_id)
        output_path = target_proteins_dir / f"{safe_model_name}_proteins.mfa"

        if not protein_ids:
            output_path.write_text("")
            continue

        # Remove duplicates while preserving order
        unique_protein_ids = []
        seen = set()
        for pid in protein_ids:
            if pid not in seen:
                unique_protein_ids.append(pid)
                seen.add(pid)

        # Use seqkit to extract the proteins
        tmp_ids_file = output_path.parent / f"{output_path.name}.ids.tmp"
        tmp_ids_file.write_text("\n".join(unique_protein_ids) + "\n")

        cmd = [seqkit_bin, "grep", "-f", str(tmp_ids_file), str(proteins_fa)]

        with output_path.open("w") as out:
            result = subprocess.run(cmd, stdout=out, text=True, check=False)

        if result.returncode != 0:
            print(f"Warning: seqkit failed to extract proteins for {model_id}")
        else:
            ko_meta = kofam_metadata_by_model.get(model_id)
            if ko_meta is not None and ko_meta.definition:
                label = f"{ko_meta.ko_id}|{ko_meta.definition}"
                rewritten_path = output_path.parent / f"{output_path.name}.tmp"
                with output_path.open("r") as src, rewritten_path.open("w") as dst:
                    for line in src:
                        if line.startswith(">"):
                            dst.write(f"{line.rstrip()}|{label}\n")
                        else:
                            dst.write(line)
                rewritten_path.replace(output_path)

            print(
                f"    Extracted {len(unique_protein_ids)} proteins for {model_id} (from screened contigs)"
            )

        # Cleanup temporary file
        tmp_ids_file.unlink()


def _build_target_hmms(
    target_proteins_dir: Path,
    outdir: Path,
    threads: int = 1,
    aligner_bin: str = "mafft",
) -> None:
    """
    Build HMM models from target protein sequences using pyHMMER.
    Creates one HMM file per model from the corresponding protein FASTA files.

    For single sequences, builds HMM directly using builder.build().
    For multiple sequences, uses the configured external aligner before MSA creation.
    """
    target_hmms_dir = outdir / "target_hmms"
    target_hmms_dir.mkdir(parents=True, exist_ok=True)

    protein_files = list(target_proteins_dir.glob("*_proteins.mfa"))

    if not protein_files:
        print("    No target protein files found for HMM building")
        return

    print(f"    Building HMMs for {len(protein_files)} protein sets...")

    def _build_single_hmm(protein_file: Path) -> None:
        model_name = protein_file.stem.replace("_proteins", "")
        hmm_output_path = target_hmms_dir / f"{model_name}.hmm"

        try:
            alphabet = pyhmmer.easel.Alphabet.amino()
            builder = pyhmmer.plan7.Builder(alphabet)
            background = pyhmmer.plan7.Background(alphabet)

            sequences = [
                pyhmmer.easel.TextSequence(name=seq_id.encode(), sequence=seq_str)
                for seq_id, seq_str in _read_fasta(protein_file)
                if seq_str
            ]

            if len(sequences) == 0:
                print(f"      Skipping {model_name}: no valid sequences found")
                return
            elif len(sequences) == 1:
                # For single sequence, use builder.build() method
                digital_seq = sequences[0].digitize(alphabet)
                hmm, _, _ = builder.build(digital_seq, background)
                hmm.name = model_name.encode()
                print(f"      Built HMM from 1 sequence: {model_name}")
            else:
                aligned_path = protein_file.with_suffix(".aligned.faa")
                try:
                    with aligned_path.open("w") as aligned_output:
                        result = subprocess.run(
                            [
                                aligner_bin,
                                "--auto",
                                "--thread",
                                str(1 if threads and threads > 1 else max(1, threads)),
                                str(protein_file),
                            ],
                            stdout=aligned_output,
                            stderr=subprocess.PIPE,
                            text=True,
                            check=False,
                        )
                    if result.returncode != 0:
                        raise RuntimeError(
                            f"{aligner_bin} failed with code {result.returncode}: "
                            f"{result.stderr.strip()}"
                        )

                    aligned_sequences = [
                        pyhmmer.easel.TextSequence(
                            name=seq_id.encode(), sequence=seq_str
                        )
                        for seq_id, seq_str in _read_fasta(aligned_path)
                    ]
                    text_msa = pyhmmer.easel.TextMSA(
                        name=model_name.encode(), sequences=aligned_sequences
                    )
                    digital_msa = text_msa.digitize(alphabet)
                    hmm, _, _ = builder.build_msa(digital_msa, background)
                    print(
                        f"      Built HMM from {len(aligned_sequences)} aligned sequences: {model_name}"
                    )
                finally:
                    aligned_path.unlink(missing_ok=True)

            if not hmm.name:
                hmm.name = model_name.encode()
            hmm.command_line = None

            with hmm_output_path.open("wb") as f:
                hmm.write(f)

        except (OSError, RuntimeError, ValueError) as exc:
            print(f"      Warning: Failed to build HMM for {model_name}: {exc}")
            import traceback

            traceback.print_exc()

    # Use threads if requested; otherwise fall back to sequential processing
    if threads and threads > 1:
        print(f"    Using {threads} threads for HMM building")
        with ThreadPool(processes=threads) as pool:
            pool.map(_build_single_hmm, protein_files)
    else:
        for protein_file in protein_files:
            _build_single_hmm(protein_file)


def _screen(cfg: ScreenConfig) -> ScreenPlan:
    if cfg.all_cazymes or cfg.all_puls:
        _clean_bulk_output_artifacts(cfg.outdir)
        with tempfile.TemporaryDirectory(prefix="phu-dbcan-") as work_dir:
            return _screen_impl(cfg, Path(work_dir))
    return _screen_impl(cfg, None)


def _clean_bulk_output_artifacts(output_dir: Path) -> None:
    """Remove only known PHU-generated bulk intermediates and stale reports."""
    for name in (
        "resolved_dbcan_hmms",
        "cazyme_matches.tsv",
        "cazyme_matches.tsv.gz",
        "cazyme_summary.tsv",
        "cazyme_class_summary.tsv",
        "pul_matches.tsv",
        "pul_matches.tsv.gz",
        "pul_family_support.tsv.gz",
        "pul_summary.tsv",
        "substrate_summary.tsv",
        "query_manifest.json",
    ):
        path = output_dir / name
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()
    for path in output_dir.glob("hits_*.domtblout"):
        path.unlink()


def _screen_impl(cfg: ScreenConfig, bulk_work_dir: Path | None) -> ScreenPlan:
    """
    Screen contigs for protein families using pyHMMER.

    Main workflow:
    1. Predict proteins with pyrodigal
    2. Search proteins against HMMs with pyhmmer.hmmsearch
    3. Parse results, combine, and select best hits per contig
    4. Extract screened contigs
    5. Optionally extract matched proteins per HMM model from screened contigs only
    """
    if not cfg.input_contigs.exists():
        raise FileNotFoundError(f"Input file not found: {cfg.input_contigs}")

    # Resolve positional inputs: local HMM paths, PFAM accessions, and KO IDs.
    if cfg.all_cazymes:
        if cfg.hmms or cfg.all_puls:
            raise ValueError("--all-cazymes cannot be combined with other screen queries")
        cfg.hmms, query = _resolve_all_cazymes_inputs(
            bulk_work_dir or cfg.outdir
        )
        pul_rules = ()
        kofam_metadata_by_model = {}
    elif cfg.all_puls:
        if cfg.hmms:
            raise ValueError("--all-puls cannot be combined with positional queries")
        cfg.hmms, query, pul_rules = _resolve_all_puls_inputs(
            bulk_work_dir or cfg.outdir
        )
        kofam_metadata_by_model = {}
    else:
        cfg.hmms, kofam_metadata_by_model, query = _resolve_hmm_inputs(
            cfg.hmms, cfg.outdir, return_query=True
        )
        pul_rules = ()

    # Check all resolved HMM files exist.
    for hmm in cfg.hmms:
        if not hmm.exists():
            raise FileNotFoundError(f"HMM file not found: {hmm}")

    queried_model_ids = _query_model_ids(cfg.hmms, cfg.hmm_mode)

    # Discover binaries required by the selected output paths.
    seqkit_bin = _binaries()
    aligner_bin = _executable(["mafft"]) if cfg.save_target_hmms else ""
    plan = cfg.plan()
    if bulk_work_dir is not None:
        plan.domtbl_paths = {
            hmm.stem: bulk_work_dir / f"hits_{hmm.stem}.domtblout"
            for hmm in plan.hmms
        }
    if query.all_cazymes or query.all_puls:
        plan.keep_domtbl = False
    if query.pul_id is not None and cfg.combine_mode == "any":
        plan.combine_mode = "all"
    plan.hmmer_bin = ""  # Not used with pyHMMER
    plan.seqkit_bin = seqkit_bin

    # Create output directory
    plan.outdir.mkdir(parents=True, exist_ok=True)

    # Use cache-aware protein prediction
    pred_inputs = PredictionInputs(
        input_contigs=plan.input_contigs,
        mode=plan.mode,
        min_protein_len_aa=plan.min_protein_len_aa,
        translation_table=plan.translation_table,
    )
    cache_enabled = os.environ.get("PHU_CACHE", "on") != "off"
    cache_artifact = run_click_task(
        "Predicting proteins",
        get_or_predict_proteins,
        pred_inputs,
        use_cache=cache_enabled,
        threads=plan.threads,
    )

    print(
        "Predicting proteins with pyrodigal…"
        + (" [cache hit]" if cache_artifact.cache_hit else "")
    )
    print(f"  Proteins predicted: {cache_artifact.protein_count}")

    n_prot = cache_artifact.protein_count
    proteins_fa = cache_artifact.proteins_path

    if n_prot == 0:
        print("No proteins predicted. Exiting with empty outputs.")
        plan.out_contigs.write_text("")
        plan.kept_ids.write_text("")
        if query.all_puls:
            write_pul_outputs([], [], [], plan.outdir)
            write_cazyme_evidence([], [], plan.outdir / "evidence" / "cazyme_hits.tsv.gz")
            run_path = plan.outdir / ".phu" / "run.json"
            run_path.parent.mkdir(parents=True, exist_ok=True)
            run_path.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "phu_version": "0.9.0.dev0",
                        "command": "screen",
                        "mode": "all-puls",
                        "input": str(plan.input_contigs),
                        "database_name": "dbcan",
                        "database_hashes": _dbcan_hashes(),
                        "searched_family_count": len(query.normalized_models),
                        "raw_hmm_hit_count": 0,
                        "qualifying_hmm_hit_count": 0,
                        "contigs_with_cazyme_evidence": 0,
                        "complete_pul_match_count": 0,
                        "retained_contig_count": 0,
                        "generated_files": [
                            "screened_contigs.fasta", "kept_contigs.txt",
                            "pul_matches.tsv.gz", "pul_family_support.tsv.gz",
                            "pul_summary.tsv", "substrate_summary.tsv",
                            "evidence/cazyme_hits.tsv.gz", ".phu/run.json",
                        ],
                    },
                    indent=2,
                )
                + "\n"
            )
        if query.all_cazymes:
            write_cazyme_outputs([], [], [], plan.outdir)
            run_path = plan.outdir / ".phu" / "run.json"
            run_path.parent.mkdir(parents=True, exist_ok=True)
            run_path.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "phu_version": "0.9.0.dev0",
                        "command": "screen",
                        "mode": "all-cazymes",
                        "input": str(plan.input_contigs),
                        "database_name": "dbcan",
                        "database_version": query.database_version,
                        "database_digest": query.database_digest,
                        "database_hashes": _dbcan_hashes(),
                        "thresholds": {
                            "max_evalue": plan.max_evalue,
                            "min_bitscore": plan.min_bitscore,
                            "dbcan_independent_evalue": "<1e-15",
                            "dbcan_hmm_coverage": ">0.35",
                        },
                        "included_cazyme_classes": list(CAZYME_CLASSES),
                        "counts": {
                            "input_contigs": sum(1 for _ in _read_fasta(plan.input_contigs)),
                            "retained_contigs": 0,
                            "qualifying_hits": 0,
                            "matched_proteins": 0,
                            "matched_families": 0,
                            "matched_classes": 0,
                        },
                        "generated_files": [
                            "screened_contigs.fasta",
                            "kept_contigs.txt",
                            "cazyme_matches.tsv.gz",
                            "cazyme_summary.tsv",
                            "cazyme_class_summary.tsv",
                            "evidence/cazyme_hits.tsv.gz",
                            ".phu/run.json",
                        ],
                    },
                    indent=2,
                )
                + "\n"
            )
        if cache_artifact.temp_dir is not None:
            shutil.rmtree(cache_artifact.temp_dir, ignore_errors=True)
        return plan

    print(
        f"Running pyhmmer.hmmsearch for {len(plan.hmms)} HMM file(s) (mode: {plan.hmm_mode})…"
    )

    # Use pyHMMER for all searches
    all_hits = run_click_task(
        "Running HMM searches",
        lambda: list(
            _hmmsearch(
                hmm_paths=plan.hmms,
                proteins_fa=proteins_fa,
                domtbl_paths=plan.domtbl_paths,
                threads=plan.threads,
                hmm_mode=plan.hmm_mode,
                keep_domtbl=plan.keep_domtbl,
                cut_ga=plan.cut_ga,
            )
        ),
    )

    total_models = len(queried_model_ids)
    print(
        f"  {plan.hmm_mode.capitalize()} HMM mode: {total_models} models "
        f"from {len(plan.hmms)} files"
    )

    print(f"    Found {len(all_hits)} hits")

    print(
        f"Parsing results and selecting best hits per contig (combine_mode: {plan.combine_mode})…"
    )
    pul_matches: list[PULMatch] = []
    qualifying_cazyme_hits: list[Hit] = []
    if query.all_puls:
        qualifying_cazyme_hits = _filter_qualifying_hits(
            all_hits,
            plan.min_bitscore,
            plan.max_evalue,
            kofam_metadata_by_model,
            plan.use_kofam_thresholds,
            query.normalized_models,
        )
        filtered_hits, _ = _choose_best_contigs(
            all_hits,
            min_bitscore=plan.min_bitscore,
            max_evalue=plan.max_evalue,
            top_per_contig=plan.top_per_contig,
            combine_mode="any",
            min_hmm_hits=plan.min_hmm_hits,
            total_hmm_models=total_models,
            queried_model_ids=queried_model_ids,
            kofam_metadata_by_model=kofam_metadata_by_model,
            use_kofam_thresholds=plan.use_kofam_thresholds,
            dbcan_model_ids=query.normalized_models,
        )
        pul_matches, contig_ids = evaluate_pul_signatures(filtered_hits, pul_rules)
        kept_contig_set = set(contig_ids)
        kept_hits = [hit for hit in filtered_hits if hit.contig in kept_contig_set]
    elif query.all_cazymes:
        qualifying_cazyme_hits = _filter_qualifying_hits(
            all_hits,
            plan.min_bitscore,
            plan.max_evalue,
            kofam_metadata_by_model,
            plan.use_kofam_thresholds,
            query.normalized_models,
        )
        kept_hits, contig_ids = _choose_best_contigs(
            all_hits,
            min_bitscore=plan.min_bitscore,
            max_evalue=plan.max_evalue,
            top_per_contig=plan.top_per_contig,
            combine_mode="any",
            min_hmm_hits=plan.min_hmm_hits,
            total_hmm_models=total_models,
            queried_model_ids=queried_model_ids,
            kofam_metadata_by_model=kofam_metadata_by_model,
            use_kofam_thresholds=plan.use_kofam_thresholds,
            dbcan_model_ids=query.normalized_models,
        )
    else:
        kept_hits, contig_ids = _choose_best_contigs(
            all_hits,
            min_bitscore=plan.min_bitscore,
            max_evalue=plan.max_evalue,
            top_per_contig=plan.top_per_contig,
            combine_mode=plan.combine_mode,
            min_hmm_hits=plan.min_hmm_hits,
            total_hmm_models=total_models,
            queried_model_ids=queried_model_ids,
            kofam_metadata_by_model=kofam_metadata_by_model,
            use_kofam_thresholds=plan.use_kofam_thresholds,
            dbcan_model_ids=query.normalized_models,
        )

    if query.all_puls:
        write_pul_outputs(
            pul_matches,
            qualifying_cazyme_hits,
            contig_ids,
            plan.outdir,
        )
        write_cazyme_evidence(
            qualifying_cazyme_hits,
            contig_ids,
            plan.outdir / "evidence" / "cazyme_hits.tsv.gz",
        )
    if query.all_cazymes:
        write_cazyme_outputs(
            qualifying_cazyme_hits,
            kept_hits,
            contig_ids,
            plan.outdir,
        )

    plan.kept_ids.write_text("\n".join(contig_ids) + ("\n" if contig_ids else ""))

    # Extract target proteins per HMM if requested
    if plan.save_target_proteins:
        print("Extracting matched proteins per HMM model from screened contigs…")
        _extract_target_proteins(
            kept_hits,
            contig_ids,
            proteins_fa,
            plan.outdir,
            plan.hmm_mode,
            plan.seqkit_bin,
            kofam_metadata_by_model=kofam_metadata_by_model,
        )

        # Build HMMs from target proteins if requested
        if plan.save_target_hmms:
            print("Building HMMs from target protein sequences…")
            target_proteins_dir = plan.outdir / "target_proteins"
            _build_target_hmms(
                target_proteins_dir,
                plan.outdir,
                threads=plan.threads,
                aligner_bin=aligner_bin,
            )

    print(f"Extracting {len(contig_ids)} contig(s) with seqkit…")
    run_click_task(
        "Extracting contigs",
        _seqkit_extract,
        input_fa=plan.input_contigs,
        ids=contig_ids,
        output_fa=plan.out_contigs,
        seqkit_bin=plan.seqkit_bin,
    )

    # Output handling: copy proteins to output folder if requested
    if plan.keep_proteins:
        shutil.copy(proteins_fa, plan.outdir / "proteins.faa")
        write_prediction_metadata(
            plan.outdir / ".phu_prediction_metadata.json",
            cache_hit=cache_artifact.cache_hit,
            cache_key=cache_artifact.cache_key,
            cache_dir=cache_artifact.cache_dir,
        )

    # Clean up temp prediction directory (only set when caching is disabled)
    if cache_artifact.temp_dir is not None:
        shutil.rmtree(cache_artifact.temp_dir, ignore_errors=True)

    # Clean up domtblout if not requested
    if not plan.keep_domtbl:
        for domtbl_path in plan.domtbl_paths.values():
            if domtbl_path.exists():
                domtbl_path.unlink()

    database_hashes = _dbcan_hashes() if query.all_cazymes else {}
    query_manifest = {
        "original_query": list(query.original_inputs),
        "normalized_models": list(query.normalized_models),
        "pul_id": query.pul_id,
        "pul_substrate": query.pul_substrate,
        "pul_raw_rule": query.pul_raw_rule,
        "matching_rule": "threshold" if query.pul_id and cfg.combine_mode == "threshold" else query.matching_rule,
        "database_hashes": database_hashes,
    }
    if query.all_cazymes:
        by_class = defaultdict(int)
        for model in query.normalized_models:
            by_class[model.split("_", 1)[0].rstrip("0123456789")] += 1
        query_manifest.update(
            {
                "query_mode": "all_cazymes",
                "total_hmm_profiles": query.total_hmm_profiles,
                "canonical_cazyme_profiles_selected": len(query.normalized_models),
                "ancillary_profiles_excluded": list(query.excluded_ancillary_profiles),
                "cazyme_profiles_by_class": dict(by_class),
                "profiles_searched": list(query.normalized_models),
                "profiles_with_retained_hits": sorted({hit.model for hit in kept_hits}),
                "retained_contig_count": len(contig_ids),
                "retained_hit_count": len(kept_hits),
            }
        )
    if query.all_puls:
        query_manifest.update(
            {
                "all_puls": True,
                "total_pul_count": query.total_pul_count,
                "resolvable_pul_count": query.resolvable_pul_count,
                "skipped_unresolved_puls": list(query.skipped_puls),
                "unique_hmm_families_searched": len(query.normalized_models),
                "matched_pul_count": len(pul_matches),
                "retained_contig_count": len(contig_ids),
            }
        )
    if query.all_cazymes:
        run_manifest = {
            "schema_version": 1,
            "phu_version": "0.9.0.dev0",
            "command": "screen",
            "mode": "all-cazymes",
            "input": str(plan.input_contigs),
            "database_name": "dbcan",
            "database_version": query.database_version,
            "database_digest": query.database_digest,
            "database_hashes": database_hashes,
            "thresholds": {
                "max_evalue": plan.max_evalue,
                "min_bitscore": plan.min_bitscore,
                "dbcan_independent_evalue": "<1e-15",
                "dbcan_hmm_coverage": ">0.35",
            },
            "included_cazyme_classes": list(CAZYME_CLASSES),
            "counts": {
                "input_contigs": sum(1 for _ in _read_fasta(plan.input_contigs)),
                "retained_contigs": len(contig_ids),
                "qualifying_hits": len(qualifying_cazyme_hits),
                "matched_proteins": len({hit.prot_id for hit in kept_hits}),
                "matched_families": len({hit.model for hit in kept_hits}),
                "matched_classes": len({cazyme_class(hit.model) for hit in kept_hits}),
            },
            "generated_files": [
                "screened_contigs.fasta",
                "kept_contigs.txt",
                "cazyme_matches.tsv.gz",
                "cazyme_summary.tsv",
                "cazyme_class_summary.tsv",
                "evidence/cazyme_hits.tsv.gz",
                ".phu/run.json",
            ],
            "query": query_manifest,
        }
        manifest_path = plan.outdir / ".phu" / "run.json"
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(json.dumps(run_manifest, indent=2) + "\n")
    elif query.all_puls:
        run_manifest = {
            "schema_version": 1,
            "phu_version": "0.9.0.dev0",
            "command": "screen",
            "mode": "all-puls",
            "input": str(plan.input_contigs),
            "database_name": "dbcan",
            "database_hashes": _dbcan_hashes(),
            "searched_family_count": len(query.normalized_models),
            "raw_hmm_hit_count": len(all_hits),
            "qualifying_hmm_hit_count": len(qualifying_cazyme_hits),
            "contigs_with_cazyme_evidence": len({hit.contig for hit in qualifying_cazyme_hits}),
            "complete_pul_match_count": len(pul_matches),
            "retained_contig_count": len(contig_ids),
            "skipped_unresolved_puls": list(query.skipped_puls),
            "generated_files": [
                "screened_contigs.fasta", "kept_contigs.txt", "pul_matches.tsv.gz",
                "pul_family_support.tsv.gz", "pul_summary.tsv", "substrate_summary.tsv",
                "evidence/cazyme_hits.tsv.gz", ".phu/run.json",
            ],
        }
        manifest_path = plan.outdir / ".phu" / "run.json"
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(json.dumps(run_manifest, indent=2) + "\n")
    else:
        (plan.outdir / "query_manifest.json").write_text(
            json.dumps(query_manifest, indent=2) + "\n"
        )

    print(f"Done. Output FASTA: {plan.out_contigs}")
    files_msg = f"Also wrote: {plan.kept_ids.name} (contig IDs)"
    if plan.keep_domtbl:
        files_msg += f" and {len(plan.domtbl_paths)} domtblout files"
    if query.all_cazymes:
        files_msg += (
            " and cazyme_matches.tsv.gz, cazyme_summary.tsv, "
            "cazyme_class_summary.tsv, and evidence/cazyme_hits.tsv.gz"
        )
    if query.all_puls:
        files_msg += (
            " and pul_matches.tsv.gz, pul_family_support.tsv.gz, pul_summary.tsv, "
            "substrate_summary.tsv, and evidence/cazyme_hits.tsv.gz"
        )
    if plan.keep_proteins:
        proteins_msg = (
            "proteins.faa (cached)" if cache_artifact.cache_hit else "proteins.faa"
        )
        files_msg += f" and {proteins_msg}"
    if plan.save_target_proteins:
        files_msg += " and target proteins in target_proteins/ folder"
    if plan.save_target_hmms:
        files_msg += " and target HMMs in target_hmms/ folder"
    print(f"{files_msg}.")

    return plan
