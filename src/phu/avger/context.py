"""Genomic context of avger candidates.

For every candidate gene (predicted or filtered) avger takes a window of up to
``flank`` genes on each side, summarises the V/VL-scores of the flanking genes
and groups the windows into gene cluster families (GCFs) with IGUA
(Larralde et al.), run on MMseqs2 as bundled by pymmseqs.
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass, field
from importlib import resources
from importlib.metadata import version
from pathlib import Path
from typing import Iterable, Mapping

from phu.core.gene_prediction import PredictedGene

DATABASES = ("pfam", "kofam", "phrog")

NEIGHBORHOOD_COLUMNS = (
    "protein_id",
    "contig",
    "gene_number",
    "decision",
    "window",
    "window_start",
    "window_end",
    "genes_left",
    "genes_right",
    "gcf_id",
    "gcf_windows",
    "gcf_contigs",
    "gcf_predictions",
)

FLANK_COLUMNS = (
    "protein_id",
    "side",
    "database",
    "genes",
    "scored_genes",
    "max_v",
    "mean_vl",
)


@dataclass
class ContextResult:
    """Context tables plus the per-anchor family, for the main tables."""

    neighborhoods: list[dict[str, object]]
    flanks: list[dict[str, object]]
    gcf_by_protein: dict[str, dict[str, object]]
    versions: dict[str, str] = field(default_factory=dict)
    gcfs: int = 0


def mmseqs_binary() -> str:
    """Path of the MMseqs2 binary bundled with pymmseqs."""
    return str(resources.files("pymmseqs") / "bin" / "mmseqs")


def context_versions() -> dict[str, str]:
    return {name: version(name) for name in ("igua", "pymmseqs", "gb-io")}


def _contig_sequences(path: Path, wanted: set[str]) -> dict[str, str]:
    from phu.core.gene_prediction import _read_fasta_python

    return {name: seq for name, seq in _read_fasta_python(path) if name in wanted}


def _flank_summary(
    genes: list[PredictedGene], scores: Mapping[str, Mapping[str, tuple]]
) -> dict[str, dict[str, object]]:
    """Genes, scored genes, max V and mean VL per database for one flank.

    A flank without genes reports 0 genes and empty statistics; genes without
    a score count as genes but not as scored genes.
    """
    out = {}
    for db in DATABASES:
        pairs = [scores[g.gene_id][db] for g in genes if db in scores.get(g.gene_id, {})]
        v = [p[0] for p in pairs]
        vl = [p[1] for p in pairs]
        out[db] = {
            "genes": len(genes),
            "scored_genes": len(pairs),
            "max_v": max(v) if v else "",
            "mean_vl": round(sum(vl) / len(vl), 4) if vl else "",
        }
    return out


def _run_igua(records: list, workdir: Path, threads: int):
    """Run IGUA on GenBank records; returns its GCF membership table."""
    import gb_io
    import rich.console
    import rich.progress
    from igua.dataset.genbank import GenBankDataset
    from igua.mmseqs import MMSeqs
    from igua.pipeline import Pipeline

    gbk = workdir / "windows.gbk"
    with gbk.open("wb") as handle:
        gb_io.dump(records, handle)
    # IGUA's readers need a progress object; give it a disabled one. Its
    # MMseqs2 driver does not inherit the pipeline's threads or workdir.
    quiet = rich.progress.Progress(
        console=rich.console.Console(quiet=True), disable=True
    )
    tmp = workdir / "mmseqs-tmp"
    tmp.mkdir()
    result = Pipeline(
        jobs=threads,
        progress=quiet,
        mmseqs=MMSeqs(
            binary=mmseqs_binary(), progress=None, threads=threads, tempdir=str(tmp)
        ),
        workdir=workdir,
    ).run(GenBankDataset(gbk))
    return result.gcfs


def build_context(
    *,
    contigs_path: Path,
    genes: Iterable[PredictedGene],
    anchors: Mapping[str, str],
    scores: Mapping[str, Mapping[str, tuple]],
    flank: int,
    threads: int,
    workdir: Path,
) -> tuple[ContextResult, list]:
    """Windows, flank summaries and GCFs for the anchor genes.

    ``anchors`` maps protein id -> decision (an avg class or a filter
    status); ``scores`` maps protein id -> database -> (V, VL). Returns the
    result and the GenBank records of the windows.
    """
    import gb_io

    by_contig: dict[str, list[PredictedGene]] = {}
    for gene in genes:
        by_contig.setdefault(gene.contig_id, []).append(gene)
    for contig_genes in by_contig.values():
        contig_genes.sort(key=lambda g: g.ordinal)
    position = {
        g.gene_id: (g.contig_id, i)
        for contig_genes in by_contig.values()
        for i, g in enumerate(contig_genes)
    }
    anchor_ids = [a for a in sorted(anchors) if a in position]
    sequences = _contig_sequences(
        contigs_path, {position[a][0] for a in anchor_ids}
    )

    neighborhoods: list[dict[str, object]] = []
    flanks: list[dict[str, object]] = []
    records = []
    for n, anchor in enumerate(anchor_ids, start=1):
        contig, idx = position[anchor]
        contig_genes = by_contig[contig]
        first = max(0, idx - flank)
        window = contig_genes[first : idx + flank + 1]
        left, right = contig_genes[first:idx], contig_genes[idx + 1 : idx + flank + 1]
        start = min(min(g.start, g.end) for g in window)  # 1-based, inclusive
        end = max(max(g.start, g.end) for g in window)
        name = f"W{n:07d}"
        features = []
        for g in window:
            lo, hi = min(g.start, g.end), max(g.start, g.end)
            location = gb_io.Range(lo - start, hi - start + 1)
            if g.strand < 0:
                location = gb_io.Complement(location)
            features.append(
                gb_io.Feature(
                    "CDS",
                    location,
                    qualifiers=[
                        gb_io.Qualifier("locus_tag", g.gene_id),
                        gb_io.Qualifier(
                            "translation", g.amino_acid_sequence.rstrip("*")
                        ),
                    ],
                )
            )
        records.append(
            gb_io.Record(
                sequences[contig][start - 1 : end].encode(),
                name=name,
                molecule_type="DNA",
                features=features,
            )
        )
        neighborhoods.append(
            {
                "protein_id": anchor,
                "contig": contig,
                "gene_number": contig_genes[idx].ordinal,
                "decision": anchors[anchor],
                "window": name,
                "window_start": start,
                "window_end": end,
                "genes_left": len(left),
                "genes_right": len(right),
            }
        )
        for side, side_genes in (("left", left), ("right", right)):
            for db, stats in _flank_summary(side_genes, scores).items():
                flanks.append(
                    {"protein_id": anchor, "side": side, "database": db, **stats}
                )

    gcf_of: dict[str, str] = {}
    if len(records) == 1:
        gcf_of[records[0].name] = "GCF0000001"
    elif records:
        table = _run_igua(records, workdir, threads)
        gcf_of = dict(zip(table["cluster_id"], table["gcf_id"]))
        missing = {r.name for r in records} - set(gcf_of)
        if missing:
            raise RuntimeError(
                f"IGUA returned no family for {len(missing)} of {len(records)} windows"
            )

    members: dict[str, list[dict[str, object]]] = {}
    for row in neighborhoods:
        row["gcf_id"] = gcf_of.get(row["window"], "")
        members.setdefault(row["gcf_id"], []).append(row)
    gcf_by_protein = {}
    for gcf, rows in members.items():
        stats = {
            "gcf_windows": len(rows),
            "gcf_contigs": len({r["contig"] for r in rows}),
            "gcf_predictions": sum(r["decision"].startswith("putative_") for r in rows),
        }
        for r in rows:
            r.update(stats)
            gcf_by_protein[r["protein_id"]] = {"gcf_id": gcf, **stats}
    return (
        ContextResult(
            neighborhoods=neighborhoods,
            flanks=flanks,
            gcf_by_protein=gcf_by_protein,
            versions=context_versions(),
            gcfs=len(members),
        ),
        records,
    )


def write_genbank(records: list, path: Path) -> None:
    import gb_io

    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(delete=False, dir=path.parent) as handle:
        gb_io.dump(records, handle)
    Path(handle.name).replace(path)

