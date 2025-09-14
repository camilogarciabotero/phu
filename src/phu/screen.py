from __future__ import annotations
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import typer
from pyrodigal import GeneFinder   # pyrodigal>=3

from ._exec import run, _executable, CmdNotFound

app = typer.Typer(help="Screen contigs for a protein family using HMMER on predicted CDS.")


# ---------- Utilities ----------

def _cmd_exists(exe: str) -> bool:
    return shutil.which(exe) is not None

@dataclass
class ScreenConfig:
    """Configuration for screening contigs for protein families."""
    input_contigs: Path
    hmms: List[Path]  # Changed from hmm: Path
    outdir: Path = Path("phu-screen")
    mode: str = "meta"  # pyrodigal mode: meta|single
    threads: int = 1
    min_bitscore: Optional[float] = None
    max_evalue: Optional[float] = 1e-5
    top_per_contig: int = 1
    min_gene_len: int = 90
    translation_table: int = 11
    keep_proteins: bool = False
    keep_domtbl: bool = True
    combine_mode: str = "any"  # New: "any", "all", "threshold"
    min_hmm_hits: int = 1  # New: for threshold mode
    
    def __post_init__(self):
        """Validate configuration parameters."""
        if self.threads < 0:
            raise ValueError("threads must be >= 0")
        if self.threads == 0:
            # HMMER interprets --cpu 0 as "turn off multithreading"
            # This is valid, so we allow it
            pass
        if not self.hmms:
            raise ValueError("at least one HMM must be provided")
        for hmm in self.hmms:
            if not hmm.exists():
                raise FileNotFoundError(f"HMM file not found: {hmm}")
        if self.combine_mode not in {"any", "all", "threshold"}:
            raise ValueError("combine_mode must be 'any', 'all', or 'threshold'")
        if self.combine_mode == "threshold" and self.min_hmm_hits < 1:
            raise ValueError("min_hmm_hits must be >= 1 for threshold mode")
    
    def plan(self) -> "ScreenPlan":
        """Create execution plan from configuration."""
        if self.mode not in {"meta", "single"}:
            raise ValueError("mode must be 'meta' or 'single'")
        
        effective_threads = self.threads
        
        # Generate per-HMM domtbl paths
        domtbl_paths = {hmm.name: self.outdir / f"hits_{hmm.name}.domtblout" for hmm in self.hmms}
        
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
            top_per_contig=self.top_per_contig,
            min_gene_len=self.min_gene_len,
            translation_table=self.translation_table,
            keep_proteins=self.keep_proteins,
            keep_domtbl=self.keep_domtbl,
            proteins_fa=self.outdir / "proteins.faa",
            domtbl_paths=domtbl_paths,  # Changed from domtbl: Path
            kept_json=self.outdir / "kept_hits.json",
            kept_ids=self.outdir / "kept_contigs.txt",
            out_contigs=self.outdir / "screened_contigs.fasta",
            combine_mode=self.combine_mode,
            min_hmm_hits=self.min_hmm_hits,
        )


@dataclass
class ScreenPlan:
    """Execution plan for screening operation."""
    hmmer_bin: str
    seqkit_bin: str
    input_contigs: Path
    hmms: List[Path]  # Changed from hmm: Path
    outdir: Path
    mode: str
    threads: int
    min_bitscore: Optional[float]
    max_evalue: Optional[float]
    top_per_contig: int
    min_gene_len: int
    translation_table: int
    keep_proteins: bool
    keep_domtbl: bool
    proteins_fa: Path
    domtbl_paths: Dict[str, Path]  # Changed from domtbl: Path
    kept_json: Path
    kept_ids: Path
    out_contigs: Path
    combine_mode: str  # New
    min_hmm_hits: int  # New

def _binaries() -> tuple[str, str]:
    """
    Discover required binaries for screening.
    """
    hmmer = _executable(["hmmsearch"])
    seqkit = _executable(["seqkit"])
    return hmmer, seqkit

def _read_fasta(fp: Path) -> Iterable[Tuple[str, str]]:
    """Tiny FASTA reader (header up to first whitespace is the id)."""
    with fp.open() as f:
        seq_id, chunks = None, []
        for line in f:
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(chunks)
                seq_id = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if seq_id is not None:
            yield seq_id, "".join(chunks)

@dataclass
class Hit:
    contig: str
    prot_id: str
    target: str
    bitscore: float
    evalue: float
    hmm: str  # New: source HMM name

def _parse_domtblout(domtbl_path: Path, hmm_name: str) -> Iterable[Hit]:  # Added hmm_name param
    """
    Parse HMMER --domtblout (hmmsearch). Returns domain-level hits.
    Format spec: https://hmmer-web-docs.readthedocs.io/en/latest/output-files/domtab.html
    """
    with domtbl_path.open() as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            # HMMER domtblout is whitespace-separated with many columns.
            # We'll capture the essentials by positions (0-based):
            # 0=tname 3=qlen? (varies by tool), 17=i-Evalue, 13=bitscore, 18=env-from, 19=env-to, etc.
            # Safer: split and map by index.
            cols = re.split(r"\s+", line.strip())
            try:
                target_name = cols[0]     # HMM/target name
                query_name  = cols[3]     # our protein sequence ID
                i_evalue    = float(cols[12])  # independent E-value (index differs by versions; try 12 first)
                bits        = float(cols[13])
            except Exception:
                # Fallback for slight index shifts (older HMMER versions)
                try:
                    i_evalue = float(cols[11])
                    bits = float(cols[12])
                except Exception:
                    continue
            # Our protein IDs encode contig|gene index; extract original contig ID
            prot_id = target_name
            # Split only on the first "|" to handle contig names that might contain "|"
            if "|" in prot_id:
                contig_id = prot_id.split("|", 1)[0]
            else:
                # Fallback if no "|" found (shouldn't happen with our naming scheme)
                contig_id = prot_id
            yield Hit(contig=contig_id, prot_id=prot_id, target=target_name, bitscore=bits, evalue=i_evalue, hmm=hmm_name)  # Added hmm


# ---------- Core pipeline ----------

def _predict_proteins_pyrodigal(
    contigs_fa: Path,
    output_prot_fa: Path,
    mode: str = "meta",
    min_len: int = 90,
    translation_table: int = 11,
) -> int:
    """
    Use pyrodigal to predict CDS and write protein FASTA.
    Headers encode contig and CDS index as: contig|gene<idx>
    Returns number of proteins written.
    """
    # Initialize GeneFinder according to the API
    gf = GeneFinder(meta=(mode == "meta"), min_gene=min_len)
    
    # Note: translation_table is only used in train() method for single mode
    # In meta mode, pre-trained profiles are used, so no translation table needed
    
    n_prot = 0
    with output_prot_fa.open("w") as out:
        for contig_id, seq in _read_fasta(contigs_fa):
            genes = gf.find_genes(seq)
            # Translate CDS; pyrodigal returns Genes object that is directly iterable
            for i, gene in enumerate(genes, start=1):
                aa = gene.translate()
                if not aa:
                    continue
                prot_id = f"{contig_id}|gene{i}"
                out.write(f">{prot_id}\n{aa}\n")
                n_prot += 1
    return n_prot

def _hmmsearch(
    hmm: Path,
    proteins_fa: Path,
    domtbl_path: Path,
    threads: int = 1,
    hmmer_bin: str = "hmmsearch",
    extra_args: Optional[List[str]] = None,
) -> None:
    """Run hmmsearch with proper HMMER command structure."""
    # Validate inputs before running
    if not hmm.exists():
        raise FileNotFoundError(f"HMM file not found: {hmm}")
    if not proteins_fa.exists():
        raise FileNotFoundError(f"Protein FASTA file not found: {proteins_fa}")
    if proteins_fa.stat().st_size == 0:
        raise ValueError(f"Protein FASTA file is empty: {proteins_fa}")
    
    # Build command according to HMMER documentation:
    # hmmsearch [options] <hmmfile> <seqfile>
    cmd = [hmmer_bin]
    
    # Add all options first
    # Note: --cpu <n> specifies number of worker threads
    # HMMER will spawn <n>+1 total threads (workers + master)
    cmd.extend(["--cpu", str(threads)])
    cmd.extend(["--domtblout", str(domtbl_path)])
    
    if extra_args:
        cmd.extend(extra_args)
    
    # Add positional arguments: HMM file first, then sequence file
    cmd.append(str(hmm))
    cmd.append(str(proteins_fa))
    
    # Use subprocess directly for better error handling
    try:
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            check=True
        )
        # Log successful completion for debugging
        if result.returncode == 0:
            print(f"  hmmsearch completed successfully using {threads} worker threads")
    except subprocess.CalledProcessError as e:
        print(f"hmmsearch failed with exit code {e.returncode}")
        print(f"Command: {' '.join(cmd)}")
        if e.stdout:
            print(f"STDOUT:\n{e.stdout}")
        if e.stderr:
            print(f"STDERR:\n{e.stderr}")
        raise RuntimeError(f"hmmsearch failed: {e.stderr}") from e

def _choose_best_contigs(
    hits: Iterable[Hit],
    min_bitscore: Optional[float],
    max_evalue: Optional[float],
    top_per_contig: int = 1,
    combine_mode: str = "any",
    min_hmm_hits: int = 1,
) -> Tuple[List[Hit], List[str]]:
    """
    Filter by thresholds, combine hits per contig based on mode, then pick top N hits per contig by bitscore.
    Returns (kept_hits, list_of_contig_ids).
    """
    from collections import defaultdict
    per_contig: Dict[str, List[Hit]] = defaultdict(list)

    for h in hits:
        if min_bitscore is not None and h.bitscore < min_bitscore:
            continue
        if max_evalue is not None and h.evalue > max_evalue:
            continue
        per_contig[h.contig].append(h)

    kept: List[Hit] = []
    kept_contigs: List[str] = []
    for contig, lst in per_contig.items():
        # Group hits by HMM
        hmm_hits = defaultdict(list)
        for h in lst:
            hmm_hits[h.hmm].append(h)
        
        # Apply combine logic
        if combine_mode == "any":
            if hmm_hits:  # At least one HMM hit
                kept.extend(lst)
                kept_contigs.append(contig)
        elif combine_mode == "all":
            if len(hmm_hits) == len(set(h.hmm for h in lst)):  # Hit all HMMs
                kept.extend(lst)
                kept_contigs.append(contig)
        elif combine_mode == "threshold":
            if len(hmm_hits) >= min_hmm_hits:
                kept.extend(lst)
                kept_contigs.append(contig)
        
        # If kept, sort and slice top_per_contig (but keep all for now, as per original)
        if contig in kept_contigs:
            lst.sort(key=lambda x: (x.bitscore, -x.evalue), reverse=True)
            kept = [h for h in kept if h.contig == contig][:max(1, top_per_contig)]  # Adjust kept list

    return kept, kept_contigs

def _seqkit_extract(
    input_fa: Path,
    ids: List[str],
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
        p = subprocess.run(cmd, stdout=out, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"seqkit grep failed with code {p.returncode}")
    tmp.unlink()  # cleanup


def _screen(cfg: ScreenConfig) -> ScreenPlan:
    """
    Screen contigs for a protein family using HMMER.
    
    Main workflow:
    1. Predict proteins with pyrodigal
    2. Search proteins against each HMM with hmmsearch
    3. Parse results, combine, and select best hits per contig
    4. Extract screened contigs
    """
    if not cfg.input_contigs.exists():
        raise FileNotFoundError(f"Input file not found: {cfg.input_contigs}")
    
    # Discover binaries
    hmmer_bin, seqkit_bin = _binaries()
    plan = cfg.plan()
    plan.hmmer_bin = hmmer_bin
    plan.seqkit_bin = seqkit_bin
    
    # Create output directory
    plan.outdir.mkdir(parents=True, exist_ok=True)
    
    print("Predicting proteins with pyrodigal…")
    n_prot = _predict_proteins_pyrodigal(
        plan.input_contigs,
        plan.proteins_fa,
        mode=plan.mode,
        min_len=plan.min_gene_len,
        translation_table=plan.translation_table,
    )
    print(f"  Proteins predicted: {n_prot}")
    
    if n_prot == 0:
        print("No proteins predicted. Exiting with empty outputs.")
        plan.out_contigs.write_text("")
        plan.kept_ids.write_text("")
        plan.kept_json.write_text(json.dumps([], indent=2))
        if not plan.keep_proteins and plan.proteins_fa.exists():
            plan.proteins_fa.unlink()
        return plan
    
    print("Running hmmsearch for each HMM…")
    all_hits = []
    for hmm in plan.hmms:
        domtbl_path = plan.domtbl_paths[hmm.name]
        _hmmsearch(
            hmm=hmm,
            proteins_fa=plan.proteins_fa,
            domtbl_path=domtbl_path,
            threads=plan.threads,
            hmmer_bin=plan.hmmer_bin,
            extra_args=["--noali"],
        )
        hits = list(_parse_domtblout(domtbl_path, hmm.name))
        all_hits.extend(hits)
    
    print("Parsing results and selecting best hits per contig…")
    kept_hits, contig_ids = _choose_best_contigs(
        all_hits,
        min_bitscore=plan.min_bitscore,
        max_evalue=plan.max_evalue,
        top_per_contig=plan.top_per_contig,
        combine_mode=plan.combine_mode,
        min_hmm_hits=plan.min_hmm_hits,
    )
    
    plan.kept_ids.write_text("\n".join(contig_ids) + ("\n" if contig_ids else ""))
    plan.kept_json.write_text(
        json.dumps([h.__dict__ for h in kept_hits], indent=2)
    )
    
    print(f"Extracting {len(contig_ids)} contig(s) with seqkit…")
    _seqkit_extract(
        input_fa=plan.input_contigs,
        ids=contig_ids,
        output_fa=plan.out_contigs,
        seqkit_bin=plan.seqkit_bin,
    )
    
    # Clean up if requested
    if not plan.keep_proteins and plan.proteins_fa.exists():
        plan.proteins_fa.unlink()
    if not plan.keep_domtbl:
        for path in plan.domtbl_paths.values():
            if path.exists():
                path.unlink()
    
    print(f"Done. Output FASTA: {plan.out_contigs}")
    files_msg = f"Also wrote: {plan.kept_ids.name} (contig IDs), {plan.kept_json.name} (hits JSON)"
    if plan.keep_domtbl:
        files_msg += f" and {len(plan.hmms)} domtblout files"
    print(f"{files_msg}.")
    
    return plan