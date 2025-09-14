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
from multiprocessing.pool import ThreadPool

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
    hmms: List[Path]  # Changed from hmm: Path to support multiple HMMs
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
    combine_mode: str = "any"  # New: how to combine hits from multiple HMMs
    min_hmm_hits: int = 1  # New: minimum number of HMMs that must hit a contig
    save_target_proteins: bool = False  # New: save matched proteins per HMM
    
    def __post_init__(self):
        """Validate configuration parameters."""
        if self.threads < 0:
            raise ValueError("threads must be >= 0")
        if self.threads == 0:
            # HMMER interprets --cpu 0 as "turn off multithreading"
            # This is valid, so we allow it
            pass
        if not self.hmms:
            raise ValueError("At least one HMM file must be provided")
        if self.combine_mode not in {"any", "all", "threshold"}:
            raise ValueError("combine_mode must be 'any', 'all', or 'threshold'")
    
    def plan(self) -> "ScreenPlan":
        """Create execution plan from configuration."""
        if self.mode not in {"meta", "single"}:
            raise ValueError("mode must be 'meta' or 'single'")
        
        effective_threads = self.threads
        
        # Create domtbl paths for each HMM
        domtbl_paths = {}
        target_protein_paths = {}
        for hmm in self.hmms:
            hmm_name = hmm.stem  # filename without extension
            domtbl_paths[hmm_name] = self.outdir / f"hits_{hmm_name}.domtblout"
            target_protein_paths[hmm_name] = self.outdir / "target_proteins" / f"{hmm_name}_proteins.mfa"
        
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
            combine_mode=self.combine_mode,
            min_hmm_hits=self.min_hmm_hits,
            save_target_proteins=self.save_target_proteins,
            proteins_fa=self.outdir / "proteins.faa",
            domtbl_paths=domtbl_paths,
            target_protein_paths=target_protein_paths,  # New
            kept_ids=self.outdir / "kept_contigs.txt",
            out_contigs=self.outdir / "screened_contigs.fasta",
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
    combine_mode: str
    min_hmm_hits: int
    save_target_proteins: bool  # New
    proteins_fa: Path
    domtbl_paths: Dict[str, Path]  # Changed from domtbl: Path
    target_protein_paths: Dict[str, Path]  # New
    kept_ids: Path
    out_contigs: Path

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

def _parse_domtblout(domtbl_path: Path, hmm_name: str) -> Iterable[Hit]:
    """
    Parse HMMER --domtblout (hmmsearch). Returns domain-level hits.
    Format spec: https://hmmer-web-docs.readthedocs.io/en/latest/output-files/domtab.html
    """
    with domtbl_path.open() as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = re.split(r"\s+", line.strip())
            try:
                target_name = cols[0]     # HMM/target name
                query_name  = cols[3]     # our protein sequence ID
                i_evalue    = float(cols[12])  # independent E-value
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
                contig_id = prot_id
            
            yield Hit(contig=contig_id, prot_id=prot_id, target=f"{hmm_name}:{target_name}", bitscore=bits, evalue=i_evalue)


# ---------- Core pipeline ----------

def _predict_proteins_pyrodigal(
    contigs_fa: Path,
    output_prot_fa: Path,
    mode: str = "meta",
    min_len: int = 90,
    translation_table: int = 11,
    threads: int = 1,
) -> int:
    """
    Use pyrodigal to predict CDS and write protein FASTA.
    Headers encode contig and CDS index as: contig|gene<idx>
    Returns number of proteins written.
    
    Uses ThreadPool for parallel processing of contigs when threads > 1.
    """
    # Initialize GeneFinder according to the API
    gf = GeneFinder(meta=(mode == "meta"), min_gene=min_len)
    
    # Read all contigs first
    contigs = list(_read_fasta(contigs_fa))
    
    if not contigs:
        return 0
    
    # Process contigs in parallel if threads > 1
    if threads > 1:
        # Extract sequences for parallel processing
        sequences = [seq for _, seq in contigs]
        
        # Use ThreadPool to process sequences in parallel
        with ThreadPool(processes=threads) as pool:
            genes_results = pool.map(gf.find_genes, sequences)
        
        # Write results
        n_prot = 0
        with output_prot_fa.open("w") as out:
            for (contig_id, _), genes in zip(contigs, genes_results):
                for i, gene in enumerate(genes, start=1):
                    aa = gene.translate()
                    if not aa:
                        continue
                    prot_id = f"{contig_id}|gene{i}"
                    out.write(f">{prot_id}\n{aa}\n")
                    n_prot += 1
    else:
        # Single-threaded processing (original logic)
        n_prot = 0
        with output_prot_fa.open("w") as out:
            for contig_id, seq in contigs:
                genes = gf.find_genes(seq)
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
    total_hmms: int = 1,
) -> Tuple[List[Hit], List[str]]:
    """
    Filter by thresholds, then pick top N hits per contig by bitscore.
    For multiple HMMs, apply combination logic.
    Returns (kept_hits, list_of_contig_ids).
    """
    from collections import defaultdict
    per_contig: Dict[str, List[Hit]] = defaultdict(list)

    # Filter hits by thresholds
    for h in hits:
        if min_bitscore is not None and h.bitscore < min_bitscore:
            continue
        if max_evalue is not None and h.evalue > max_evalue:
            continue
        per_contig[h.contig].append(h)

    # Apply combination logic
    kept: List[Hit] = []
    kept_contigs: List[str] = []
    
    for contig, contig_hits in per_contig.items():
        if combine_mode == "any":
            # Keep contigs that have hits from any HMM
            if contig_hits:
                contig_hits.sort(key=lambda x: (x.bitscore, -x.evalue), reverse=True)
                slice_ = contig_hits[:max(1, top_per_contig)]
                kept.extend(slice_)
                kept_contigs.append(contig)
        
        elif combine_mode == "all":
            # Keep contigs that have hits from ALL HMMs
            hmm_names = set(hit.target.split(':', 1)[0] for hit in contig_hits)
            if len(hmm_names) == total_hmms:  # Must have hits from ALL HMMs
                # For "all" mode, ensure we get exactly one hit per HMM per contig
                hits_per_hmm = defaultdict(list)
                for hit in contig_hits:
                    hmm_name = hit.target.split(':', 1)[0]
                    hits_per_hmm[hmm_name].append(hit)
                
                # Take the best hit for each HMM
                for hmm_name, hmm_hits in hits_per_hmm.items():
                    hmm_hits.sort(key=lambda x: (x.bitscore, -x.evalue), reverse=True)
                    kept.extend(hmm_hits[:1])  # Take only the best hit per HMM
                
                kept_contigs.append(contig)
        
        elif combine_mode == "threshold":
            # Keep contigs that have hits from at least min_hmm_hits HMMs
            hmm_names = set(hit.target.split(':', 1)[0] for hit in contig_hits)
            if len(hmm_names) >= min_hmm_hits:
                contig_hits.sort(key=lambda x: (x.bitscore, -x.evalue), reverse=True)
                slice_ = contig_hits[:max(1, top_per_contig)]
                kept.extend(slice_)
                kept_contigs.append(contig)

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


def _extract_target_proteins(
    kept_hits: List[Hit],
    kept_contig_ids: List[str],  # Add this parameter
    proteins_fa: Path,
    target_protein_paths: Dict[str, Path],
    seqkit_bin: str = "seqkit",
) -> None:
    """
    Extract matched proteins per HMM model from the final kept contigs only.
    This ensures the proteins match those from contigs in screened_contigs.fasta
    and respects the combine_mode filtering.
    """
    from collections import defaultdict
    
    # Create a set of kept contig IDs for fast lookup
    kept_contig_set = set(kept_contig_ids)
    
    # Group protein IDs by HMM model - only from kept hits AND kept contigs
    proteins_per_hmm: Dict[str, List[str]] = defaultdict(list)
    
    for hit in kept_hits:
        # Only include proteins from contigs that actually made it to the final output
        if hit.contig in kept_contig_set:
            hmm_name = hit.target.split(':', 1)[0]  # Extract HMM name from "hmm_name:target_name"
            proteins_per_hmm[hmm_name].append(hit.prot_id)
    
    # Create target_proteins directory
    if target_protein_paths:
        target_proteins_dir = next(iter(target_protein_paths.values())).parent
        target_proteins_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract proteins for each HMM
    for hmm_name, protein_ids in proteins_per_hmm.items():
        if hmm_name not in target_protein_paths:
            continue
            
        output_path = target_protein_paths[hmm_name]
        
        if not protein_ids:
            # Write empty file
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
            result = subprocess.run(cmd, stdout=out, text=True)
            
        if result.returncode != 0:
            print(f"Warning: seqkit failed to extract proteins for {hmm_name}")
        else:
            print(f"    Extracted {len(unique_protein_ids)} proteins for {hmm_name} (from screened contigs)")
            
        # Cleanup temporary file
        tmp_ids_file.unlink()


def _screen(cfg: ScreenConfig) -> ScreenPlan:
    """
    Screen contigs for protein families using HMMER.
    
    Main workflow:
    1. Predict proteins with pyrodigal
    2. Search proteins against each HMM with hmmsearch
    3. Parse results, combine, and select best hits per contig
    4. Extract screened contigs
    5. Optionally extract matched proteins per HMM model from screened contigs only
    """
    if not cfg.input_contigs.exists():
        raise FileNotFoundError(f"Input file not found: {cfg.input_contigs}")
    
    # Check all HMM files exist
    for hmm in cfg.hmms:
        if not hmm.exists():
            raise FileNotFoundError(f"HMM file not found: {hmm}")
    
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
        threads=plan.threads,  # Added threads parameter
    )
    print(f"  Proteins predicted: {n_prot}")
    
    if n_prot == 0:
        print("No proteins predicted. Exiting with empty outputs.")
        plan.out_contigs.write_text("")
        plan.kept_ids.write_text("")
        if not plan.keep_proteins and plan.proteins_fa.exists():
            plan.proteins_fa.unlink()
        return plan
    
    print(f"Running hmmsearch for {len(plan.hmms)} HMM(s)…")
    all_hits = []
    
    for hmm in plan.hmms:
        hmm_name = hmm.stem
        domtbl_path = plan.domtbl_paths[hmm_name]
        print(f"  Searching with {hmm.name}...")
        
        _hmmsearch(
            hmm=hmm,
            proteins_fa=plan.proteins_fa,
            domtbl_path=domtbl_path,
            threads=plan.threads,
            hmmer_bin=plan.hmmer_bin,
            extra_args=["--noali"],
        )
        
        hits = list(_parse_domtblout(domtbl_path, hmm_name))
        all_hits.extend(hits)
        print(f"    Found {len(hits)} hits")
    
    print(f"Parsing results and selecting best hits per contig (combine_mode: {plan.combine_mode})…")
    kept_hits, contig_ids = _choose_best_contigs(
        all_hits,
        min_bitscore=plan.min_bitscore,
        max_evalue=plan.max_evalue,
        top_per_contig=plan.top_per_contig,
        combine_mode=plan.combine_mode,
        min_hmm_hits=plan.min_hmm_hits,
        total_hmms=len(plan.hmms),  # Pass the total number of HMMs
    )
    
    plan.kept_ids.write_text("\n".join(contig_ids) + ("\n" if contig_ids else ""))
    
    # Extract target proteins per HMM if requested
    if plan.save_target_proteins:
        print("Extracting matched proteins per HMM model from screened contigs…")
        _extract_target_proteins(
            kept_hits,
            contig_ids,  # Pass the final list of kept contig IDs
            plan.proteins_fa,
            plan.target_protein_paths,
            plan.seqkit_bin,
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
        for domtbl_path in plan.domtbl_paths.values():
            if domtbl_path.exists():
                domtbl_path.unlink()
    
    print(f"Done. Output FASTA: {plan.out_contigs}")
    files_msg = f"Also wrote: {plan.kept_ids.name} (contig IDs)"
    if plan.keep_domtbl:
        files_msg += f" and {len(plan.domtbl_paths)} domtblout files"
    if plan.save_target_proteins:
        files_msg += f" and target proteins in target_proteins/ folder"
    print(f"{files_msg}.")
    
    return plan