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
from collections import defaultdict

import typer
from pyrodigal import GeneFinder   # pyrodigal>=3
from pyrodigal_gv import ViralGeneFinder  # New import for viral gene prediction

import pyhmmer
import pyhmmer.plan7
import pyhmmer.easel

from ._exec import run, _executable, CmdNotFound

app = typer.Typer(help="Screen contigs for a protein family using pyHMMER on predicted CDS.")


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
    save_target_hmms: bool = False  # New: build and save HMMs from target proteins
    hmm_mode: str = "pure"  # New: "pure" for single models, "mixed" for pressed/concatenated HMMs
    
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
        if self.hmm_mode not in {"pure", "mixed"}:
            raise ValueError("hmm_mode must be 'pure' or 'mixed'")
        if self.save_target_hmms and not self.save_target_proteins:
            raise ValueError("save_target_hmms requires save_target_proteins to be True")

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
            top_per_contig=self.top_per_contig,
            min_gene_len=self.min_gene_len,
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
    save_target_hmms: bool  # New
    hmm_mode: str  # New
    proteins_fa: Path
    domtbl_paths: Dict[str, Path]  # Changed from domtbl: Path
    kept_ids: Path
    out_contigs: Path

def _binaries() -> str:
    """
    Discover required binaries for screening.
    Only seqkit is needed since we use pyHMMER instead of HMMER binary.
    """
    seqkit = _executable(["seqkit"])
    return seqkit

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
    model: str
    bitscore: float
    evalue: float


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
    gf = ViralGeneFinder(meta=(mode == "meta"), min_gene=min_len)
    
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
    hmm_paths: List[Path],
    proteins_fa: Path,
    domtbl_paths: Dict[str, Path],
    threads: int = 1,
    hmm_mode: str = "pure",
    keep_domtbl: bool = True,
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
    hits_list = list(pyhmmer.hmmsearch(hmms, proteins, cpus=threads, bit_cutoffs=None))
    
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
        model_name = top_hits.query.name.decode()
        
        # Determine model identifier based on hmm_mode
        if hmm_mode == "pure":
            # Use filename as model ID for pure mode
            model_id = hmm_names[i] if i < len(hmm_names) else model_name
        else:
            # Use actual HMM name for mixed mode
            model_id = model_name
        
        for hit in top_hits:
            if hit.included:  # pyHMMER's inclusion check
                prot_id = hit.name.decode()
                
                # Extract contig from prot_id with robust handling of multiple "|" characters
                # Expected format: "contig_name|gene<idx>" where contig_name may contain "|"
                # Use regex to find the last "|gene" pattern
                import re
                gene_pattern = r'\|gene\d+$'
                match = re.search(gene_pattern, prot_id)
                if match:
                    # Split at the position where "|gene" starts
                    contig = prot_id[:match.start()]
                else:
                    # Fallback: if no "|gene" pattern found, try simple split
                    if "|" in prot_id:
                        # Take everything before the last "|" as contig ID
                        contig = prot_id.rsplit("|", 1)[0]
                    else:
                        # No "|" found, use entire protein ID as contig ID
                        contig = prot_id
                
                yield Hit(
                    contig=contig,
                    prot_id=prot_id,
                    model=model_id,
                    bitscore=hit.score,
                    evalue=hit.evalue
                )

def _choose_best_contigs(
    hits: Iterable[Hit],
    min_bitscore: Optional[float],
    max_evalue: Optional[float],
    top_per_contig: int = 1,
    combine_mode: str = "any",
    min_hmm_hits: int = 1,
    total_hmm_models: int = 1,
) -> Tuple[List[Hit], List[str]]:
    """
    Filter by thresholds, then pick top N hits per contig by bitscore.
    For multiple HMMs, apply combination logic.
    
    Returns (kept_hits, list_of_contig_ids).
    """
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
            if contig_hits:
                hits_per_model = defaultdict(list)
                for hit in contig_hits:
                    hits_per_model[hit.model].append(hit)

                # For each model, take the best hit(s) (by bitscore, then evalue)
                for model_hits in hits_per_model.values():
                    model_hits.sort(key=lambda x: (x.bitscore, -x.evalue), reverse=True)
                    kept.extend(model_hits[:max(1, top_per_contig)])
                kept_contigs.append(contig)
        
        elif combine_mode == "all":
            # Keep contigs that have hits from ALL models
            model_names = set(hit.model for hit in contig_hits)
            if len(model_names) == total_hmm_models:
                hits_per_model = defaultdict(list)
                for hit in contig_hits:
                    hits_per_model[hit.model].append(hit)
                
                # Take the best hit for each model
                for model_hits in hits_per_model.values():
                    model_hits.sort(key=lambda x: (x.bitscore, -x.evalue), reverse=True)
                    kept.extend(model_hits[:1])
                
                kept_contigs.append(contig)
        
        elif combine_mode == "threshold":
            # Keep contigs that have hits from at least min_hmm_hits models
            model_names = set(hit.model for hit in contig_hits)
            if len(model_names) >= min_hmm_hits:
                contig_hits.sort(key=lambda x: (x.bitscore, -x.evalue), reverse=True)
                kept.extend(contig_hits[:max(1, top_per_contig)])
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
    kept_contig_ids: List[str],
    proteins_fa: Path,
    outdir: Path,
    hmm_mode: str,
    seqkit_bin: str = "seqkit",
) -> None:
    """
    Extract matched proteins per HMM model from the final kept contigs only.
    """
    # Create a set of kept contig IDs for fast lookup
    kept_contig_set = set(kept_contig_ids)
    
    # Group protein IDs by model - only from kept hits AND kept contigs
    proteins_per_model: Dict[str, List[str]] = defaultdict(list)
    
    for hit in kept_hits:
        if hit.contig in kept_contig_set:
            proteins_per_model[hit.model].append(hit.prot_id)
    
    # Create target_proteins directory
    target_proteins_dir = outdir / "target_proteins"
    target_proteins_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract proteins for each model
    for model_id, protein_ids in proteins_per_model.items():
        # Create a safe filename from the model identifier
        safe_model_name = re.sub(r'[^\w\-_.]', '_', model_id)
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
            result = subprocess.run(cmd, stdout=out, text=True)
            
        if result.returncode != 0:
            print(f"Warning: seqkit failed to extract proteins for {model_id}")
        else:
            print(f"    Extracted {len(unique_protein_ids)} proteins for {model_id} (from screened contigs)")
            
        # Cleanup temporary file
        tmp_ids_file.unlink()

def _build_target_hmms(
    target_proteins_dir: Path,
    outdir: Path,
    threads: int = 1,
) -> None:
    """
    Build HMM models from target protein sequences using pyHMMER.
    Creates one HMM file per model from the corresponding protein FASTA files.
    
    For single sequences, builds HMM directly using builder.build().
    For multiple sequences, aligns them by padding to the same length before MSA creation.
    """
    # Create target_hmms directory
    target_hmms_dir = outdir / "target_hmms"
    target_hmms_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all protein FASTA files in target_proteins directory
    protein_files = list(target_proteins_dir.glob("*_proteins.mfa"))
    
    if not protein_files:
        print("    No target protein files found for HMM building")
        return
    
    # Initialize HMM builder and background
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    
    print(f"    Building HMMs for {len(protein_files)} protein sets...")
    
    for protein_file in protein_files:
        model_name = protein_file.stem.replace("_proteins", "")
        hmm_output_path = target_hmms_dir / f"{model_name}.hmm"
        
        try:
            # Read protein sequences as text first
            sequences = []
            with open(protein_file, 'r') as f:
                seq_id, seq_chunks = None, []
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if seq_id is not None and seq_chunks:
                            sequence = ''.join(seq_chunks)
                            if sequence:  # Only add non-empty sequences
                                text_seq = pyhmmer.easel.TextSequence(name=seq_id.encode(), sequence=sequence)
                                sequences.append(text_seq)
                        seq_id = line[1:].split()[0]
                        seq_chunks = []
                    else:
                        seq_chunks.append(line)
                
                # Don't forget the last sequence
                if seq_id is not None and seq_chunks:
                    sequence = ''.join(seq_chunks)
                    if sequence:
                        text_seq = pyhmmer.easel.TextSequence(name=seq_id.encode(), sequence=sequence)
                        sequences.append(text_seq)
            
            if len(sequences) == 0:
                print(f"      Skipping {model_name}: no valid sequences found")
                continue
            elif len(sequences) == 1:
                # For single sequence, use builder.build() method
                digital_seq = sequences[0].digitize(alphabet)
                hmm = builder.build(digital_seq, background)
                hmm.name = model_name.encode()
                print(f"      Built HMM from 1 sequence: {model_name}")
            else:
                # For multiple sequences, align by padding to same length
                max_len = max(len(seq.sequence) for seq in sequences)
                
                # Pad sequences to same length with gaps
                aligned_sequences = []
                for seq in sequences:
                    padded_seq = seq.sequence + '-' * (max_len - len(seq.sequence))
                    aligned_seq = pyhmmer.easel.TextSequence(name=seq.name, sequence=padded_seq)
                    aligned_sequences.append(aligned_seq)
                
                # Create MSA from aligned sequences
                text_msa = pyhmmer.easel.TextMSA(name=model_name.encode(), sequences=aligned_sequences)
                digital_msa = text_msa.digitize(alphabet)
                hmm, _, _ = builder.build_msa(digital_msa, background)
                print(f"      Built HMM from {len(sequences)} aligned sequences: {model_name}")
            
            # Ensure HMM has proper name
            if not hmm.name:
                hmm.name = model_name.encode()
            
            # Write HMM to file
            with hmm_output_path.open("wb") as f:
                hmm.write(f)
                
        except Exception as e:
            print(f"      Warning: Failed to build HMM for {model_name}: {e}")
            import traceback
            traceback.print_exc()
            continue

def _screen(cfg: ScreenConfig) -> ScreenPlan:
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
    
    # Check all HMM files exist
    for hmm in cfg.hmms:
        if not hmm.exists():
            raise FileNotFoundError(f"HMM file not found: {hmm}")
    
    # Discover binaries (only seqkit needed)
    seqkit_bin = _binaries()
    plan = cfg.plan()
    plan.hmmer_bin = ""  # Not used with pyHMMER
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
        threads=plan.threads,
    )
    print(f"  Proteins predicted: {n_prot}")
    
    if n_prot == 0:
        print("No proteins predicted. Exiting with empty outputs.")
        plan.out_contigs.write_text("")
        plan.kept_ids.write_text("")
        if not plan.keep_proteins and plan.proteins_fa.exists():
            plan.proteins_fa.unlink()
        return plan
    
    print(f"Running pyhmmer.hmmsearch for {len(plan.hmms)} HMM file(s) (mode: {plan.hmm_mode})…")
    
    # Use pyHMMER for all searches
    all_hits = list(_hmmsearch(
        hmm_paths=plan.hmms,
        proteins_fa=plan.proteins_fa,
        domtbl_paths=plan.domtbl_paths,
        threads=plan.threads,
        hmm_mode=plan.hmm_mode,
        keep_domtbl=plan.keep_domtbl,
    ))
    
    unique_model_ids = set(hit.model for hit in all_hits)
    if plan.hmm_mode == "pure":
        total_models = len(plan.hmms)  # Each file is one model
        print(f"  Pure HMM mode: {total_models} models from {len(plan.hmms)} files")
    else:
        total_models = len(unique_model_ids)  # Count actual models found
        print(f"  Mixed HMM mode: {total_models} unique models found from {len(plan.hmms)} files")
    
    print(f"    Found {len(all_hits)} hits")
    
    print(f"Parsing results and selecting best hits per contig (combine_mode: {plan.combine_mode})…")
    kept_hits, contig_ids = _choose_best_contigs(
        all_hits,
        min_bitscore=plan.min_bitscore,
        max_evalue=plan.max_evalue,
        top_per_contig=plan.top_per_contig,
        combine_mode=plan.combine_mode,
        min_hmm_hits=plan.min_hmm_hits,
        total_hmm_models=total_models,
    )
    
    plan.kept_ids.write_text("\n".join(contig_ids) + ("\n" if contig_ids else ""))
    
    # Extract target proteins per HMM if requested
    if plan.save_target_proteins:
        print("Extracting matched proteins per HMM model from screened contigs…")
        _extract_target_proteins(
            kept_hits,
            contig_ids,
            plan.proteins_fa,
            plan.outdir,
            plan.hmm_mode,
            plan.seqkit_bin,
        )
        
        # Build HMMs from target proteins if requested
        if plan.save_target_hmms:
            print("Building HMMs from target protein sequences…")
            target_proteins_dir = plan.outdir / "target_proteins"
            _build_target_hmms(
                target_proteins_dir,
                plan.outdir,
                threads=plan.threads,
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
    if plan.save_target_hmms:
        files_msg += f" and target HMMs in target_hmms/ folder"
    print(f"{files_msg}.")
    
    return plan