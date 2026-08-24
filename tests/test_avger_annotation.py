from __future__ import annotations

from pathlib import Path

import pytest

from phu import avger_annotation as ann
from phu.kofam_db import KOFamMetadata


class FakeAlignment:
    def __init__(self, hmm_from=1, hmm_to=10, target_from=2, target_to=20):
        self.hmm_from = hmm_from
        self.hmm_to = hmm_to
        self.target_from = target_from
        self.target_to = target_to


class FakeDomain:
    def __init__(self, score: float, included: bool = True, alignment: FakeAlignment | None = None):
        self.score = score
        self.included = included
        self.alignment = alignment or FakeAlignment()


class FakeHit:
    def __init__(self, name: str, score: float, evalue: float, domains: list[FakeDomain]):
        self.name = name.encode()
        self.score = score
        self.evalue = evalue
        self.domains = domains


class FakeCutoffs:
    def __init__(self, available: bool, ga1=0.0, ga2=0.0):
        self._available = available
        self.gathering1 = ga1
        self.gathering2 = ga2

    def gathering_available(self):
        return self._available


class FakeHMM:
    def __init__(self, name: str, accession: str, cutoffs: FakeCutoffs | None = None):
        self.name = name.encode()
        self.accession = accession.encode()
        self.cutoffs = cutoffs


class FakeTopHits:
    def __init__(self, query: FakeHMM, hits: list[FakeHit]):
        self.query = query
        self._hits = hits

    def __iter__(self):
        return iter(self._hits)


def _patch_hmmsearch(monkeypatch, outputs):
    def fake_hmmsearch(*args, **kwargs):
        return iter(outputs)

    monkeypatch.setattr(ann.pyhmmer, "hmmsearch", fake_hmmsearch)


def _patch_hmmfile(monkeypatch):
    class FakeHMMFile:
        def __init__(self, path):
            self.path = path

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr(ann.pyhmmer.plan7, "HMMFile", FakeHMMFile)


def _patch_sequencefile(monkeypatch):
    class FakeSequenceFile:
        def __init__(self, path, digital=True):
            self.path = path
            self.digital = digital

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def read_block(self, sequences=None):
            if getattr(self, "_read", False):
                return []
            self._read = True
            return [object()]

    monkeypatch.setattr(ann.pyhmmer.easel, "SequenceFile", FakeSequenceFile)


def test_parse_contig_id_from_protein_id_pipe_and_gene_suffix():
    assert ann.parse_contig_id_from_protein_id("NODE_1|len100|gene3") == "NODE_1|len100"
    assert ann.parse_contig_id_from_protein_id("NODE_1|len100|protA") == "NODE_1|len100"
    assert ann.parse_contig_id_from_protein_id("contig42") == "contig42"


def test_pfam_skips_missing_ga_by_default(monkeypatch):
    _patch_hmmfile(monkeypatch)
    _patch_sequencefile(monkeypatch)
    _patch_hmmsearch(
        monkeypatch,
        [
            FakeTopHits(
                query=FakeHMM("MISSINGGA", "PF00001.1", cutoffs=FakeCutoffs(False)),
                hits=[
                    FakeHit(
                        "ctgA|gene1",
                        score=120.0,
                        evalue=1e-20,
                        domains=[FakeDomain(80.0)],
                    )
                ],
            )
        ],
    )

    cfg = ann.AnnotationConfig(max_evalue=1e-5, pfam_require_ga=True, pfam_missing_ga_policy="skip_model")
    best, passing, models, skipped = ann._search_and_collect(
        database="pfam",
        hmm_path=Path("pfam.hmm"),
        proteins_path=Path("proteins.faa"),
        cfg=cfg,
    )

    assert models == 1
    assert skipped == 1
    assert passing == 0
    assert best == {}


def test_pfam_include_without_ga_policy(monkeypatch):
    _patch_hmmfile(monkeypatch)
    _patch_sequencefile(monkeypatch)
    _patch_hmmsearch(
        monkeypatch,
        [
            FakeTopHits(
                query=FakeHMM("NOGA", "PF00002.9", cutoffs=FakeCutoffs(False)),
                hits=[FakeHit("ctgA|gene1", score=55.0, evalue=1e-12, domains=[FakeDomain(30.0)])],
            )
        ],
    )

    cfg = ann.AnnotationConfig(
        max_evalue=1e-5,
        pfam_require_ga=True,
        pfam_missing_ga_policy="include_without_ga",
    )
    best, passing, models, skipped = ann._search_and_collect(
        database="pfam",
        hmm_path=Path("pfam.hmm"),
        proteins_path=Path("proteins.faa"),
        cfg=cfg,
    )

    assert models == 1
    assert skipped == 0
    assert passing == 1
    rec = best["ctgA|gene1"]
    assert rec.model_accession == "PF00002"
    assert rec.threshold_value is None


def test_pfam_best_hit_is_highest_effective_score_then_evalue_then_accession(monkeypatch):
    _patch_hmmfile(monkeypatch)
    _patch_sequencefile(monkeypatch)
    _patch_hmmsearch(
        monkeypatch,
        [
            FakeTopHits(
                query=FakeHMM("A", "PF00010.1", cutoffs=FakeCutoffs(True, ga1=20.0, ga2=10.0)),
                hits=[FakeHit("ctgA|gene1", score=40.0, evalue=1e-8, domains=[FakeDomain(15.0)])],
            ),
            FakeTopHits(
                query=FakeHMM("B", "PF00009.1", cutoffs=FakeCutoffs(True, ga1=20.0, ga2=10.0)),
                hits=[FakeHit("ctgA|gene1", score=40.0, evalue=1e-8, domains=[FakeDomain(12.0)])],
            ),
        ],
    )

    cfg = ann.AnnotationConfig(max_evalue=1e-5)
    best, passing, _, _ = ann._search_and_collect(
        database="pfam",
        hmm_path=Path("pfam.hmm"),
        proteins_path=Path("proteins.faa"),
        cfg=cfg,
    )

    assert passing == 2
    assert best["ctgA|gene1"].model_accession == "PF00009"


def test_kofam_domain_vs_full_threshold_logic(monkeypatch):
    _patch_hmmfile(monkeypatch)
    _patch_sequencefile(monkeypatch)
    _patch_hmmsearch(
        monkeypatch,
        [
            FakeTopHits(
                query=FakeHMM("K00001", "K00001"),
                hits=[FakeHit("ctgA|gene1", score=40.0, evalue=1e-30, domains=[FakeDomain(80.0)])],
            ),
            FakeTopHits(
                query=FakeHMM("K00002", "K00002"),
                hits=[FakeHit("ctgA|gene2", score=70.0, evalue=1e-30, domains=[FakeDomain(45.0)])],
            ),
        ],
    )

    cfg = ann.AnnotationConfig(max_evalue=1e-5)
    kofam_meta = {
        "K00001": KOFamMetadata(
            ko_id="K00001",
            threshold=50.0,
            score_type="domain",
            profile_type="all families",
            definition="a",
        ),
        "K00002": KOFamMetadata(
            ko_id="K00002",
            threshold=60.0,
            score_type="full",
            profile_type="all families",
            definition="b",
        ),
    }

    best, passing, models, skipped = ann._search_and_collect(
        database="kofam",
        hmm_path=Path("kofam.hmm"),
        proteins_path=Path("proteins.faa"),
        cfg=cfg,
        kofam_meta_by_model=kofam_meta,
    )

    assert models == 2
    assert skipped == 0
    assert passing == 2
    assert best["ctgA|gene1"].effective_score == pytest.approx(80.0)
    assert best["ctgA|gene1"].score_type == "domain"
    assert best["ctgA|gene2"].effective_score == pytest.approx(70.0)
    assert best["ctgA|gene2"].score_type == "full"


def test_kofam_rejects_evalue_above_limit(monkeypatch):
    _patch_hmmfile(monkeypatch)
    _patch_sequencefile(monkeypatch)
    _patch_hmmsearch(
        monkeypatch,
        [
            FakeTopHits(
                query=FakeHMM("K00003", "K00003"),
                hits=[FakeHit("ctgA|gene1", score=999.0, evalue=1e-2, domains=[FakeDomain(999.0)])],
            )
        ],
    )

    cfg = ann.AnnotationConfig(max_evalue=1e-5)
    kofam_meta = {
        "K00003": KOFamMetadata(
            ko_id="K00003",
            threshold=1.0,
            score_type="full",
            profile_type="all families",
            definition="c",
        )
    }

    best, passing, _, _ = ann._search_and_collect(
        database="kofam",
        hmm_path=Path("kofam.hmm"),
        proteins_path=Path("proteins.faa"),
        cfg=cfg,
        kofam_meta_by_model=kofam_meta,
    )

    assert passing == 0
    assert best == {}


def test_all_hits_output_written_with_header(tmp_path, monkeypatch):
    _patch_hmmfile(monkeypatch)
    _patch_sequencefile(monkeypatch)
    _patch_hmmsearch(
        monkeypatch,
        [
            FakeTopHits(
                query=FakeHMM("A", "PF00021.1", cutoffs=FakeCutoffs(True, ga1=20.0, ga2=10.0)),
                hits=[FakeHit("ctgA|gene1", score=30.0, evalue=1e-9, domains=[FakeDomain(12.0)])],
            )
        ],
    )

    out = tmp_path / "all_hits.tsv"
    cfg = ann.AnnotationConfig(max_evalue=1e-5, keep_all_hits=True, all_hits_path=out)
    with pytest.raises(ValueError):
        ann.annotate_proteins_complete_databases(tmp_path / "proteins.faa", ann.AnnotationConfig(keep_all_hits=True))

    handle, writer = ann._open_all_hits_writer(out)
    try:
        best, passing, _, _ = ann._search_and_collect(
            database="pfam",
            hmm_path=Path("pfam.hmm"),
            proteins_path=Path("proteins.faa"),
            cfg=cfg,
            all_hits_writer=writer,
        )
    finally:
        handle.close()

    text = out.read_text()
    lines = [line for line in text.splitlines() if line.strip()]
    assert lines[0].startswith("protein_id\tcontig_id\tdatabase")
    assert "ctgA|gene1" in lines[1]
    assert passing == 1
    assert "ctgA|gene1" in best