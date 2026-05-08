from pathlib import Path

import phu.screen as screen
from phu.kofam_db import KOFamMetadata


def test_choose_best_contigs_uses_kofam_full_thresholds():
    hits = [
        screen.Hit(
            contig="c1",
            prot_id="c1|gene1",
            model="K00001",
            bitscore=49.0,
            domain_bitscore=80.0,
            evalue=1e-20,
        ),
        screen.Hit(
            contig="c2",
            prot_id="c2|gene1",
            model="K00001",
            bitscore=61.0,
            domain_bitscore=10.0,
            evalue=1e-20,
        ),
    ]

    kept_hits, contigs = screen._choose_best_contigs(
        hits,
        min_bitscore=None,
        max_evalue=1e-5,
        kofam_metadata_by_model={
            "K00001": KOFamMetadata(
                ko_id="K00001",
                threshold=60.0,
                score_type="full",
                profile_type="trim",
                definition="demo",
            )
        },
        use_kofam_thresholds=True,
    )

    assert len(kept_hits) == 1
    assert contigs == ["c2"]


def test_choose_best_contigs_uses_kofam_domain_thresholds():
    hits = [
        screen.Hit(
            contig="c1",
            prot_id="c1|gene1",
            model="K00002",
            bitscore=200.0,
            domain_bitscore=40.0,
            evalue=1e-20,
        ),
        screen.Hit(
            contig="c2",
            prot_id="c2|gene1",
            model="K00002",
            bitscore=20.0,
            domain_bitscore=65.0,
            evalue=1e-20,
        ),
    ]

    kept_hits, contigs = screen._choose_best_contigs(
        hits,
        min_bitscore=None,
        max_evalue=1e-5,
        kofam_metadata_by_model={
            "K00002": KOFamMetadata(
                ko_id="K00002",
                threshold=50.0,
                score_type="domain",
                profile_type="trim",
                definition="demo",
            )
        },
        use_kofam_thresholds=True,
    )

    assert len(kept_hits) == 1
    assert contigs == ["c2"]


def test_extract_target_proteins_adds_ko_definition_to_headers(tmp_path: Path, monkeypatch):
    proteins_fa = tmp_path / "proteins.faa"
    proteins_fa.write_text(
        ">contigA|gene1\nMPEPTIDE\n"
        ">contigB|gene1\nMSEQUENCE\n"
    )

    class _FakeCompleted:
        def __init__(self, returncode: int):
            self.returncode = returncode

    def _fake_seqkit_run(cmd, stdout, text=True):
        ids_file = Path(cmd[3])
        wanted = set(ids_file.read_text().split())

        records = proteins_fa.read_text().strip().split("\n>")
        for i, raw in enumerate(records):
            record = raw if i == 0 else ">" + raw
            lines = record.splitlines()
            header = lines[0].lstrip(">")
            if header in wanted:
                stdout.write(record + "\n")

        return _FakeCompleted(returncode=0)

    monkeypatch.setattr(screen.subprocess, "run", _fake_seqkit_run)

    hit = screen.Hit(
        contig="contigA",
        prot_id="contigA|gene1",
        model="K00001",
        bitscore=100.0,
        domain_bitscore=100.0,
        evalue=1e-20,
    )

    screen._extract_target_proteins(
        kept_hits=[hit],
        kept_contig_ids=["contigA"],
        proteins_fa=proteins_fa,
        outdir=tmp_path,
        hmm_mode="pure",
        seqkit_bin="seqkit",
        kofam_metadata_by_model={
            "K00001": KOFamMetadata(
                ko_id="K00001",
                threshold=50.0,
                score_type="full",
                profile_type="trim",
                definition="hexokinase",
            )
        },
    )

    out = tmp_path / "target_proteins" / "K00001_proteins.mfa"
    content = out.read_text()
    assert ">contigA|gene1|K00001|hexokinase" in content


def test_hmmsearch_retries_without_gathering_cutoffs_on_missing_cutoffs(
    tmp_path: Path, monkeypatch
):
    hmm_path = tmp_path / "K22014.hmm"
    hmm_path.write_text("HMMER3/f\nNAME  K22014\n//\n")

    proteins_fa = tmp_path / "proteins.faa"
    proteins_fa.write_text(">contig1|gene1\nMPEPTIDE\n")

    class _FakeHMMFile:
        def __init__(self, _path):
            self._hmm = object()

        def __enter__(self):
            return [self._hmm]

        def __exit__(self, exc_type, exc, tb):
            return False

    class _FakeSequenceFile:
        def __init__(self, _path, digital=True):
            self._proteins = object()

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def read_block(self):
            return self._proteins

    class _FakeMissingCutoffs(Exception):
        pass

    calls = []

    def _fake_hmmsearch(hmms, proteins, cpus, bit_cutoffs):
        calls.append(bit_cutoffs)
        if bit_cutoffs == "gathering":
            raise _FakeMissingCutoffs("missing gathering cutoffs")
        return []

    monkeypatch.setattr(screen.pyhmmer.plan7, "HMMFile", _FakeHMMFile)
    monkeypatch.setattr(screen.pyhmmer.easel, "SequenceFile", _FakeSequenceFile)
    monkeypatch.setattr(screen.pyhmmer.plan7, "MissingCutoffs", _FakeMissingCutoffs)
    monkeypatch.setattr(screen.pyhmmer, "hmmsearch", _fake_hmmsearch)

    hits = list(
        screen._hmmsearch(
            hmm_paths=[hmm_path],
            proteins_fa=proteins_fa,
            domtbl_paths={},
            threads=1,
            hmm_mode="pure",
            keep_domtbl=False,
            cut_ga=True,
        )
    )

    assert hits == []
    assert calls == ["gathering", None]
