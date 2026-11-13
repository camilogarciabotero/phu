from pathlib import Path

import pytest

import phu.commands.screen as screen
from phu.db.kofam import KOFamMetadata


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


def test_all_requires_queried_models_with_zero_hits():
    hit = screen.Hit(
        contig="c1",
        prot_id="c1|gene1",
        model="model-a",
        bitscore=100.0,
        evalue=1e-20,
    )

    kept_hits, contigs = screen._choose_best_contigs(
        [hit],
        min_bitscore=None,
        max_evalue=1e-5,
        combine_mode="all",
        queried_model_ids=["model-a", "model-b"],
    )

    assert kept_hits == []
    assert contigs == []


def test_threshold_counts_distinct_queried_models():
    hits = [
        screen.Hit(
            contig="c1",
            prot_id="c1|gene1",
            model="model-a",
            bitscore=100.0,
            evalue=1e-20,
        ),
        screen.Hit(
            contig="c1",
            prot_id="c1|gene2",
            model="model-a",
            bitscore=90.0,
            evalue=1e-20,
        ),
    ]

    kept_hits, contigs = screen._choose_best_contigs(
        hits,
        min_bitscore=None,
        max_evalue=1e-5,
        combine_mode="threshold",
        min_hmm_hits=2,
        queried_model_ids=["model-a", "model-b"],
    )

    assert kept_hits == []
    assert contigs == []


@pytest.mark.parametrize("model_count", [0, 2])
def test_query_model_ids_rejects_non_single_pure_files(
    tmp_path: Path, monkeypatch, model_count: int
):
    hmm_path = tmp_path / "models.hmm"

    class FakeModel:
        def __init__(self, name: str):
            self.name = name.encode()

    class FakeHMMFile:
        def __init__(self, _path):
            self.models = [FakeModel(f"model-{i}") for i in range(model_count)]

        def __enter__(self):
            return iter(self.models)

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr(screen.pyhmmer.plan7, "HMMFile", FakeHMMFile)

    with pytest.raises(ValueError, match="exactly one model"):
        screen._query_model_ids([hmm_path], "pure")


def test_generated_hmm_does_not_store_command_line(tmp_path: Path):
    alphabet = screen.pyhmmer.easel.Alphabet.amino()
    builder = screen.pyhmmer.plan7.Builder(alphabet)
    background = screen.pyhmmer.plan7.Background(alphabet)
    sequence = screen.pyhmmer.easel.TextSequence(name=b"protein", sequence="MPEPTIDE")
    hmm, _, _ = builder.build(sequence.digitize(alphabet), background)

    hmm.command_line = None
    output = tmp_path / "model.hmm"
    with output.open("wb") as handle:
        hmm.write(handle)

    with screen.pyhmmer.plan7.HMMFile(output) as handle:
        written_hmm = next(iter(handle))
    assert written_hmm.command_line is None


def test_extract_target_proteins_adds_ko_definition_to_headers(
    tmp_path: Path, monkeypatch
):
    proteins_fa = tmp_path / "proteins.faa"
    proteins_fa.write_text(">contigA|gene1\nMPEPTIDE\n>contigB|gene1\nMSEQUENCE\n")

    class _FakeCompleted:
        def __init__(self, returncode: int):
            self.returncode = returncode

    def _fake_seqkit_run(cmd, stdout, text=True, check=False):
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

    out = tmp_path / "target-proteins" / "K00001.faa"
    content = out.read_text()
    assert ">contigA|gene1|K00001|hexokinase" in content


def test_hmmsearch_retries_without_gathering_cutoffs_on_missing_cutoffs(
    tmp_path: Path, monkeypatch, capsys
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
        def __init__(self, _path, digital=True, alphabet=None):
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

    def _fake_hmmsearch(hmms, proteins, cpus, parallel, bit_cutoffs):
        assert parallel == "queries"
        calls.append(bit_cutoffs)
        if bit_cutoffs == "gathering":
            raise _FakeMissingCutoffs("missing gathering cutoffs")
        return []

    monkeypatch.setattr(screen.pyhmmer.plan7, "HMMFile", _FakeHMMFile)
    monkeypatch.setattr(screen.pyhmmer.easel, "SequenceFile", _FakeSequenceFile)
    monkeypatch.setattr(screen.pyhmmer.plan7, "MissingCutoffs", _FakeMissingCutoffs)
    monkeypatch.setattr(screen.pyhmmer, "hmmsearch", _fake_hmmsearch)

    from phu.core import log as _log

    monkeypatch.setattr(_log, "_verbose", True)
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
    # The verbose log records the search that actually produced the hits.
    retry_lines = [
        line for line in capsys.readouterr().out.splitlines() if line.startswith("[cmd]")
    ]
    assert len(retry_lines) == 1
    assert "pyhmmer.hmmsearch" in retry_lines[0]
    assert "bit_cutoffs=none" in retry_lines[0]


def _coverage_hit(model, coverage, evalue=1e-20):
    return screen.Hit(
        contig=f"{model}-{coverage}",
        prot_id=f"{model}-{coverage}|gene1",
        model=model,
        bitscore=100.0,
        domain_bitscore=100.0,
        evalue=evalue,
        domain_i_evalue=1e-30,
        hmm_coverage=coverage,
    )


def test_user_models_need_more_than_the_minimum_hmm_coverage():
    hits = [
        _coverage_hit("my_family", 0.36),
        _coverage_hit("my_family", 0.35),  # the boundary is excluded
        _coverage_hit("my_family", None),  # no usable alignment
        _coverage_hit("PF00001", 0.10),  # Pfam keeps its own thresholds
        _coverage_hit("K00001", 0.10),  # so does KOfam
    ]
    ko = {"K00001": KOFamMetadata("K00001", 50.0, "full", "trim", "demo")}

    def kept(**kwargs):
        return [
            (hit.model, hit.hmm_coverage)
            for hit in screen._filter_qualifying_hits(
                hits, None, 1e-6, ko, True, database_model_ids=("PF00001", "K00001"), **kwargs
            )
        ]

    assert kept(min_hmm_coverage=0.35) == [("my_family", 0.36), ("PF00001", 0.1), ("K00001", 0.1)]
    assert len(kept(min_hmm_coverage=None)) == 5


def test_screen_defaults_to_evalue_1e6_and_user_coverage_035(tmp_path: Path):
    contigs = tmp_path / "contigs.fna"
    contigs.write_text(">c\nATG\n")
    hmm = tmp_path / "my_family.hmm"
    hmm.write_text("")

    cfg = screen.ScreenConfig(input_contigs=contigs, hmms=[hmm])

    assert (cfg.max_evalue, cfg.min_hmm_coverage) == (1e-6, 0.35)
    plan = cfg.plan()
    assert (plan.max_evalue, plan.min_hmm_coverage) == (1e-6, 0.35)


def test_min_protein_coverage_is_off_by_default_and_applies_to_every_model():
    def hit(model, protein_coverage):
        return screen.Hit(
            contig=f"{model}-{protein_coverage}",
            prot_id=f"{model}-{protein_coverage}|gene1",
            model=model,
            bitscore=100.0,
            domain_bitscore=100.0,
            evalue=1e-20,
            hmm_coverage=0.9,
            protein_coverage=protein_coverage,
        )

    hits = [hit("my_family", 0.5), hit("my_family", 0.49), hit("my_family", None), hit("PF00001", 0.2)]

    def kept(**kwargs):
        return [
            (item.model, item.protein_coverage)
            for item in screen._filter_qualifying_hits(
                hits, None, 1e-6, database_model_ids=("PF00001",), **kwargs
            )
        ]

    assert len(kept()) == 4
    assert kept(min_protein_coverage=0.5) == [("my_family", 0.5)]
