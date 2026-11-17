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


def _model(name, sequence, gathering=None):
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    seq = pyhmmer.easel.TextSequence(name=name.encode(), sequence=sequence)
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(
        seq.digitize(alphabet), pyhmmer.plan7.Background(alphabet)
    )
    hmm.name = name.encode()
    if gathering is not None:
        hmm.cutoffs.gathering = (gathering, gathering)
    return hmm


def test_gathering_cutoffs_survive_models_without_them(tmp_path, capsys):
    # A Pfam-like model with a gathering cutoff too high for the protein, and
    # a model without cutoffs (own HMM, KOfam, PHROG, ECOD) searched with it.
    # The first must stay filtered by its cutoff: before, the missing cutoff
    # of the second made the whole search run without any.
    protein = "MKKLLAAVAGWAAAPAAACDEFGHIKLMNPQRSTVWYMKKLLAAVAGWAAAPAAACDEFGHIK"
    paths = []
    for name, gathering in (("pfam_like", 1000.0), ("no_cutoffs", None)):
        path = tmp_path / f"{name}.hmm"
        with path.open("wb") as handle:
            _model(name, protein, gathering).write(handle)
        paths.append(path)
    proteins_fa = tmp_path / "proteins.faa"
    proteins_fa.write_text(f">contig1|gene1\n{protein}\n")

    def search(cut_ga):
        hits = screen._hmmsearch(
            hmm_paths=paths,
            proteins_fa=proteins_fa,
            domtbl_paths={},
            threads=1,
            hmm_mode="pure",
            keep_domtbl=False,
            cut_ga=cut_ga,
        )
        return sorted(hit.model for hit in hits)

    assert search(cut_ga=True) == ["no_cutoffs"]
    assert search(cut_ga=False) == ["no_cutoffs", "pfam_like"]
    assert "retrying without --cut-ga" not in capsys.readouterr().out


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


def test_phrog_queries_use_the_avger_score_and_coverage_rule():
    def hit(model, bitscore, hmm_coverage, protein_coverage):
        return screen.Hit(
            contig=f"{model}-{bitscore}-{hmm_coverage}-{protein_coverage}",
            prot_id="p|gene1",
            model=model,
            bitscore=bitscore,
            domain_bitscore=bitscore,
            evalue=1e-20,
            hmm_coverage=hmm_coverage,
            protein_coverage=protein_coverage,
        )

    hits = [
        hit("phrog_7", 60.0, 0.5, 0.5),  # exactly at every limit
        hit("phrog_7", 59.9, 1.0, 1.0),
        hit("phrog_7", 300.0, 0.49, 1.0),
        hit("phrog_7", 300.0, 1.0, 0.49),
        hit("phrog_7", 300.0, None, None),
        hit("phrog_9", 20.0, 0.9, 0.1),  # a user's own model that happens to be named like one
    ]

    def kept(**kwargs):
        return [
            (item.model, item.bitscore)
            for item in screen._filter_qualifying_hits(
                hits, None, 1e-6, database_model_ids=("phrog_7",), min_hmm_coverage=0.35, **kwargs
            )
        ]

    assert kept() == [("phrog_7", 60.0), ("phrog_9", 20.0)]
    assert len(kept(use_phrog_thresholds=False)) == 6


def test_extract_target_proteins_adds_phrog_annotation_to_headers(
    tmp_path: Path, monkeypatch
):
    from phu.db.phrog import PhrogAnnotation

    proteins_fa = tmp_path / "proteins.faa"
    proteins_fa.write_text(">contigA|gene1\nMPEPTIDE\n")

    monkeypatch.setattr(
        screen,
        "get_all_phrog_annotations",
        lambda: {"phrog_7": PhrogAnnotation("phrog_7", "anti-restriction", "moron")},
    )

    screen._extract_target_proteins(
        kept_hits=[
            screen.Hit(
                contig="contigA",
                prot_id="contigA|gene1",
                model="phrog_7",
                bitscore=100.0,
                domain_bitscore=100.0,
                evalue=1e-20,
            )
        ],
        kept_contig_ids=["contigA"],
        proteins_fa=proteins_fa,
        outdir=tmp_path,
        hmm_mode="pure",
    )

    content = (tmp_path / "target-proteins" / "phrog_7.faa").read_text()
    assert ">contigA|gene1|phrog_7|anti-restriction|moron" in content


def test_ecod_queries_need_half_of_the_model_aligned():
    def hit(model, hmm_coverage, protein_coverage):
        return screen.Hit(
            contig=f"{model}-{hmm_coverage}-{protein_coverage}",
            prot_id="p|gene1",
            model=model,
            bitscore=25.0,
            domain_bitscore=25.0,
            evalue=1e-20,
            hmm_coverage=hmm_coverage,
            protein_coverage=protein_coverage,
        )

    hits = [
        hit("ecod:1.1.1.2", 0.5, 0.1),  # at the limit; a domain of a long protein
        hit("ecod:1.1.1.2", 0.49, 1.0),
        hit("ecod:1.1.1.2", None, None),
    ]

    def kept(**kwargs):
        return [
            (item.model, item.hmm_coverage)
            for item in screen._filter_qualifying_hits(
                hits, None, 1e-6, database_model_ids=("ecod:1.1.1.2",), **kwargs
            )
        ]

    assert kept() == [("ecod:1.1.1.2", 0.5)]
    assert len(kept(use_ecod_thresholds=False)) == 3


def test_ecod_ids_parse_at_every_level():
    from phu.db.ecod import ecod_level, is_ecod_id, normalize_ecod_id

    assert [is_ecod_id(token) for token in ("ecod:2004", "ECOD:2004.1.1.5", "2004.1.1.5")] == [
        True,
        True,
        False,
    ]
    assert not is_ecod_id("ecod:2004.1.1.5.3")
    assert not is_ecod_id("ecod:2004.x")
    assert normalize_ecod_id(" ECOD:02004.01.1 ") == "ecod:2004.1.1"
    assert [ecod_level(f"ecod:{i}") for i in ("1", "1.1", "1.1.1", "1.1.1.2")] == ["X", "H", "T", "F"]
