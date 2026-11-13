import json
from pathlib import Path

import pytest

import phu.commands.avger as avg
from phu.avger.tracks import AvgAnnotationResults
from phu.avger.reference import AvgReference
from phu.avger.annotation import AnnotationHit, AnnotationResults
from phu.core.gene_prediction import CacheArtifact, PredictedGene


def test_avg_config_defaults_and_validation(tmp_path: Path):
    input_contigs = tmp_path / "contigs.fna"
    input_contigs.write_text(">contig-1\nATG\n")
    config = avg.AvgConfig(input_contigs=input_contigs)

    assert config.output_folder == Path("phu-avger")
    assert config.mode == "meta"
    assert config.translation_table is None
    assert config.min_amg_weight == 0.6
    assert config.filter_mode == "standard"
    assert config.scoring_evalue == 1e-5

    with pytest.raises(ValueError, match="min_amg_weight"):
        avg.AvgConfig(input_contigs=config.input_contigs, min_amg_weight=1.1)
    with pytest.raises(ValueError, match="filter_mode"):
        avg.AvgConfig(input_contigs=config.input_contigs, filter_mode="bad")


def test_require_databases_resolves_every_database_avger_uses(monkeypatch):
    resolved = []
    monkeypatch.setattr(avg, "resolve_db", resolved.append)

    avg.require_databases()
    avg.require_databases(phrog=False)

    assert resolved == ["pfam", "kofam", "avg", "phrog", "pfam", "kofam", "avg"]


def test_run_avg_reuses_prediction_cache_and_passes_track_paths(monkeypatch, tmp_path):
    input_contigs = tmp_path / "contigs.fna"
    input_contigs.write_text(">contig-1\nATG\n")
    prediction = CacheArtifact(
        proteins_path=tmp_path / "proteins.faa",
        protein_count=1,
        cache_hit=True,
        cache_key="cache-key",
    )
    calls = {}
    monkeypatch.setattr(avg, "require_databases", lambda phrog: None)
    monkeypatch.setattr(avg, "get_or_predict_proteins", lambda *args, **kwargs: prediction)
    monkeypatch.setattr(avg, "write_avg_outputs", lambda *args, **kwargs: ())
    monkeypatch.setattr(
        avg,
        "annotate_avg_tracks",
        lambda *args, **kwargs: calls.update(kwargs) or AvgAnnotationResults(
            AnnotationResults({}, {}, 0, 0, 0), AnnotationResults({}, {}, 0, 0, 0)
        ),
    )

    result = avg.run_avg(avg.AvgConfig(input_contigs=input_contigs, keep_hits=True))

    assert result.prediction.cache_hit is True
    assert calls["relaxed_all_hits_path"].name == "hmm-hits-relaxed.tsv.gz"
    assert calls["strict_all_hits_path"].name == "hmm-hits-strict.tsv.gz"


def test_avg_outputs_preserve_pfam_and_kofam_function_metadata(monkeypatch, tmp_path):
    input_contigs = tmp_path / "contigs.fna"
    input_contigs.write_text(">contig-1\nATG\n")
    gene = PredictedGene("contig-1", "contig-1|gene1", 1, 30, 1, 1, "ATG", "M")

    def hit(database, model_id, name, description):
        return AnnotationHit(
            protein_id=gene.gene_id,
            contig_id=gene.contig_id,
            database=database,
            model_id=model_id,
            model_accession=model_id,
            score_type="full",
            effective_score=50,
            full_score=50,
            domain_score=None,
            evalue=1e-10,
            hmm_from=None,
            hmm_to=None,
            target_from=None,
            target_to=None,
            threshold_source="test",
            threshold_value=1,
            model_name=name,
            model_description=description,
        )

    relaxed_hit = hit("pfam", "PF00001", "PF test", "Pfam function")
    strict_hit = hit("kofam", "K00001", "K00001", "KOfam function")
    annotations = AvgAnnotationResults(
        AnnotationResults({relaxed_hit.protein_id: relaxed_hit}, {}, 1, 1, 0),
        AnnotationResults({}, {strict_hit.protein_id: strict_hit}, 1, 1, 0),
    )
    prediction = CacheArtifact(tmp_path / "proteins.faa", 1, True, "key", genes=[gene])
    provenance = {"source": "local", "bundle_version": "phu-db-test"}
    monkeypatch.setattr(avg, "resolve_db", lambda name: tmp_path / name)
    monkeypatch.setattr(avg, "bundle_provenance", lambda: provenance)
    monkeypatch.setattr(
        avg,
        "load_avg_reference",
        lambda root: AvgReference(
            scores=[{"database": "pfam", "accession": "PF00001", "v_score": "9", "vl_score": "2"}],
            positive=[],
            filters=[],
        ),
    )

    avg.write_avg_outputs(avg.AvgConfig(input_contigs=input_contigs, output_folder=tmp_path / "out"), prediction, annotations)

    header, row = (tmp_path / "out" / "gene-annotations.tsv").read_text().splitlines()
    assert "pfam_name" in header
    assert "pfam_description" in header
    assert "kofam_name" in header
    assert "kofam_description" in header
    assert "K00001" in row
    assert "KOfam function" in row
    evidence = (tmp_path / "out" / "avg-evidence.tsv").read_text()
    assert "PF test" in evidence
    assert "Pfam function" in evidence
    assert not (tmp_path / "out" / "run.json").exists()
    run = json.loads((tmp_path / "out" / ".phu" / "run.json").read_text())
    assert run["phu_db"] == provenance
    assert "avg_reference" not in run


def _phrog_outputs(monkeypatch, tmp_path, strict_ids):
    """One contig, two genes; gene1 is a candidate through its PHROG scores."""
    input_contigs = tmp_path / "contigs.fna"
    input_contigs.write_text(">contig-1\nATG\n")
    genes = [
        PredictedGene("contig-1", f"contig-1|gene{n}", 1, 30, 1, n, "ATG", "M") for n in (1, 2)
    ]

    def hit(gene, model_id, description=""):
        return AnnotationHit(
            protein_id=gene.gene_id,
            contig_id=gene.contig_id,
            database="phrog",
            model_id=model_id,
            model_accession=model_id,
            score_type="full",
            effective_score=80,
            full_score=80,
            domain_score=80,
            evalue=1e-20,
            hmm_from=1,
            hmm_to=50,
            target_from=1,
            target_to=50,
            threshold_source="phrog_score_coverage",
            threshold_value=60,
            model_name=model_id,
            model_description=description,
            model_coverage=0.5,
            protein_coverage=0.75,
        )

    relaxed = {genes[0].gene_id: hit(genes[0], "phrog_2"), genes[1].gene_id: hit(genes[1], "phrog_1")}
    strict = {genes[0].gene_id: hit(genes[0], strict_ids, "anti-restriction protein")}
    annotations = AvgAnnotationResults(
        AnnotationResults({}, {}, 2, 2, 0, best_phrog_by_protein=relaxed),
        AnnotationResults({}, {}, 1, 1, 0, best_phrog_by_protein=strict),
    )
    (tmp_path / "phrog").mkdir()
    (tmp_path / "phrog" / "annotations.tsv").write_text(
        "phrog\tcolor\tannot\tcategory\n1\t#f\tintegrase\tintegration and excision\n"
        "2\t#f\tanti-restriction protein\tmoron\n"
    )
    monkeypatch.setattr(avg, "resolve_db", lambda name: tmp_path / name)
    monkeypatch.setattr(avg, "bundle_provenance", lambda: {})
    monkeypatch.setattr(
        avg,
        "load_avg_reference",
        lambda root: AvgReference(
            # Contig AVL = (1.0 + 6.0) / 2 = 3.5; gene1 has VL 1.0 and V 0.8.
            scores=[
                {"database": "phrog", "accession": "phrog_2", "v_score": "0.8", "vl_score": "1.0"},
                {"database": "phrog", "accession": "phrog_1", "v_score": "10", "vl_score": "6.0"},
            ],
            positive=[
                {"database": "phrog", "accession": "phrog_2", "proposed_class": "putative_apg", "name": "x"}
            ],
            filters=[
                {
                    "database": "phrog",
                    "accession": "phrog_1",
                    "proposed_class": "putative_apg",
                    "filter_category": "filter_essential",
                }
            ],
        ),
    )
    prediction = CacheArtifact(tmp_path / "proteins.faa", 2, True, "key", genes=genes)
    out = tmp_path / "out"
    avg.write_avg_outputs(avg.AvgConfig(input_contigs=input_contigs, output_folder=out), prediction, annotations)
    lines = (out / "gene-annotations.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]], out


def test_phrog_alone_can_make_and_classify_a_candidate(monkeypatch, tmp_path):
    rows, out = _phrog_outputs(monkeypatch, tmp_path, "phrog_2")

    gene1, gene2 = rows
    assert gene1["zhou_phrog_candidate"] == "true"
    assert gene1["zhou_avg_candidate"] == "true"
    assert gene1["phrog_scaffold_avl"] == "3.5"
    assert gene1["phrog_scored_gene_count"] == "2"
    assert (gene1["phrog_id"], gene1["phrog_category"]) == ("phrog_2", "moron")
    assert (gene1["phrog_model_coverage"], gene1["phrog_protein_coverage"]) == ("0.5", "0.75")
    assert (gene1["classification_status"], gene1["avg_class"]) == ("classified", "putative_apg")
    assert gene2["zhou_phrog_candidate"] == "false"
    assert gene2["phrog_id"] == ""
    assert len((out / "avg-predictions.tsv").read_text().splitlines()) == 2
    run = json.loads((out / ".phu" / "run.json").read_text())
    assert run["parameters"]["phrog"] is True
    assert run["output_schema_version"] == 2


def test_phrog_columns_stay_in_the_tables_when_phrog_is_off(monkeypatch, tmp_path):
    input_contigs = tmp_path / "contigs.fna"
    input_contigs.write_text(">contig-1\nATG\n")
    gene = PredictedGene("contig-1", "contig-1|gene1", 1, 30, 1, 1, "ATG", "M")
    monkeypatch.setattr(avg, "resolve_db", lambda name: tmp_path / name)
    monkeypatch.setattr(avg, "bundle_provenance", lambda: {})
    monkeypatch.setattr(avg, "load_avg_reference", lambda root: AvgReference([], [], []))
    empty = AnnotationResults({}, {}, 0, 0, 0)

    avg.write_avg_outputs(
        avg.AvgConfig(input_contigs=input_contigs, output_folder=tmp_path / "out", phrog=False),
        CacheArtifact(tmp_path / "proteins.faa", 1, True, "key", genes=[gene]),
        AvgAnnotationResults(empty, empty),
    )

    header = (tmp_path / "out" / "gene-annotations.tsv").read_text().splitlines()[0].split("\t")
    assert {"phrog_id", "phrog_v_score", "phrog_scaffold_avl", "zhou_phrog_candidate"} <= set(header)


@pytest.mark.parametrize(
    ("phrog_hit", "status", "avg_class"),
    [(None, "classified", "putative_apg"), ("phrog_1", "filtered", "")],
)
def test_a_phrog_filter_hit_vetoes_a_class_another_database_supports(
    monkeypatch, tmp_path, phrog_hit, status, avg_class
):
    """Pfam alone makes the gene an APG; a PHROG hit on the filter list blocks it."""
    input_contigs = tmp_path / "contigs.fna"
    input_contigs.write_text(">contig-1\nATG\n")
    gene = PredictedGene("contig-1", "contig-1|gene1", 1, 30, 1, 1, "ATG", "M")

    def hit(database, model_id):
        return AnnotationHit(
            protein_id=gene.gene_id,
            contig_id=gene.contig_id,
            database=database,
            model_id=model_id,
            model_accession=model_id,
            score_type="full",
            effective_score=80,
            full_score=80,
            domain_score=80,
            evalue=1e-20,
            hmm_from=1,
            hmm_to=50,
            target_from=1,
            target_to=50,
            threshold_source="test",
            threshold_value=1,
        )

    pfam = {gene.gene_id: hit("pfam", "PF00001")}
    phrog = {} if phrog_hit is None else {gene.gene_id: hit("phrog", phrog_hit)}
    annotations = AvgAnnotationResults(
        AnnotationResults(dict(pfam), {}, 1, 1, 0),
        AnnotationResults(dict(pfam), {}, 1, 1, 0, best_phrog_by_protein=phrog),
    )
    (tmp_path / "phrog").mkdir()
    (tmp_path / "phrog" / "annotations.tsv").write_text(
        "phrog\tcolor\tannot\tcategory\n1\t#f\tmajor head protein\thead and packaging\n"
    )
    monkeypatch.setattr(avg, "resolve_db", lambda name: tmp_path / name)
    monkeypatch.setattr(avg, "bundle_provenance", lambda: {})
    monkeypatch.setattr(
        avg,
        "load_avg_reference",
        lambda root: AvgReference(
            scores=[{"database": "pfam", "accession": "PF00001", "v_score": "1", "vl_score": "1"}],
            positive=[
                {"database": "pfam", "accession": "PF00001", "proposed_class": "putative_apg", "name": "x"}
            ],
            filters=[
                {
                    "database": "phrog",
                    "accession": "phrog_1",
                    "proposed_class": "putative_apg",
                    "filter_category": "filter_essential",
                }
            ],
        ),
    )
    out = tmp_path / "out"
    avg.write_avg_outputs(
        # A lone scored gene has AVL equal to its own VL, so lower the contig cutoff.
        avg.AvgConfig(input_contigs=input_contigs, output_folder=out, scaffold_avl_cutoff=0.5),
        CacheArtifact(tmp_path / "proteins.faa", 1, True, "key", genes=[gene]),
        annotations,
    )

    header, row = (line.split("\t") for line in (out / "gene-annotations.tsv").read_text().splitlines())
    result = dict(zip(header, row))
    assert result["zhou_avg_candidate"] == "true"
    assert (result["classification_status"], result["avg_class"]) == (status, avg_class)
    assert len((out / "avg-predictions.tsv").read_text().splitlines()) == (2 if avg_class else 1)
