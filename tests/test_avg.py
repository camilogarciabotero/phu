from pathlib import Path

import pytest

import phu.avg as avg
from phu.avg_annotation import AvgAnnotationResults
from phu.avger_annotation import AnnotationHit, AnnotationResults
from phu.gene_prediction_core import CacheArtifact, PredictedGene


def test_avg_config_defaults_and_validation(tmp_path: Path):
    input_contigs = tmp_path / "contigs.fna"
    input_contigs.write_text(">contig-1\nATG\n")
    config = avg.AvgConfig(input_contigs=input_contigs)

    assert config.output_folder == Path("phu-avger")
    assert config.mode == "meta"
    assert config.translation_table == 11
    assert config.min_amg_weight == 0.6
    assert config.filter_mode == "standard"
    assert config.scoring_evalue == 1e-5

    with pytest.raises(ValueError, match="min_amg_weight"):
        avg.AvgConfig(input_contigs=config.input_contigs, min_amg_weight=1.1)
    with pytest.raises(ValueError, match="filter_mode"):
        avg.AvgConfig(input_contigs=config.input_contigs, filter_mode="bad")


def test_database_readiness_and_actionable_missing_message(monkeypatch):
    status = {
        "pfam": {"downloaded": True, "indexed": True},
        "kofam": {"downloaded": False, "indexed": False},
        "avg": {"downloaded": True, "indexed": False},
    }
    monkeypatch.setattr(avg, "get_pfam_database_status", lambda: status["pfam"])
    monkeypatch.setattr(avg, "get_kofam_database_status", lambda: status["kofam"])
    monkeypatch.setattr(avg, "get_avg_database_status", lambda: status["avg"])

    assert avg.missing_databases() == ["kofam", "avg"]
    with pytest.raises(FileNotFoundError, match="phu dbs prepare pfam kofam avg"):
        avg.require_databases()


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
    monkeypatch.setattr(avg, "require_databases", lambda: None)
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
    assert calls["relaxed_all_hits_path"].name == "relaxed_hits.tsv.gz"
    assert calls["strict_all_hits_path"].name == "strict_hits.tsv.gz"


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
    (tmp_path / "manifest.json").write_text("{}")
    monkeypatch.setattr(avg, "get_avg_database_status", lambda: {"root": str(tmp_path)})
    monkeypatch.setattr(avg, "get_pfam_database_status", lambda: {"name": "pfam"})
    monkeypatch.setattr(avg, "get_kofam_database_status", lambda: {"name": "kofam"})
    monkeypatch.setattr(
        avg,
        "_reference_rows",
        lambda path: (
            [{"database": "pfam", "accession": "PF00001", "v_score": "9", "vl_score": "2"}]
            if path.name == "v_scores.tsv"
            else []
        ),
    )

    avg.write_avg_outputs(avg.AvgConfig(input_contigs=input_contigs, output_folder=tmp_path / "out"), prediction, annotations)

    header, row = (tmp_path / "out" / "genes.tsv").read_text().splitlines()
    assert "pfam_name" in header
    assert "pfam_description" in header
    assert "kofam_name" in header
    assert "kofam_description" in header
    assert "K00001" in row
    assert "KOfam function" in row
    evidence = (tmp_path / "out" / "evidence.tsv").read_text()
    assert "PF test" in evidence
    assert "Pfam function" in evidence
    assert not (tmp_path / "out" / "run.json").exists()
    assert (tmp_path / "out" / ".phu" / "run.json").exists()
