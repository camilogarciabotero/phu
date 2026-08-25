from pathlib import Path

import pytest

from phu.avg_annotation import annotate_avg_tracks
from phu.avger_annotation import AnnotationConfig, AnnotationResults


def empty_results():
    return AnnotationResults({}, {}, 0, 0, 0)


def test_annotation_config_validates_track():
    assert AnnotationConfig(track="relaxed").track == "relaxed"
    assert AnnotationConfig(kofam_require_thresholds=False).kofam_require_thresholds is False
    with pytest.raises(ValueError, match="track"):
        AnnotationConfig(track="other")


def test_avg_annotation_runs_distinct_relaxed_and_strict_configs(monkeypatch, tmp_path: Path):
    calls = []

    def fake_annotate(proteins_path, config):
        calls.append((proteins_path, config))
        return empty_results()

    monkeypatch.setattr("phu.avg_annotation.annotate_proteins_complete_databases", fake_annotate)

    result = annotate_avg_tracks(tmp_path / "proteins.faa", threads=4)

    assert result.relaxed.passing_hit_count == 0
    assert len(calls) == 2
    relaxed = calls[0][1]
    strict = calls[1][1]
    assert relaxed.track == "relaxed"
    assert relaxed.pfam_require_ga is False
    assert relaxed.kofam_require_thresholds is False
    assert strict.track == "strict"
    assert strict.pfam_require_ga is True
    assert strict.kofam_require_thresholds is True
    assert all(config.threads == 4 for _, config in calls)


def test_avg_annotation_requires_both_raw_hit_paths():
    with pytest.raises(ValueError, match="both all-hit paths"):
        annotate_avg_tracks(
            Path("proteins.faa"),
            keep_all_hits=True,
            relaxed_all_hits_path=Path("relaxed.tsv"),
        )
