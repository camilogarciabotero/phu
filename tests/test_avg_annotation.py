from pathlib import Path

import pytest

from phu.avger.tracks import annotate_avg_tracks
from phu.avger.annotation import AnnotationConfig, AnnotationResults


def empty_results():
    return AnnotationResults({}, {}, 0, 0, 0)


def test_annotation_config_validates_track():
    assert AnnotationConfig(track="relaxed").track == "relaxed"
    assert AnnotationConfig(kofam_require_thresholds=False).kofam_require_thresholds is False
    with pytest.raises(ValueError, match="track"):
        AnnotationConfig(track="other")


def test_avg_annotation_searches_once_for_both_tracks(monkeypatch, tmp_path: Path):
    calls = []

    def fake_annotate(proteins_path, configs, *, phrog):
        assert phrog is True
        calls.append((proteins_path, list(configs)))
        return [empty_results() for _ in configs]

    monkeypatch.setattr("phu.avger.tracks.annotate_proteins_tracks", fake_annotate)

    result = annotate_avg_tracks(tmp_path / "proteins.faa", threads=4)

    assert result.relaxed.passing_hit_count == 0
    assert len(calls) == 1
    relaxed, strict = calls[0][1]
    assert relaxed.track == "relaxed"
    assert relaxed.pfam_require_ga is False
    assert relaxed.kofam_require_thresholds is False
    assert strict.track == "strict"
    assert strict.pfam_require_ga is True
    assert strict.kofam_require_thresholds is True
    assert relaxed.threads == strict.threads == 4


def test_avg_annotation_requires_both_raw_hit_paths():
    with pytest.raises(ValueError, match="both all-hit paths"):
        annotate_avg_tracks(
            Path("proteins.faa"),
            keep_all_hits=True,
            relaxed_all_hits_path=Path("relaxed.tsv"),
        )
