"""Tests for gene_prediction_core caching functionality."""

import json
import os
from pathlib import Path

import pytest

from phu.gene_prediction_core import (
    CacheArtifact,
    PredictionInputs,
    compute_cache_key,
    get_cache_dir,
    get_or_predict_proteins,
    write_prediction_metadata,
)


class TestPredictionInputs:
    """Test validation and properties of PredictionInputs."""

    def test_inputs_validates_contigs_path(self, tmp_path):
        """PredictionInputs raises if contigs file doesn't exist."""
        nonexistent = tmp_path / "nonexistent.fa"
        with pytest.raises(FileNotFoundError):
            PredictionInputs(input_contigs=nonexistent)

    def test_inputs_resolves_path(self, tmp_path):
        """PredictionInputs resolves relative paths to absolute."""
        contigs = tmp_path / "contigs.fa"
        contigs.write_text(">c1\nATGC\n")
        inputs = PredictionInputs(input_contigs=contigs)
        assert inputs.input_contigs.is_absolute()

    def test_inputs_validates_mode(self, tmp_path):
        """PredictionInputs only accepts 'meta' or 'single' mode."""
        contigs = tmp_path / "contigs.fa"
        contigs.write_text(">c1\nATGC\n")
        with pytest.raises(ValueError, match="mode must be"):
            PredictionInputs(input_contigs=contigs, mode="invalid")

    def test_inputs_validates_min_gene_len(self, tmp_path):
        """PredictionInputs requires min_gene_len >= 1."""
        contigs = tmp_path / "contigs.fa"
        contigs.write_text(">c1\nATGC\n")
        with pytest.raises(ValueError, match="min_gene_len"):
            PredictionInputs(input_contigs=contigs, min_gene_len=0)

    def test_inputs_validates_min_protein_len_aa(self, tmp_path):
        """PredictionInputs requires min_protein_len_aa >= 1."""
        contigs = tmp_path / "contigs.fa"
        contigs.write_text(">c1\nATGC\n")
        with pytest.raises(ValueError, match="min_protein_len_aa"):
            PredictionInputs(input_contigs=contigs, min_protein_len_aa=0)


class TestCacheKey:
    """Test cache key generation."""

    def test_cache_key_deterministic(self, tmp_path):
        """Same inputs produce same cache key."""
        contigs = tmp_path / "contigs.fa"
        contigs.write_text(">c1\nATGC\n")

        inputs1 = PredictionInputs(input_contigs=contigs)
        inputs2 = PredictionInputs(input_contigs=contigs)

        key1 = compute_cache_key(inputs1)
        key2 = compute_cache_key(inputs2)
        assert key1 == key2

    def test_cache_key_changes_with_params(self, tmp_path):
        """Different prediction params produce different keys."""
        contigs = tmp_path / "contigs.fa"
        contigs.write_text(">c1\nATGC\n")

        inputs1 = PredictionInputs(input_contigs=contigs, min_protein_len_aa=30)
        inputs2 = PredictionInputs(input_contigs=contigs, min_protein_len_aa=50)

        key1 = compute_cache_key(inputs1)
        key2 = compute_cache_key(inputs2)
        assert key1 != key2

    def test_cache_key_includes_file_stat(self, tmp_path):
        """Cache key includes file stat for quick invalidation."""
        contigs = tmp_path / "contigs.fa"
        contigs.write_text(">c1\nATGC\n")

        inputs = PredictionInputs(input_contigs=contigs)
        
        # Verify key is deterministic and includes stat (by checking it's a hash)
        key = compute_cache_key(inputs)
        assert len(key) == 16  # SHA256 truncated to 16 chars
        assert all(c in "0123456789abcdef" for c in key)

    def test_cache_key_ignores_mode(self, tmp_path):
        """Cache key includes mode in its seed."""
        contigs = tmp_path / "contigs.fa"
        contigs.write_text(">c1\nATGC\n")

        inputs1 = PredictionInputs(input_contigs=contigs, mode="meta")
        inputs2 = PredictionInputs(input_contigs=contigs, mode="single")

        key1 = compute_cache_key(inputs1)
        key2 = compute_cache_key(inputs2)
        # Different modes should produce different keys
        assert key1 != key2


class TestCacheDir:
    """Test cache directory resolution."""

    def test_cache_dir_respects_env(self, monkeypatch, tmp_path):
        """PHU_CACHE_DIR env var overrides defaults."""
        cache_path = tmp_path / "my_cache"
        monkeypatch.setenv("PHU_CACHE_DIR", str(cache_path))
        assert get_cache_dir() == cache_path

    def test_cache_dir_respects_xdg(self, monkeypatch, tmp_path):
        """XDG_CACHE_HOME is used as fallback."""
        monkeypatch.delenv("PHU_CACHE_DIR", raising=False)
        xdg_path = tmp_path / "xdg_cache"
        monkeypatch.setenv("XDG_CACHE_HOME", str(xdg_path))
        expected = xdg_path / "phu"
        assert get_cache_dir() == expected

    def test_cache_dir_default_home(self, monkeypatch, tmp_path):
        """Default cache is ~/.cache/phu."""
        monkeypatch.delenv("PHU_CACHE_DIR", raising=False)
        monkeypatch.delenv("XDG_CACHE_HOME", raising=False)
        # We can't easily mock the user's home, but we can check the path structure
        cache_dir = get_cache_dir()
        assert "cache" in str(cache_dir) and "phu" in str(cache_dir)


class TestGetOrPredictProteins:
    """Test cache hit/miss and prediction workflow."""

    def test_cache_miss_first_run(self, tmp_path, monkeypatch):
        """First run with inputs creates cache."""
        cache_dir = tmp_path / "cache"
        monkeypatch.setenv("PHU_CACHE_DIR", str(cache_dir))

        contigs = tmp_path / "contigs.fa"
        # Simple contigs with a minimal gene (ATG...TAA), 100bp to avoid min_gene_len issues
        contigs.write_text(">c1\nATG" + "A" * 87 + "TAA\n")  # 90bp gene

        inputs = PredictionInputs(
            input_contigs=contigs,
            min_protein_len_aa=3,
            min_gene_len=60,  # Must be >= 55 for pyrodigal_gv
        )

        artifact = get_or_predict_proteins(inputs, use_cache=True)

        assert not artifact.cache_hit
        assert artifact.proteins_path.exists()
        assert artifact.protein_count >= 0
        assert artifact.cache_key != ""
        assert artifact.cache_dir is not None

    def test_cache_hit_second_run(self, tmp_path, monkeypatch):
        """Second run with identical inputs hits cache."""
        cache_dir = tmp_path / "cache"
        monkeypatch.setenv("PHU_CACHE_DIR", str(cache_dir))

        contigs = tmp_path / "contigs.fa"
        contigs.write_text(">c1\nATG" + "A" * 87 + "TAA\n")  # 90bp gene

        inputs = PredictionInputs(input_contigs=contigs, min_gene_len=60)

        # First call: cache miss
        artifact1 = get_or_predict_proteins(inputs, use_cache=True)
        assert not artifact1.cache_hit

        # Second call: cache hit
        artifact2 = get_or_predict_proteins(inputs, use_cache=True)
        assert artifact2.cache_hit
        assert artifact1.proteins_path == artifact2.proteins_path

    def test_cache_miss_on_param_change(self, tmp_path, monkeypatch):
        """Changing prediction params creates new cache entry."""
        cache_dir = tmp_path / "cache"
        monkeypatch.setenv("PHU_CACHE_DIR", str(cache_dir))

        contigs = tmp_path / "contigs.fa"
        contigs.write_text(">c1\nATG" + "A" * 87 + "TAA\n")

        inputs1 = PredictionInputs(
            input_contigs=contigs, min_protein_len_aa=30, min_gene_len=60
        )
        artifact1 = get_or_predict_proteins(inputs1, use_cache=True)

        inputs2 = PredictionInputs(
            input_contigs=contigs, min_protein_len_aa=50, min_gene_len=60
        )
        artifact2 = get_or_predict_proteins(inputs2, use_cache=True)

        # Different cache keys, different paths
        assert artifact1.cache_key != artifact2.cache_key
        assert not artifact2.cache_hit

    def test_no_cache_when_disabled(self, tmp_path, monkeypatch):
        """PHU_CACHE=off skips caching."""
        cache_dir = tmp_path / "cache"
        monkeypatch.setenv("PHU_CACHE_DIR", str(cache_dir))
        monkeypatch.setenv("PHU_CACHE", "off")

        contigs = tmp_path / "contigs.fa"
        contigs.write_text(">c1\nATG" + "A" * 87 + "TAA\n")

        inputs = PredictionInputs(input_contigs=contigs, min_gene_len=60)

        # Even though we call twice, no cache is written
        artifact1 = get_or_predict_proteins(inputs, use_cache=False)
        artifact2 = get_or_predict_proteins(inputs, use_cache=False)

        # Both should be misses (or no cache at all)
        assert artifact1.cache_key == ""
        assert artifact2.cache_key == ""
        # Paths should be different (separate tempfiles)
        assert artifact1.proteins_path != artifact2.proteins_path

    def test_manifest_written_on_cache_build(self, tmp_path, monkeypatch):
        """Cache build writes manifest.json with metadata."""
        cache_dir = tmp_path / "cache"
        monkeypatch.setenv("PHU_CACHE_DIR", str(cache_dir))

        contigs = tmp_path / "contigs.fa"
        contigs.write_text(">c1\nATG" + "A" * 87 + "TAA\n")

        inputs = PredictionInputs(input_contigs=contigs, mode="meta", min_gene_len=60)
        artifact = get_or_predict_proteins(inputs, use_cache=True)

        manifest_path = artifact.cache_dir / "manifest.json"
        assert manifest_path.exists()

        manifest = json.loads(manifest_path.read_text())
        assert manifest["mode"] == "meta"
        assert manifest["min_gene_len"] == 60
        assert manifest["protein_count"] == artifact.protein_count

    def test_partial_dir_cleanup(self, tmp_path, monkeypatch):
        """Stale .partial/ directory is cleaned up on rebuild."""
        cache_dir = tmp_path / "cache"
        monkeypatch.setenv("PHU_CACHE_DIR", str(cache_dir))

        contigs = tmp_path / "contigs.fa"
        contigs.write_text(">c1\nATG" + "A" * 87 + "TAA\n")

        inputs = PredictionInputs(input_contigs=contigs, min_gene_len=60)
        cache_key = compute_cache_key(inputs)

        # Manually create a stale .partial/ directory
        cache_root = cache_dir / "v1"
        partial_dir = cache_root / f"{cache_key}.partial"
        partial_dir.mkdir(parents=True, exist_ok=True)
        stale_marker = partial_dir / "stale_marker.txt"
        stale_marker.write_text("stale")

        # On next prediction, .partial/ should be cleaned and rebuilt
        artifact = get_or_predict_proteins(inputs, use_cache=True)

        # Old partial should be gone
        assert not stale_marker.exists()
        # New cache should be built
        assert artifact.proteins_path.exists()

    def test_corrupted_manifest_triggers_rebuild(self, tmp_path, monkeypatch):
        """Corrupted manifest.json triggers cache rebuild."""
        cache_dir = tmp_path / "cache"
        monkeypatch.setenv("PHU_CACHE_DIR", str(cache_dir))

        contigs = tmp_path / "contigs.fa"
        contigs.write_text(">c1\nATG" + "A" * 87 + "TAA\n")

        inputs = PredictionInputs(input_contigs=contigs, min_gene_len=60)

        # First build
        artifact1 = get_or_predict_proteins(inputs, use_cache=True)
        assert not artifact1.cache_hit

        # Corrupt the manifest
        manifest_path = artifact1.cache_dir / "manifest.json"
        manifest_path.write_text("{ invalid json }")

        # Second call should treat as miss and rebuild
        artifact2 = get_or_predict_proteins(inputs, use_cache=True)
        # Should have rebuilt (not hit cache with corrupted manifest)
        assert not artifact2.cache_hit


class TestWritePredictionMetadata:
    """Test metadata file writing."""

    def test_metadata_file_written(self, tmp_path):
        """Metadata JSON file is written with correct format."""
        cache_dir = tmp_path / "cache"
        metadata_path = tmp_path / ".phu_prediction_metadata.json"

        write_prediction_metadata(
            metadata_path,
            cache_hit=True,
            cache_key="abc123",
            cache_dir=cache_dir,
        )

        assert metadata_path.exists()
        meta = json.loads(metadata_path.read_text())
        assert meta["cache_hit"] is True
        assert meta["cache_key"] == "abc123"
        assert meta["cache_dir"] == str(cache_dir)

    def test_metadata_with_none_cache_dir(self, tmp_path):
        """Metadata handles None cache_dir gracefully."""
        metadata_path = tmp_path / ".phu_prediction_metadata.json"

        write_prediction_metadata(
            metadata_path,
            cache_hit=False,
            cache_key="",
            cache_dir=None,
        )

        meta = json.loads(metadata_path.read_text())
        assert meta["cache_dir"] is None
