from __future__ import annotations

from pathlib import Path

import pytest

import phu.commands.screen as screen
import phu.core.memory as memory
from test_hmm_parallel_reproducibility import (
    _build_hmm,
    _hit_signature,
)


@pytest.fixture
def no_cgroup(monkeypatch):
    monkeypatch.setattr(memory, "_cgroup_room", lambda: None)


def _set_available(monkeypatch, available):
    monkeypatch.setattr(memory, "available_memory", lambda: available)


def test_fits_when_file_is_a_small_share(tmp_path, monkeypatch):
    path = tmp_path / "p.faa"
    path.write_bytes(b"x" * 100)
    _set_available(monkeypatch, 2_000)  # 10% = 200 = 100 bytes x 2 overhead
    assert memory.fits_in_memory(path)
    _set_available(monkeypatch, 1_999)
    assert not memory.fits_in_memory(path)


def test_fits_when_memory_is_unknown(tmp_path, monkeypatch):
    path = tmp_path / "p.faa"
    path.write_bytes(b"x" * 10**6)
    _set_available(monkeypatch, None)
    assert memory.fits_in_memory(path)


def test_cgroup_limit_caps_available_memory(monkeypatch):
    class VM:
        available = 10**12

    monkeypatch.setattr(memory.psutil, "virtual_memory", lambda: VM)
    monkeypatch.setattr(memory, "_cgroup_room", lambda: 4 * 10**9)
    assert memory.available_memory() == 4 * 10**9


def test_cgroup_room_reads_v2_and_ignores_no_limit(tmp_path, monkeypatch):
    limit, usage = tmp_path / "max", tmp_path / "current"
    usage.write_text("1000\n")
    monkeypatch.setattr(memory, "_CGROUP_V2", (str(limit), str(usage)))
    monkeypatch.setattr(memory, "_CGROUP_V1", ("/nonexistent", "/nonexistent"))
    limit.write_text("max\n")
    assert memory._cgroup_room() is None
    limit.write_text("5000\n")
    assert memory._cgroup_room() == 4000


def test_protein_batch_size(tmp_path, monkeypatch):
    path = tmp_path / "p.faa"
    path.write_bytes(b"x" * 100)
    _set_available(monkeypatch, 10**9)
    assert memory.protein_batch_size(path, 50_000, 10_000) == 50_000
    assert memory.protein_batch_size(path, 5_000, 10_000) == 10_000
    _set_available(monkeypatch, 10)
    assert memory.protein_batch_size(path, 50_000, 10_000) == 10_000


@pytest.mark.parametrize("cut_ga", [False, True])
def test_screen_streamed_targets_give_same_hits(tmp_path, monkeypatch, cut_ga):
    """Reading proteins from disk finds what loading them finds, for every model."""
    hmm_path = tmp_path / "two.hmm"
    with hmm_path.open("wb") as handle:
        _build_hmm("model_a", "MKKLLAAVAGAAAAPAAA").write(handle)
        _build_hmm("model_b", "WWHHCCYYWWHHCCYYWW").write(handle)
    proteins = tmp_path / "proteins.faa"
    proteins.write_text(
        ">c1_1\nMKKLLAAVAGAAAAPAAA\n"
        ">c1_2\nWWHHCCYYWWHHCCYYWW\n"
        ">c2_1\nGKKLLAAVAGAAAAPAAA\n"
    )

    def run():
        return list(
            screen._hmmsearch(
                [hmm_path],
                proteins,
                {},
                hmm_mode="mixed",
                keep_domtbl=False,
                cut_ga=cut_ga,
            )
        )

    monkeypatch.setattr(screen, "fits_in_memory", lambda path: True)
    loaded = run()
    monkeypatch.setattr(screen, "fits_in_memory", lambda path: False)
    streamed = run()
    # The second model must still see every protein after the first read the file.
    assert {(hit.model, hit.prot_id) for hit in streamed} == {
        ("model_a", "c1_1"),
        ("model_a", "c2_1"),
        ("model_b", "c1_2"),
    }
    assert _hit_signature(streamed) == _hit_signature(loaded)
