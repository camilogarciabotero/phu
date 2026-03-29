from pathlib import Path

import phu.jack as jack


class _FakeHit:
    def __init__(self, name: str, included: bool, score: float, evalue: float):
        self.name = name
        self.included = included
        self.score = score
        self.evalue = evalue


class _FakeSequenceFile:
    def __init__(self, *_args, **_kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def read_block(self):
        return []


def test_run_jackhmmer_accepts_iteration_objects(monkeypatch, tmp_path):
    class _Iteration:
        def __init__(self, iteration, converged, hits):
            self.iteration = iteration
            self.converged = converged
            self.hits = hits

    hit1 = _FakeHit("contigA|gene1", True, 50.0, 1e-20)
    hit2 = _FakeHit("contigA|gene2", False, 10.0, 1e-2)
    iters = [_Iteration(1, False, [hit1, hit2])]

    monkeypatch.setattr(jack.pyhmmer.easel, "SequenceFile", _FakeSequenceFile)
    monkeypatch.setattr(jack.pyhmmer.hmmer, "jackhmmer", lambda *a, **k: iters)

    kept, summary = jack._run_jackhmmer(
        query=object(),
        alphabet=object(),
        proteins_fa=tmp_path / "proteins.faa",
        iterations=2,
        inc_evalue=1e-3,
        max_evalue=1e-5,
        threads=1,
    )

    assert len(kept) == 1
    assert kept[0].contig == "contigA"
    assert summary[0]["iteration"] == 1
    assert summary[0]["n_hits"] == 2
    assert summary[0]["n_included"] == 1


def test_run_jackhmmer_accepts_list_iterations(monkeypatch, tmp_path):
    # Simulates pyhmmer versions where jackhmmer yields list-like iterations.
    iter1 = [
        _FakeHit("contigB|gene1", True, 42.0, 1e-12),
        _FakeHit("contigC|gene1", True, 41.0, 1e-8),
    ]
    iter2 = [
        _FakeHit("contigB|gene2", True, 60.0, 1e-40),
        _FakeHit("contigC|gene3", False, 9.0, 1e-1),
    ]

    monkeypatch.setattr(jack.pyhmmer.easel, "SequenceFile", _FakeSequenceFile)
    monkeypatch.setattr(jack.pyhmmer.hmmer, "jackhmmer", lambda *a, **k: [iter1, iter2])

    kept, summary = jack._run_jackhmmer(
        query=object(),
        alphabet=object(),
        proteins_fa=tmp_path / "proteins.faa",
        iterations=3,
        inc_evalue=1e-3,
        max_evalue=1e-5,
        threads=1,
    )

    assert len(kept) == 1
    assert kept[0].contig == "contigB"
    assert summary[0]["iteration"] == 1
    assert summary[1]["iteration"] == 2
    assert summary[0]["n_hits"] == 2
    assert summary[1]["n_included"] == 1
