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


class _FakeHMM:
    def __init__(self, payload: bytes):
        self.payload = payload

    def write(self, fh):
        fh.write(self.payload)


def test_run_jackhmmer_accepts_iteration_objects(monkeypatch, tmp_path):
    class _Iteration:
        def __init__(self, iteration, converged, hits, hmm=None):
            self.iteration = iteration
            self.converged = converged
            self.hits = hits
            self.hmm = hmm

    hit1 = _FakeHit("contigA|gene1", True, 50.0, 1e-20)
    hit2 = _FakeHit("contigA|gene2", False, 10.0, 1e-2)
    iters = [_Iteration(1, False, [hit1, hit2])]

    monkeypatch.setattr(jack.pyhmmer.easel, "SequenceFile", _FakeSequenceFile)
    monkeypatch.setattr(jack.pyhmmer.hmmer, "jackhmmer", lambda *a, **k: [iters])

    kept, summary = jack._run_jackhmmer(
        query=object(),
        alphabet=object(),
        seed_id="seedA",
        proteins_fa=tmp_path / "proteins.faa",
        iterations=2,
        inc_evalue=1e-3,
        max_evalue=1e-5,
        threads=1,
    )

    assert len(kept) == 1
    assert kept[0].contig == "contigA"
    assert kept[0].model == "seedA"
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
    monkeypatch.setattr(jack.pyhmmer.hmmer, "jackhmmer", lambda *a, **k: [[iter1, iter2]])

    kept, summary = jack._run_jackhmmer(
        query=object(),
        alphabet=object(),
        seed_id="seedB",
        proteins_fa=tmp_path / "proteins.faa",
        iterations=3,
        inc_evalue=1e-3,
        max_evalue=1e-5,
        threads=1,
    )

    assert len(kept) == 1
    assert kept[0].contig == "contigB"
    assert kept[0].model == "seedB"
    assert summary[0]["iteration"] == 1
    assert summary[1]["iteration"] == 2
    assert summary[0]["n_hits"] == 2
    assert summary[1]["n_included"] == 1


def test_run_jackhmmer_saves_last_iteration_hmm(monkeypatch, tmp_path):
    class _Iteration:
        def __init__(self, iteration, converged, hits, hmm):
            self.iteration = iteration
            self.converged = converged
            self.hits = hits
            self.hmm = hmm

    hit = _FakeHit("contigD|gene1", True, 77.0, 1e-30)
    final_hmm_path = tmp_path / "last_iteration.hmm"
    iter1 = _Iteration(1, False, [hit], _FakeHMM(b"HMM-1\n"))
    iter2 = _Iteration(2, True, [hit], _FakeHMM(b"HMM-2\n"))

    monkeypatch.setattr(jack.pyhmmer.easel, "SequenceFile", _FakeSequenceFile)
    monkeypatch.setattr(jack.pyhmmer.hmmer, "jackhmmer", lambda *a, **k: [[iter1, iter2]])

    kept, summary = jack._run_jackhmmer(
        query=object(),
        alphabet=object(),
        seed_id="seedC",
        proteins_fa=tmp_path / "proteins.faa",
        iterations=3,
        inc_evalue=1e-3,
        max_evalue=1e-5,
        threads=1,
        hmm_output_path=final_hmm_path,
    )

    assert len(kept) == 1
    assert len(summary) == 2
    assert final_hmm_path.read_bytes() == b"HMM-2\n"


def test_run_jackhmmer_skip_hmm_export_when_unavailable(monkeypatch, tmp_path):
    iter1 = [_FakeHit("contigE|gene1", True, 43.0, 1e-15)]
    final_hmm_path = tmp_path / "last_iteration.hmm"

    monkeypatch.setattr(jack.pyhmmer.easel, "SequenceFile", _FakeSequenceFile)
    monkeypatch.setattr(jack.pyhmmer.hmmer, "jackhmmer", lambda *a, **k: [[iter1]])

    kept, summary = jack._run_jackhmmer(
        query=object(),
        alphabet=object(),
        seed_id="seedD",
        proteins_fa=tmp_path / "proteins.faa",
        iterations=2,
        inc_evalue=1e-3,
        max_evalue=1e-5,
        threads=1,
        hmm_output_path=final_hmm_path,
    )

    assert len(kept) == 1
    assert len(summary) == 1
    assert not final_hmm_path.exists()


def test_choose_top_hits_any_mode():
    hits = [
        jack.Hit(contig="c1", prot_id="c1|gene1", model="seed1", bitscore=100.0, evalue=1e-20),
        jack.Hit(contig="c1", prot_id="c1|gene2", model="seed2", bitscore=90.0, evalue=1e-10),
        jack.Hit(contig="c2", prot_id="c2|gene1", model="seed2", bitscore=80.0, evalue=1e-8),
    ]
    kept, contigs = jack._choose_top_hits_per_contig(
        hits,
        top_per_contig=1,
        combine_mode="any",
        min_seed_hits=1,
        total_seeds=2,
    )

    assert set(contigs) == {"c1", "c2"}
    assert len(kept) == 2


def test_choose_top_hits_all_mode():
    hits = [
        jack.Hit(contig="c1", prot_id="c1|gene1", model="seed1", bitscore=100.0, evalue=1e-20),
        jack.Hit(contig="c1", prot_id="c1|gene2", model="seed2", bitscore=90.0, evalue=1e-10),
        jack.Hit(contig="c2", prot_id="c2|gene1", model="seed1", bitscore=80.0, evalue=1e-8),
    ]
    kept, contigs = jack._choose_top_hits_per_contig(
        hits,
        top_per_contig=5,
        combine_mode="all",
        min_seed_hits=1,
        total_seeds=2,
    )

    assert contigs == ["c1"]
    assert len(kept) == 2
    assert set(h.model for h in kept) == {"seed1", "seed2"}


def test_choose_top_hits_threshold_mode():
    hits = [
        jack.Hit(contig="c1", prot_id="c1|gene1", model="seed1", bitscore=100.0, evalue=1e-20),
        jack.Hit(contig="c1", prot_id="c1|gene2", model="seed2", bitscore=90.0, evalue=1e-10),
        jack.Hit(contig="c1", prot_id="c1|gene3", model="seed3", bitscore=85.0, evalue=1e-9),
        jack.Hit(contig="c2", prot_id="c2|gene1", model="seed1", bitscore=80.0, evalue=1e-8),
    ]
    kept, contigs = jack._choose_top_hits_per_contig(
        hits,
        top_per_contig=2,
        combine_mode="threshold",
        min_seed_hits=2,
        total_seeds=3,
    )

    assert contigs == ["c1"]
    assert len(kept) == 2


def test_read_seed_queries_rejects_duplicate_ids(monkeypatch, tmp_path):
    monkeypatch.setattr(
        jack,
        "_read_fasta",
        lambda _p: [("seed1", "AAAA"), ("seed1", "BBBB")],
    )

    seed_fa = tmp_path / "seeds.faa"
    seed_fa.write_text("")

    try:
        jack._read_seed_queries(seed_fa)
        assert False, "Expected ValueError for duplicate seed IDs"
    except ValueError as exc:
        assert "Duplicate seed sequence ID" in str(exc)
