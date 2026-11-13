from __future__ import annotations

from pathlib import Path

import pyhmmer

from phu.commands.screen import Hit, _hmmsearch


def _write_fixture_hmm(path: Path) -> None:
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences = [
        pyhmmer.easel.TextSequence(name=f"s{i}", sequence=sequence)
        for i, sequence in enumerate(
            (
                "MKKLLAAVAGAAAAPAAA",
                "MKKLLAAVAGAAAAPAAA",
                "GKKLLAAVAGAAAAPAAA",
            )
        )
    ]
    msa = pyhmmer.easel.TextMSA(name="tiny_model", sequences=sequences).digitize(
        alphabet
    )
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build_msa(
        msa, pyhmmer.plan7.Background(alphabet)
    )
    with path.open("wb") as handle:
        hmm.write(handle)


def _hit_signature(hits: list[Hit]) -> list[tuple[object, ...]]:
    return [
        (
            hit.contig,
            hit.prot_id,
            hit.model,
            hit.bitscore,
            hit.evalue,
            hit.domain_bitscore,
            hit.hmm_coverage,
            hit.domain_i_evalue,
            hit.hmm_from,
            hit.hmm_to,
            hit.target_from,
            hit.target_to,
        )
        for hit in hits
    ]


def test_hmmsearch_results_are_reproducible_across_cpu_counts(tmp_path: Path) -> None:
    """The production HMM search returns the same scientific hits across CPUs."""
    hmm_path = tmp_path / "tiny.hmm"
    proteins_path = tmp_path / "proteins.faa"
    _write_fixture_hmm(hmm_path)
    proteins_path.write_text(
        ">p1\nMAAAAAAAAAAAAAAAAAAA\n"
        ">p2\nGAAAAAAAAAAAAAAAAAAA\n"
        ">p3\nVVVVVVVVVVVVVVVVVVV\n"
    )

    one_cpu = list(
        _hmmsearch(
            [hmm_path],
            proteins_path,
            {},
            threads=1,
            keep_domtbl=False,
            cut_ga=False,
        )
    )
    two_cpu = list(
        _hmmsearch(
            [hmm_path],
            proteins_path,
            {},
            threads=2,
            keep_domtbl=False,
            cut_ga=False,
        )
    )

    assert _hit_signature(one_cpu) == _hit_signature(two_cpu)


def _build_hmm(name: str, sequence: str) -> pyhmmer.plan7.HMM:
    alphabet = pyhmmer.easel.Alphabet.amino()
    msa = pyhmmer.easel.TextMSA(
        name=name,
        sequences=[
            pyhmmer.easel.TextSequence(name=f"s{i}", sequence=sequence)
            for i in range(3)
        ],
    ).digitize(alphabet)
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build_msa(
        msa, pyhmmer.plan7.Background(alphabet)
    )
    return hmm


def _domtbl_rows(path: Path) -> list[list[str]]:
    return [
        line.split()
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ]


def test_mixed_mode_domtbl_keeps_hits_of_every_model_in_file(tmp_path: Path) -> None:
    """A multi-model HMM file yields one domain table with all models' hits."""
    seq_a = "MKKLLAAVAGWAAAPAAACDEFGHIKLMNPQRSTVWY"
    seq_b = "WYVTSRQPNMLKIHGFEDCWWHHYYFFMMKKRRDDEE"
    hmm_path = tmp_path / "multi.hmm"
    with hmm_path.open("wb") as handle:
        _build_hmm("model_a", seq_a).write(handle)
        _build_hmm("model_b", seq_b).write(handle)
    proteins_path = tmp_path / "proteins.faa"
    proteins_path.write_text(f">pa\n{seq_a}\n>pb\n{seq_b}\n")
    domtbl_path = tmp_path / "multi-hmmsearch-domains.tbl"

    hits = list(
        _hmmsearch(
            [hmm_path],
            proteins_path,
            {"multi": domtbl_path},
            threads=1,
            hmm_mode="mixed",
            keep_domtbl=True,
            cut_ga=False,
        )
    )

    assert {hit.model for hit in hits} == {"model_a", "model_b"}
    text = domtbl_path.read_text()
    # Column 4 of a domain table is the query (model) name.
    assert {row[3] for row in _domtbl_rows(domtbl_path)} == {"model_a", "model_b"}
    # A single header block, at the top of the file.
    assert text.count("# target name") == 1
    assert text.startswith("#")


def test_pure_mode_single_model_domtbl_has_header_and_rows(tmp_path: Path) -> None:
    """Single-model files still get a plain headed domain table."""
    seq_a = "MKKLLAAVAGWAAAPAAACDEFGHIKLMNPQRSTVWY"
    hmm_path = tmp_path / "single.hmm"
    with hmm_path.open("wb") as handle:
        _build_hmm("model_a", seq_a).write(handle)
    proteins_path = tmp_path / "proteins.faa"
    proteins_path.write_text(f">pa\n{seq_a}\n")
    domtbl_path = tmp_path / "single-hmmsearch-domains.tbl"

    list(
        _hmmsearch(
            [hmm_path],
            proteins_path,
            {"single": domtbl_path},
            threads=1,
            keep_domtbl=True,
            cut_ga=False,
        )
    )

    text = domtbl_path.read_text()
    assert text.count("# target name") == 1
    assert [row[3] for row in _domtbl_rows(domtbl_path)] == ["model_a"]
