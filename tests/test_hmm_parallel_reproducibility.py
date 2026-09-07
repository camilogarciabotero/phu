from __future__ import annotations

from pathlib import Path

import pyhmmer

from phu.screen import Hit, _hmmsearch


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
        ">p1\nMKKLLAAVAGAAAAPAAA\n"
        ">p2\nGKKLLAAVAGAAAAPAAA\n"
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
