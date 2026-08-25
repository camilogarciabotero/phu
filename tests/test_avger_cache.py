from pathlib import Path

from phu.avger_cache import AvgerInputs, compute_avger_cache_key, get_or_run_avger


def test_compute_avger_cache_key_changes_for_relevant_inputs(tmp_path: Path):
    fasta = tmp_path / "a.fa"
    fasta.write_text(">ctg1\nATGGGAAATTT\n")

    a = AvgerInputs(
        input_contigs=fasta,
        mode="meta",
        threads=2,
        max_evalue=1e-5,
        min_gene_len=90,
        min_protein_len_aa=30,
        translation_table=11,
        keep_all_hits=False,
        use_vscore=True,
        require_flank_support=False,
    )
    b = AvgerInputs(
        input_contigs=fasta,
        mode="meta",
        threads=2,
        max_evalue=1e-4,
        min_gene_len=90,
        min_protein_len_aa=30,
        translation_table=11,
        keep_all_hits=False,
        use_vscore=True,
        require_flank_support=False,
    )

    assert compute_avger_cache_key(a) != compute_avger_cache_key(b)


def test_get_or_run_avger_reuses_cached_result(tmp_path: Path):
    contigs = tmp_path / "contigs.fa"
    contigs.write_text(">ctg1\nATGGGAAATTT\n")

    calls: list[Path] = []

    def fake_run(inputs: AvgerInputs, output_folder: Path) -> Path:
        calls.append(output_folder)
        out = output_folder / "best_hits.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("protein_id\tcontig_id\nctg1|gene1\tctg1\n")
        return out

    first_out, first_hit = get_or_run_avger(
        AvgerInputs(
            input_contigs=contigs,
            mode="meta",
            threads=1,
            max_evalue=1e-5,
            min_gene_len=90,
            min_protein_len_aa=30,
            translation_table=11,
            keep_all_hits=False,
            use_vscore=True,
            require_flank_support=False,
        ),
        tmp_path / "run1",
        fake_run,
    )

    second_out, second_hit = get_or_run_avger(
        AvgerInputs(
            input_contigs=contigs,
            mode="meta",
            threads=1,
            max_evalue=1e-5,
            min_gene_len=90,
            min_protein_len_aa=30,
            translation_table=11,
            keep_all_hits=False,
            use_vscore=True,
            require_flank_support=False,
        ),
        tmp_path / "run2",
        fake_run,
    )

    assert first_hit is False
    assert second_hit is True
    assert first_out.read_text() == "protein_id\tcontig_id\nctg1|gene1\tctg1\n"
    assert second_out.read_text() == "protein_id\tcontig_id\nctg1|gene1\tctg1\n"
    assert len(calls) == 1
