"""The documented quickstart gives the documented result."""

from pathlib import Path

from phu.commands import screen

QUICKSTART = Path(__file__).parent.parent / "docs" / "quickstart"


def test_quickstart_keeps_phix174_and_drops_the_decoy(tmp_path, monkeypatch):
    monkeypatch.setenv("PHU_CACHE", "off")
    out = tmp_path / "qs-out"
    screen._screen(
        screen.ScreenConfig(
            input_contigs=QUICKSTART / "contigs.fna",
            hmms=[QUICKSTART / "major_head.hmm"],
            outdir=out,
        )
    )

    assert (out / "matched-contig-ids.txt").read_text().split() == ["NC_001422.1"]
    assert sorted(p.name for p in out.iterdir()) == [
        ".phu",
        "major_head-hmmsearch-domains.tbl",
        "matched-contig-ids.txt",
        "matched-contigs.fna",
    ]
