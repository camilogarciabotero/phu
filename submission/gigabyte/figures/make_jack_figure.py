"""Draw the phu jack workflow figure as SVG.

Usage:
    python make_jack_figure.py > jack-workflow.svg
    rsvg-convert -f pdf -o jack-workflow.pdf jack-workflow.svg

The contigs, proteins and hits in the figure are illustrative, not results.
"""

from figure_kit import (
    BLUE, GREEN, INK, LINE, MUTED, ORANGE, arrow, box, chip, circle, contig,
    decision_table, finish, gene, outputs, panel, path, pill, rect, start, text,
)

SEEDS = [("S1", BLUE), ("S2", GREEN)]

start("phu jack: selecting contigs by iterative search from seed proteins")

# ---- A ------------------------------------------------------------------
ax, ay = 30, 64
panel(ax, ay, 655, 330, "A", "Seed proteins and gene prediction")
text(ax + 28, ay + 66, "Seed markers (protein FASTA, one or more sequences)", 12, "bold")
for i, (name, color) in enumerate(SEEDS):
    chip(ax + 28 + i * 60, ay + 78, 52, name, color, 24, 12)
text(ax + 160, ay + 95, "used when a trusted marker protein exists but no profile HMM does", 11.5, fill=MUTED)
text(ax + 28, ay + 150, "Contigs (FASTA)", 12, "bold")
contig(ax + 28, ay + 171, 250)
for i, w in enumerate([70, 58, 76]):
    gene(ax + 40 + sum([70, 58, 76][:i]) + i * 8, ay + 160, w)
arrow(ax + 290, ay + 171, ax + 338, ay + 171)
box(ax + 342, ay + 148, 286, 46, "pyrodigal gene prediction", "meta mode by default; --mode, --ttable")
arrow(ax + 485, ay + 196, ax + 485, ay + 222)
for i in range(5):
    pill(ax + 348 + i * 56, ay + 226, 50, f"p{i + 1}")
text(ax + 28, ay + 244, "Predicted proteins (≥ 30 aa)", 12)
text(ax + 328, ay + 296, "The prediction cache is shared with phu screen and phu avger", 11.5, "normal", "middle", MUTED, "italic")

# ---- B ------------------------------------------------------------------
bx, by = 715, 64
panel(bx, by, 655, 330, "B", "Iterative search (PyHMMER jackhmmer), one seed at a time")
DOTS = [3, 6, 9]
for k, n in enumerate(DOTS):
    x = bx + 40 + k * 205
    rect(x, by + 62, 170, 150, "#ffffff", LINE)
    text(x + 85, by + 82, f"iteration {k + 1}", 12, "bold", "middle")
    chip(x + 55, by + 92, 60, "S1" if k == 0 else "profile", BLUE, 18, 10)
    for d in range(9):
        cxp = x + 32 + (d % 3) * 52
        cyp = by + 134 + (d // 3) * 26
        new = DOTS[k - 1] <= d < n if k else d < n
        if d < n:
            circle(cxp, cyp, 9, BLUE if not new else ORANGE)
        else:
            circle(cxp, cyp, 9, "#ffffff", LINE)
    if k < 2:
        arrow(x + 172, by + 137, x + 203, by + 137)
circle(bx + 48, by + 236, 7, ORANGE)
text(bx + 62, by + 240, "newly included hit (E-value ≤ 1e-3, --inc-evalue)", 11.5, fill=MUTED)
circle(bx + 388, by + 236, 7, BLUE)
text(bx + 402, by + 240, "already in the profile", 11.5, fill=MUTED)
text(bx + 28, by + 270, "After each round the included hits are aligned into a profile, which is searched again.", 11.5, fill=MUTED)
text(bx + 28, by + 288, "The search runs for up to 5 rounds (--iterations) and so reaches more distant relatives", 11.5, fill=MUTED)
text(bx + 28, by + 306, "than the seed alone. Each seed is searched independently.", 11.5, fill=MUTED)

# ---- C ------------------------------------------------------------------
cx, cy = 30, 414
panel(cx, cy, 655, 440, "C", "Final hits and contig selection")
box(cx + 28, cy + 54, 290, 46, "Final hits: E-value ≤ 1e-5", "--max-evalue")
box(cx + 338, cy + 54, 290, 46, "Best hit per contig by bit score", "--top-per-contig 1")
text(cx + 28, cy + 124, "With several seeds, --combine-mode decides which contigs are kept", 11.5, fill=MUTED)
rows = [("c1", {0, 1}), ("c2", {0}), ("c3", set()), ("c4", {1})]
rules = [
    ("any (default)", lambda p: len(p) >= 1),
    ("all", lambda p: len(p) == 2),
    ("threshold ≥ 2", lambda p: len(p) >= 2),
]
end = decision_table(cx + 60, cy + 138, [s for s, _ in SEEDS], [c for _, c in SEEDS], rows, rules, col_w=70)
for i, (name, meaning) in enumerate([
    ("any", "at least one seed's search hits the contig"),
    ("all", "every seed's search hits the contig"),
    ("threshold", "at least --min-seed-hits seeds hit the contig"),
]):
    y = end + 24 + i * 22
    text(cx + 28, y, name, 12, "bold")
    text(cx + 110, y, meaning, 12, fill=MUTED)
arrow(cx + 328, end + 80, cx + 328, end + 98, ORANGE)
chip(cx + 228, end + 102, 200, "matched contigs", ORANGE, 28, 13)

# ---- D ------------------------------------------------------------------
dx, dy = 715, 414
panel(dx, dy, 655, 440, "D", "Optional: keep the final profile")
text(dx + 28, dy + 66, "The profile of the last round can be saved and used as a query elsewhere", 11.5, fill=MUTED)
chip(dx + 60, dy + 100, 60, "S1", BLUE, 26, 12)
arrow(dx + 124, dy + 113, dx + 176, dy + 113)
box(dx + 180, dy + 90, 200, 46, "jackhmmer rounds", "panel B")
arrow(dx + 384, dy + 113, dx + 436, dy + 113)
chip(dx + 440, dy + 99, 178, "final-iteration.hmm", INK, 28, 12)
text(dx + 529, dy + 146, "--save-hmm", 11, "normal", "middle", MUTED)
path(f"M{dx + 529},{dy + 156} V{dy + 212}", MUTED, 1.6, True, "5,4")
rect(dx + 180, dy + 216, 438, 92, "#ffffff", LINE)
text(dx + 399, dy + 240, "phu screen -i other_contigs.fna final-iteration.hmm", 12, "bold", "middle")
text(dx + 399, dy + 262, "the saved profile screens other assemblies directly,", 11.5, "normal", "middle", MUTED)
text(dx + 399, dy + 280, "without repeating the iterative search", 11.5, "normal", "middle", MUTED)
text(dx + 28, dy + 346, "With one seed the file is final-iteration.hmm; with several, one profile per seed", 11.5, fill=MUTED)
text(dx + 28, dy + 364, "is written under final-iteration-hmms/.", 11.5, fill=MUTED)
text(dx + 28, dy + 398, "jackhmmer-iterations.tsv records how many hits each round added, which shows", 11.5, fill=MUTED)
text(dx + 28, dy + 416, "whether the search converged or kept growing.", 11.5, fill=MUTED)

outputs([
    ("matched-contigs.fna", "contigs that passed"),
    ("matched-contig-ids.txt", "their names"),
    ("jackhmmer-hits.tsv", "final hits and filters applied"),
    ("jackhmmer-iterations.tsv", "summary of each round"),
    ("final-iteration.hmm", "optional"),
])
finish(
    "Hits gained over the rounds are sequence-similarity evidence, not a validated functional or taxonomic assignment.",
    "Contigs, proteins and hits shown are illustrative",
)
