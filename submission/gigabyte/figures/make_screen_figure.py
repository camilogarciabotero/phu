"""Draw the phu screen workflow figure as SVG.

Usage:
    python make_screen_figure.py > screen-workflow.svg
    rsvg-convert -f pdf -o screen-workflow.pdf screen-workflow.svg

The contigs, genes and hits in the figure are illustrative, not results.
"""

from figure_kit import (
    BLUE, GREEN, INK, LINE, MUTED, ORANGE, PURPLE, arrow, box, chip, contig,
    decision_table, finish, gene, outputs, panel, pill, rect, start, tag, text,
)

MODELS = [("M1", BLUE), ("M2", GREEN), ("M3", PURPLE)]
# contig name, gene widths, {gene index: model index}
CONTIGS = [
    ("c1", [70, 58, 76, 64, 72], {0: 0, 2: 1, 4: 2}),
    ("c2", [64, 80, 60, 70], {1: 0}),
    ("c3", [74, 62, 68], {}),
    ("c4", [60, 72, 66, 78], {1: 1, 3: 2}),
]

start("phu screen: selecting contigs that carry chosen protein families")

# ---- A ------------------------------------------------------------------
ax, ay = 30, 64
panel(ax, ay, 655, 330, "A", "Query models and gene prediction")
text(ax + 28, ay + 66, "Query: any mix of", 12, "bold")
for i, (label, detail) in enumerate([
    ("HMM files", "capsid.hmm"),
    ("Pfam accessions", "PF00589"),
    ("KOfam orthologs", "K07497"),
    ("dbCAN, PULs", "GH13, PUL0001"),
]):
    x = ax + 28 + i * 152
    rect(x, ay + 76, 144, 44, "#ffffff", LINE)
    text(x + 72, ay + 94, label, 11.5, "bold", "middle")
    text(x + 72, ay + 111, detail, 10.5, "normal", "middle", MUTED)
text(ax + 28, ay + 142, "Database names are resolved from the pinned phu-db archive into models:", 11.5, fill=MUTED)
for i, (name, color) in enumerate(MODELS):
    chip(ax + 500 + i * 46, ay + 128, 40, name, color)
contig(ax + 28, ay + 189, 250)
for i, w in enumerate([70, 58, 76]):
    gene(ax + 40 + sum([70, 58, 76][:i]) + i * 8, ay + 178, w)
text(ax + 28, ay + 168, "Contigs (FASTA)", 12, "bold")
arrow(ax + 290, ay + 189, ax + 338, ay + 189)
box(ax + 342, ay + 166, 286, 46, "pyrodigal / pyrodigal-gv", "meta mode by default; --mode, --ttable")
arrow(ax + 485, ay + 214, ax + 485, ay + 240)
for i in range(5):
    pill(ax + 348 + i * 56, ay + 244, 50, f"p{i + 1}")
text(ax + 28, ay + 262, "Predicted proteins (≥ 30 aa)", 12)
text(ax + 328, ay + 304, "Predictions are cached: a new screen of the same contigs with other models skips this step", 11.5, "normal", "middle", MUTED, "italic")

# ---- B ------------------------------------------------------------------
bx, by = 715, 64
panel(bx, by, 655, 330, "B", "Profile search (PyHMMER hmmsearch)")
for i, (name, color) in enumerate(MODELS):
    chip(bx + 490 + i * 46, by + 16, 40, name, color)
for r, (name, widths, hits) in enumerate(CONTIGS):
    y = by + 72 + r * 46
    text(bx + 28, y + 15, name, 12, "bold")
    contig(bx + 56, y + 11, 30 + sum(widths) + 8 * len(widths))
    x = bx + 66
    for g, w in enumerate(widths):
        if g in hits:
            label, color = MODELS[hits[g]]
            tag(x + 4, y - 15, label, color)
            gene(x, y, w, fill=color, stroke=color, color="#ffffff", label=f"g{g + 1}")
        else:
            gene(x, y, w, label=f"g{g + 1}")
        x += w + 8
text(bx + 28, by + 270, "A hit counts when it passes the model's gathering cutoff (--cut-ga, KOfam thresholds)", 11.5, fill=MUTED)
text(bx + 28, by + 288, "or, for models without one, E-value ≤ 1e-5 (--max-evalue) and --min-bitscore if given.", 11.5, fill=MUTED)
text(bx + 28, by + 306, "The best hit per model is kept on each contig (--top-per-contig 1).", 11.5, fill=MUTED)

# ---- C ------------------------------------------------------------------
cx, cy = 30, 414
panel(cx, cy, 655, 440, "C", "Combine hits per contig (--combine-mode)")
text(cx + 28, cy + 66, "Which models hit each contig decides whether the contig is kept", 11.5, fill=MUTED)
rows = [(name, set(hits.values())) for name, _, hits in CONTIGS]
rules = [
    ("any (default)", lambda p: len(p) >= 1),
    ("all", lambda p: len(p) == 3),
    ("threshold ≥ 2", lambda p: len(p) >= 2),
]
end = decision_table(cx + 28, cy + 88, [m for m, _ in MODELS], [c for _, c in MODELS], rows, rules)
for i, (name, meaning) in enumerate([
    ("any", "at least one model hits the contig"),
    ("all", "every queried model hits the contig"),
    ("threshold", "at least --min-hmm-hits models hit the contig"),
]):
    y = end + 34 + i * 26
    text(cx + 28, y, name, 12, "bold")
    text(cx + 110, y, meaning, 12, fill=MUTED)
arrow(cx + 328, end + 108, cx + 328, end + 132, ORANGE)
chip(cx + 228, end + 136, 200, "matched contigs", ORANGE, 28, 13)
text(cx + 328, end + 184, "The query and the rule used are written to .phu/query.json", 11.5, "normal", "middle", MUTED, "italic")

# ---- D ------------------------------------------------------------------
dx, dy = 715, 414
panel(dx, dy, 655, 440, "D", "Optional: turn the hits into new models")
text(dx + 28, dy + 66, "Proteins that hit each model, from kept contigs only, can seed a model of your own data", 11.5, fill=MUTED)
box(dx + 28, dy + 84, 280, 46, "target-proteins/{model}.faa", "--save-target-proteins")
arrow(dx + 168, dy + 132, dx + 168, dy + 156)
box(dx + 28, dy + 160, 280, 46, "Align with FAMSA", "in-process (pyfamsa); no aligner to install")
arrow(dx + 168, dy + 208, dx + 168, dy + 232)
rect(dx + 28, dy + 236, 280, 46, "#ffffff", LINE, 6, 1.2, "5,3")
text(dx + 168, dy + 256, "Trim with trimAl (optional)", 12.5, "bold", "middle")
text(dx + 168, dy + 273, "--trim-target-msa gappyout | strict | strictplus", 10.5, "normal", "middle", MUTED)
arrow(dx + 168, dy + 284, dx + 168, dy + 308)
box(dx + 28, dy + 312, 280, 46, "Build the HMM (PyHMMER)", "--save-target-hmms")
arrow(dx + 168, dy + 360, dx + 168, dy + 384)
chip(dx + 58, dy + 388, 220, "target-hmms/{model}.hmm", INK, 28, 12)

# alignment sketch
sx, sy = dx + 372, dy + 150
text(sx, sy - 14, "alignment", 11, "bold", fill=MUTED)
ROWS = ["1111101111110011", "1111111101110001", "0011101111110111", "1111100111111011", "1101111111100011"]
for r, row in enumerate(ROWS):
    for c, v in enumerate(row):
        rect(sx + c * 14, sy + r * 14, 12, 11, BLUE if v == "1" else "#ffffff", LINE if v == "0" else BLUE, 2, 0.6)
rect(sx + 14 * 13 - 1, sy - 7, 14 * 2 + 1, 1.6, ORANGE, "none", 0, 0)
text(sx + 14 * 14, sy - 12, "gappy", 10, "normal", "middle", ORANGE)
text(sx, sy + 96, "trimmed", 11, "bold", fill=MUTED)
KEEP = [c for c in range(16) if c not in (13, 14)]
for r, row in enumerate(ROWS):
    for i, c in enumerate(KEEP):
        v = row[c]
        rect(sx + i * 14, sy + 108 + r * 14, 12, 11, BLUE if v == "1" else "#ffffff", LINE if v == "0" else BLUE, 2, 0.6)
text(sx, sy + 204, "Trimming shortens the model and", 11, fill=MUTED)
text(sx, sy + 220, "lowers bit scores; off by default.", 11, fill=MUTED)
text(dx + 300, dy + 400, "reusable as a query in", 11, fill=MUTED, style="italic")
text(dx + 300, dy + 416, "a later phu screen", 11, fill=MUTED, style="italic")

outputs([
    ("matched-contigs.fna", "contigs that passed"),
    ("matched-contig-ids.txt", "their names"),
    ("{hmm}-hmmsearch-domains.tbl", "per-domain hits, one per HMM file"),
    (".phu/query.json", "query, rule, database versions"),
    ("target-proteins/, target-hmms/", "optional"),
])
finish(
    "A matched contig carries sequence similarity to the queried families; it is not a functional or taxonomic assignment.",
    "Contigs, genes and hits shown are illustrative",
)
