"""Draw the phu cluster and phu normalize-lineage figure as SVG.

Usage:
    python make_cluster_lineage_figure.py > cluster-lineage-workflow.svg
    rsvg-convert -f pdf -o cluster-lineage-workflow.pdf cluster-lineage-workflow.svg

The sequences, clusters and lineage values in the figure are illustrative.
"""

from figure_kit import (
    BLUE, GREEN, INK, LINE, MUTED, ORANGE, PURPLE, arrow, box, chip, circle,
    contig, finish, panel, pill, rect, start, text,
)

MONO = "Liberation Mono, DejaVu Sans Mono, Consolas, monospace"


def code(x, y, s, size=11.5, fill=INK, anchor="start", weight="normal"):
    from figure_kit import out
    from xml.sax.saxutils import escape
    out.append(
        f'<text x="{x}" y="{y}" font-size="{size}" font-family="{MONO}" '
        f'font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">{escape(s)}</text>'
    )


start("phu cluster and phu normalize-lineage: grouping genomes and tidying their names")

# ---- A: cluster -----------------------------------------------------------
ax, ay = 30, 64
panel(ax, ay, 1340, 414, "A", "phu cluster: three vclust recipes")

text(ax + 28, ay + 66, "Viral genomes (FASTA)", 12, "bold")
contig(ax + 28, ay + 88, 190)
text(ax + 28, ay + 112, "all-vs-all, in one command", 11, fill=MUTED)

arrow(ax + 232, ay + 88, ax + 282, ay + 88)
box(ax + 286, ay + 50, 206, 46, "1 · prefilter", "k-mer pairs worth aligning")
arrow(ax + 492, ay + 73, ax + 524, ay + 73)
box(ax + 528, ay + 50, 206, 46, "2 · align", "ANI / TANI per pair (Lempel-Ziv)")
arrow(ax + 734, ay + 73, ax + 766, ay + 73)
box(ax + 770, ay + 50, 206, 46, "3 · cluster", "one algorithm per recipe")
text(ax + 1000, ay + 78, "vclust [10]", 12, "bold", "start", PURPLE)

# Recipe cards.
cards = (
    (
        "dereplication", BLUE,
        ("cd-hit algorithm", "ANI ≥ 95% · coverage ≥ 85%", "prefilter identity 95%"),
        "one representative per group of near-duplicates",
        "dereplicated-representatives.fna",
    ),
    (
        "votu", GREEN,
        ("Leiden algorithm", "ANI ≥ 95% · coverage ≥ 85%", "prefilter identity 95%"),
        "viral OTUs with MIUViG-style defaults",
        "clusters.tsv · cluster-representatives.fna",
    ),
    (
        "species", ORANGE,
        ("complete linkage", "total ANI (TANI) ≥ 95%", "prefilter identity 70%"),
        "species groups with ICTV-style defaults",
        "species-clusters.tsv · cluster-representatives.fna",
    ),
)
for i, (name, color, lines, meaning, files) in enumerate(cards):
    x = ax + 48 + i * 432
    y = ay + 148
    rect(x, y, 400, 214, "#ffffff", LINE)
    chip(x + 16, y + 16, 150, f"--mode {name}", color, 24, 12)
    for j, line in enumerate(lines):
        text(x + 16, y + 66 + j * 19, line, 11.5)
    # A small illustrative clustering: nine dots, grouped differently per card.
    groups = ([(0, 1), (2, 3, 4), (5,), (6, 7, 8)], [(0, 1, 2, 3), (4, 5), (6, 7, 8)],
              [(0, 1, 2, 3, 4, 5), (6, 7, 8)])[i]
    palette = (BLUE, GREEN, ORANGE, PURPLE)
    for g, members in enumerate(groups):
        for d in members:
            cxp = x + 240 + (d % 3) * 34
            cyp = y + 64 + (d // 3) * 30
            circle(cxp, cyp, 8, palette[g % 4] if d == members[0] else "#ffffff",
                   palette[g % 4], 1.8)
    text(x + 273, y + 158, "filled dot = representative", 10, "normal", "middle", MUTED)
    text(x + 16, y + 176, meaning, 11, fill=MUTED)
    code(x + 16, y + 198, files, 10.5, MUTED)

text(ax + 670, ay + 394,
     "anything else goes through --vclust-params, straight to vclust",
     11.5, "normal", "middle", MUTED, "italic")

# ---- B: normalize-lineage ---------------------------------------------------
bx, by = 30, 502
panel(bx, by, 1340, 350, "B", "phu normalize-lineage: compact codes for vConTACT predictions")

text(bx + 28, by + 64, "A vConTACT3 assignment table keeps its columns; only the lineage cells change:",
     12)

rows = (
    ("novel_genus_1_of_Viruses", "Viruses:NG1"),
    ("novel_genus_1_of_novel_family_2_of_unplaced_order_within_Caudoviricetes",
     "Caudoviricetes:UO:NF2:NG1"),
    ("unplaced_phylum_within_unplaced_kingdom_within_default", "DEF:UK:UP"),
)
for i, (before, after) in enumerate(rows):
    y = by + 102 + i * 42
    rect(bx + 28, y - 20, 650, 30, "#ffffff", LINE, 4)
    code(bx + 40, y, before, 11.5)
    arrow(bx + 692, y - 5, bx + 740, y - 5)
    rect(bx + 748, y - 20, 300, 30, "#ffffff", PURPLE, 4)
    code(bx + 760, y, after, 12, PURPLE, weight="bold")

legend_y = by + 238
text(bx + 28, legend_y, "How a value is rewritten", 12, "bold")
for j, line in enumerate((
    "the anchor (the last known taxon) comes first; codes follow from the highest rank to the lowest",
    "N = novel node with its number · U = unplaced node (no number) · DEF = the placeholder anchor 'default'",
    "rank letters: R K P C O F SF G (realm to genus); each node must be of a lower rank than the one before it",
)):
    text(bx + 46, legend_y + 22 + j * 19, line, 11.5, fill=MUTED)

nx = bx + 748
text(nx, legend_y, "Safety", 12, "bold")
for j, line in enumerate((
    "unparsed values are kept and counted in the QA summary; --strict makes them fatal",
    "already-compact values are left alone, so the command is idempotent",
    "--add-lineage appends one combined lineage column; nothing else is added or renamed",
)):
    text(nx + 18, legend_y + 22 + j * 19, line, 11.5, fill=MUTED)

finish(
    "Clusters and lineage values shown are illustrative.",
    "Recipes follow the vclust use cases; defaults follow MIUViG and ICTV practice "
    "· lineage grammar: vConTACT3 novel_/unplaced_ chains",
)
