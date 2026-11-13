"""Draw the phu avger workflow figure as SVG.

Usage:
    python make_avger_figure.py > avger-workflow.svg
    rsvg-convert -f pdf -o avger-workflow.pdf avger-workflow.svg

The gene values in the figure are illustrative, not results.
"""

from xml.sax.saxutils import escape

INK, MUTED, LINE, PANEL = "#1f2933", "#52606d", "#9aa5b1", "#f7f9fb"
PFAM, KOFAM, PHROG = "#2f6fb0", "#2e8b57", "#8e44ad"
CAND, WARN, GENE = "#d35400", "#c0392b", "#cbd5e0"
FONT = "Liberation Sans, Arial, Helvetica, DejaVu Sans, sans-serif"
out = []


def rect(x, y, w, h, fill="none", stroke=LINE, rx=6, sw=1.2, dash=None):
    extra = f' stroke-dasharray="{dash}"' if dash else ""
    out.append(
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" '
        f'stroke="{stroke}" stroke-width="{sw}"{extra}/>'
    )


def text(x, y, s, size=12, weight="normal", anchor="start", fill=INK, style="normal"):
    out.append(
        f'<text x="{x}" y="{y}" font-size="{size}" font-weight="{weight}" '
        f'font-style="{style}" text-anchor="{anchor}" fill="{fill}">{escape(s)}</text>'
    )


def arrow(x1, y1, x2, y2, color=MUTED, sw=1.6):
    marker = {MUTED: "head", CAND: "head-cand", WARN: "head-warn"}[color]
    out.append(
        f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" '
        f'stroke-width="{sw}" marker-end="url(#{marker})"/>'
    )


def gene(x, y, w, h=26, fill=GENE, stroke=MUTED, label="", sw=1.2, color=INK):
    tip = 10
    pts = (
        f"{x},{y} {x + w - tip},{y} {x + w},{y + h / 2} "
        f"{x + w - tip},{y + h} {x},{y + h}"
    )
    out.append(
        f'<polygon points="{pts}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>'
    )
    if label:
        text(x + (w - tip / 2) / 2, y + h / 2 + 4, label, 11, "bold", "middle", color)


def contig(x, y, width):
    out.append(
        f'<line x1="{x}" y1="{y}" x2="{x + width}" y2="{y}" stroke="{LINE}" '
        f'stroke-width="3"/>'
    )


def panel(x, y, w, h, letter, title):
    rect(x, y, w, h, PANEL, LINE, 10)
    out.append(f'<circle cx="{x + 24}" cy="{y + 26}" r="13" fill="{INK}"/>')
    text(x + 24, y + 31, letter, 14, "bold", "middle", "#ffffff")
    text(x + 46, y + 31, title, 15, "bold")


def chip(x, y, w, label, color, h=20, size=11):
    rect(x, y, w, h, color, color, 10, 0)
    text(x + w / 2, y + h / 2 + 4, label, size, "bold", "middle", "#ffffff")


WIDTHS = [78, 70, 86, 74, 82, 68, 80]
GAP = 8


def gene_xs(x0):
    xs, x = [], x0
    for w in WIDTHS:
        xs.append(x)
        x += w + GAP
    return xs


def tags(x, y, dbs):
    for i, db in enumerate(dbs):
        color = {"Pf": PFAM, "KO": KOFAM, "PH": PHROG}[db]
        rect(x + i * 24, y, 22, 13, color, color, 3, 0)
        text(x + i * 24 + 11, y + 10, db, 8.5, "bold", "middle", "#ffffff")


W, H = 1400, 1000
out.append(
    f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" '
    f'height="{H}" font-family="{FONT}">'
)
out.append("<defs>")
for name, color in (("head", MUTED), ("head-cand", CAND), ("head-warn", WARN)):
    out.append(
        f'<marker id="{name}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
        f'markerHeight="7" orient="auto-start-reverse">'
        f'<path d="M0,0 L10,5 L0,10 z" fill="{color}"/></marker>'
    )
out.append("</defs>")
rect(0, 0, W, H, "#ffffff", "none", 0, 0)
text(
    W / 2,
    40,
    "phu avger: finding auxiliary viral gene (AVG) candidates in viral contigs",
    21,
    "bold",
    "middle",
)

# ---- A: gene prediction -------------------------------------------------
ax, ay = 30, 64
panel(ax, ay, 655, 330, "A", "Gene prediction")
text(ax + 28, ay + 66, "Trusted viral contig (FASTA)", 12, fill=MUTED)
contig(ax + 28, ay + 99, 600)
for i, x in enumerate(gene_xs(ax + 40)):
    gene(x, ay + 86, WIDTHS[i], label=f"g{i + 1}")
arrow(ax + 328, ay + 124, ax + 328, ay + 156)
rect(ax + 118, ay + 160, 420, 44, "#ffffff", LINE)
text(ax + 328, ay + 179, "pyrodigal-gv, meta mode", 12.5, "bold", "middle")
text(
    ax + 328,
    ay + 196,
    "translation table chosen for each contig",
    11.5,
    "normal",
    "middle",
    MUTED,
)
arrow(ax + 328, ay + 206, ax + 328, ay + 238)
for i, x in enumerate(gene_xs(ax + 40)):
    rect(x, ay + 242, WIDTHS[i], 24, "#ffffff", MUTED, 12)
    text(x + WIDTHS[i] / 2, ay + 258, f"p{i + 1}", 11, "bold", "middle")
text(
    ax + 328,
    ay + 292,
    "Translated proteins (at least 30 aa by default)",
    12,
    "normal",
    "middle",
)
text(
    ax + 328,
    ay + 310,
    "Predictions are cached and shared with phu screen and phu jack",
    11.5,
    "normal",
    "middle",
    MUTED,
    "italic",
)

# ---- B: annotation ------------------------------------------------------
bx, by = 715, 64
panel(bx, by, 655, 330, "B", "Two HMM annotation tracks (PyHMMER)")
chip(bx + 400, by + 16, 62, "Pfam", PFAM)
chip(bx + 468, by + 16, 68, "KOfam", KOFAM)
chip(bx + 542, by + 16, 72, "PHROG", PHROG)
text(bx + 28, by + 68, "Relaxed track", 12.5, "bold")
text(
    bx + 138,
    by + 68,
    "E-value ≤ 1e-5; best hit per protein in each database",
    11.5,
    fill=MUTED,
)
hits_relaxed = [
    ["Pf", "KO", "PH"],
    ["PH"],
    ["Pf", "PH"],
    ["Pf", "KO"],
    ["KO", "PH"],
    [],
    ["Pf", "PH"],
]
hits_strict = [["Pf", "PH"], ["PH"], ["Pf"], ["Pf", "KO"], ["PH"], [], ["PH"]]
contig(bx + 28, by + 111, 600)
for i, x in enumerate(gene_xs(bx + 40)):
    tags(x + 2, by + 80, hits_relaxed[i])
    gene(x, by + 98, WIDTHS[i], label=f"g{i + 1}")
text(bx + 28, by + 146, "→ feeds V / VL scoring (panel C)", 11.5, "bold", fill=MUTED)

text(bx + 28, by + 180, "Strict track", 12.5, "bold")
text(
    bx + 120,
    by + 180,
    "Pfam GA cutoffs · KOfam KO thresholds · PHROG bit score ≥ 60, coverage ≥ 0.4",
    11.5,
    fill=MUTED,
)
contig(bx + 28, by + 223, 600)
for i, x in enumerate(gene_xs(bx + 40)):
    tags(x + 2, by + 192, hits_strict[i])
    gene(x, by + 210, WIDTHS[i], label=f"g{i + 1}")
text(
    bx + 28,
    by + 258,
    "→ feeds functional calls and reference matching (panel D)",
    11.5,
    "bold",
    fill=MUTED,
)
text(
    bx + 328,
    by + 302,
    "Pfam, KOfam and PHROG are scored and averaged independently; "
    "values are never mixed",
    11.5,
    "normal",
    "middle",
    MUTED,
    "italic",
)

# ---- C: scores ----------------------------------------------------------
cx, cy = 30, 414
panel(cx, cy, 655, 440, "C", "Viral-likeness scores (V-Score-Search tables)")
text(
    cx + 28,
    cy + 66,
    "Relaxed-track hits are looked up in the V / VL table, separately for each database",
    11.5,
    fill=MUTED,
)
VL = [4.6, 5.0, 4.2, 1.1, 4.8, 3.9, 4.4]
V = [10, 10, 9.6, 0.4, 10, 8.7, 10]
SHADES = {
    1.1: "#fdf2e9",
    3.9: "#b9a3d6",
    4.2: "#a88bcd",
    4.4: "#9a79c6",
    4.6: "#8a66bd",
    4.8: "#7b55b3",
    5.0: "#6b44a8",
}
text(cx + 28, cy + 96, "VL", 11, "bold", fill=MUTED)
text(cx + 28, cy + 158, "V", 11, "bold", fill=MUTED)
contig(cx + 48, cy + 121, 580)
xs = gene_xs(cx + 60)
for i, x in enumerate(xs):
    cand = i == 3
    mid = x + WIDTHS[i] / 2 - 3
    text(mid, cy + 96, f"{VL[i]:.1f}", 11.5, "bold", "middle", CAND if cand else INK)
    gene(
        x,
        cy + 108,
        WIDTHS[i],
        fill=SHADES[VL[i]],
        stroke=CAND if cand else MUTED,
        label=f"g{i + 1}",
        sw=2.6 if cand else 1.2,
        color=INK if cand else "#ffffff",
    )
    text(mid, cy + 158, f"{V[i]:g}", 11.5, "normal", "middle", CAND if cand else MUTED)
for first, last in ((0, 2), (4, 6)):
    left, right = xs[first], xs[last] + WIDTHS[last]
    # Keep clear of the label under the candidate gene.
    left, right = (left, right - 10) if first == 0 else (left + 10, right)
    out.append(
        f'<path d="M{left},{cy + 172} v8 H{right} v-8" fill="none" stroke="{PHROG}" '
        f'stroke-width="1.4"/>'
    )
    text(
        (left + right) / 2,
        cy + 196,
        "virus-like neighbours",
        11,
        "normal",
        "middle",
        PHROG,
    )
mid = xs[3] + WIDTHS[3] / 2 - 3
text(mid, cy + 182, "low V, low VL", 11, "bold", "middle", CAND)
text(mid, cy + 196, "in a viral context", 11, "bold", "middle", CAND)

rect(cx + 28, cy + 218, 600, 62, "#ffffff", LINE)
text(
    cx + 328,
    cy + 242,
    "contig AVL (per database) = mean VL of the contig's scored genes",
    12.5,
    "bold",
    "middle",
)
text(
    cx + 328,
    cy + 264,
    "here: (4.6 + 5.0 + 4.2 + 1.1 + 4.8 + 3.9 + 4.4) / 7 = 4.0",
    12,
    "normal",
    "middle",
    MUTED,
)

definitions = [
    ("V-score", "how exclusively a protein family occurs in viruses (10 = virus only)"),
    ("VL-score", "V-score weighted by how often the family is found in viruses"),
    ("AVL", "average VL across a contig, per database: the viral context of a gene"),
]
for i, (name, meaning) in enumerate(definitions):
    y = cy + 312 + i * 26
    text(cx + 28, y, name, 12, "bold")
    text(cx + 100, y, meaning, 12, fill=MUTED)
text(cx + 28, cy + 404, "VL shading:", 11, fill=MUTED)
for i, v in enumerate([3.9, 4.2, 4.4, 4.6, 4.8, 5.0]):
    rect(cx + 100 + i * 26, cy + 392, 24, 14, SHADES[v], "none", 2, 0)
text(cx + 100, cy + 420, "lower", 10, fill=MUTED)
text(cx + 100 + 6 * 26 - 2, cy + 420, "higher", 10, "normal", "end", MUTED)
rect(cx + 300, cy + 392, 24, 14, SHADES[1.1], CAND, 2, 1.6)
text(cx + 332, cy + 404, "below the gene cutoff", 11, fill=MUTED)

# ---- D: decision --------------------------------------------------------
dx, dy = 715, 414
panel(dx, dy, 655, 440, "D", "Candidate test and curation")
rect(dx + 28, dy + 54, 330, 106, "#ffffff", LINE)
text(dx + 42, dy + 76, "Candidate test, per database", 12.5, "bold")
for i, (lhs, rhs) in enumerate(
    [("contig AVL", "> 3"), ("gene VL", "< 3"), ("gene V", "< 10")]
):
    y = dy + 100 + i * 20
    text(dx + 42, y, lhs, 12)
    text(dx + 130, y, rhs, 12, "bold")
text(dx + 196, dy + 100, "Pfam, KOfam and", 11, fill=MUTED)
text(dx + 196, dy + 116, "PHROG are each", 11, fill=MUTED)
text(dx + 196, dy + 132, "tested on their own", 11, fill=MUTED)
arrow(dx + 360, dy + 106, dx + 468, dy + 106, WARN)
text(dx + 414, dy + 98, "fails in all", 10.5, "normal", "middle", WARN)
rect(dx + 472, dy + 92, 150, 28, "#ffffff", WARN, 14, 1.6, "5,3")
text(dx + 547, dy + 111, "not_candidate", 12, "bold", "middle", WARN)
arrow(dx + 110, dy + 162, dx + 110, dy + 186, CAND)
text(dx + 122, dy + 179, "passes in at least one database", 10.5, fill=CAND)
chip(dx + 28, dy + 190, 164, "AVG candidate", CAND, 28, 13)
arrow(dx + 110, dy + 220, dx + 110, dy + 236)

steps = [
    (
        "Match strict hits to the CheckAMG positive lists",
        "AMG: Pfam, KOfam · APG and AReG: Pfam, KOfam, PHROG",
        PHROG,
    ),
    (
        "AMG weight ≥ 0.6 (--min-amg-weight)",
        "applies to AMG evidence only",
        MUTED,
    ),
    (
        "Class filters (--filter-mode standard | strict | none)",
        "an essential or structural hit in any database, PHROG included, "
        "removes the class",
        PHROG,
    ),
]
for i, (head, detail, color) in enumerate(steps):
    y = dy + 240 + i * 50
    rect(dx + 28, y, 600, 42, "#ffffff", LINE)
    out.append(f'<circle cx="{dx + 48}" cy="{y + 21}" r="11" fill="{INK}"/>')
    text(dx + 48, y + 25, str(i + 1), 12, "bold", "middle", "#ffffff")
    text(dx + 68, y + 18, head, 12, "bold")
    text(dx + 68, y + 34, detail, 11, fill=color)

rect(dx + 28, dy + 390, 600, 40, "#ffffff", INK, 6, 1.4)
text(dx + 40, dy + 407, "Final status", 11.5, "bold")
text(
    dx + 40,
    dy + 423,
    "classified: putative_amg · putative_apg · putative_areg",
    10.5,
    fill=KOFAM,
)
text(dx + 356, dy + 407, "class_conflict · filtered · below_amg_weight", 10.5, fill=MUTED)
text(dx + 356, dy + 423, "unclassified_candidate · not_candidate", 10.5, fill=MUTED)

# ---- outputs ------------------------------------------------------------
oy = 872
text(30, oy + 32, "Outputs", 14, "bold")
files = [
    ("gene-annotations.tsv", "every predicted gene"),
    ("avg-candidates.tsv", "genes passing the test"),
    ("avg-predictions.tsv", "one final named class"),
    ("avg-evidence.tsv", "hits and reference matches"),
    (".phu/run.json", "parameters and database versions"),
]
for i, (name, what) in enumerate(files):
    x = 120 + i * 252
    rect(x, oy + 4, 240, 46, "#ffffff", LINE)
    text(x + 120, oy + 24, name, 12, "bold", "middle")
    text(x + 120, oy + 41, what, 10.5, "normal", "middle", MUTED)

text(
    W / 2,
    956,
    "All calls are computational candidates for prioritisation, "
    "not confirmed auxiliary genes.",
    11.5,
    "normal",
    "middle",
    MUTED,
    "italic",
)
text(
    W / 2,
    976,
    "Scores: Zhou et al. 2026 · reference tables: CheckAMG · terminology: "
    "Martin et al. 2025 · gene values shown are illustrative",
    10.5,
    "normal",
    "middle",
    MUTED,
)
out.append("</svg>")
print("\n".join(out))
