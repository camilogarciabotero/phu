"""Shared drawing helpers for the phu workflow figures (plain SVG output)."""

from xml.sax.saxutils import escape

INK, MUTED, LINE, PANEL = "#1f2933", "#52606d", "#9aa5b1", "#f7f9fb"
BLUE, GREEN, PURPLE = "#2f6fb0", "#2e8b57", "#8e44ad"
ORANGE, RED, GENE = "#d35400", "#c0392b", "#cbd5e0"
FONT = "Liberation Sans, Arial, Helvetica, DejaVu Sans, sans-serif"
W, H = 1400, 1000
out = []


def start(title):
    out.clear()
    out.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" '
        f'height="{H}" font-family="{FONT}">'
    )
    out.append("<defs>")
    for name, color in (("head", MUTED), ("head-o", ORANGE), ("head-r", RED)):
        out.append(
            f'<marker id="{name}" viewBox="0 0 10 10" refX="9" refY="5" '
            f'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{color}"/></marker>'
        )
    out.append("</defs>")
    rect(0, 0, W, H, "#ffffff", "none", 0, 0)
    text(W / 2, 40, title, 21, "bold", "middle")


def finish(*footer):
    for i, line in enumerate(footer):
        text(W / 2, 956 + i * 20, line, 11.5 - i, "normal", "middle", MUTED,
             "italic" if i == 0 else "normal")
    out.append("</svg>")
    print("\n".join(out))


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


def arrow(x1, y1, x2, y2, color=MUTED, sw=1.6, dash=None):
    marker = {MUTED: "head", ORANGE: "head-o", RED: "head-r"}[color]
    extra = f' stroke-dasharray="{dash}"' if dash else ""
    out.append(
        f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" '
        f'stroke-width="{sw}" marker-end="url(#{marker})"{extra}/>'
    )


def path(d, color=MUTED, sw=1.6, marker=True, dash=None):
    end = ' marker-end="url(#head)"' if marker else ""
    extra = f' stroke-dasharray="{dash}"' if dash else ""
    out.append(
        f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{sw}"{end}{extra}/>'
    )


def circle(x, y, r, fill, stroke="none", sw=1.2):
    out.append(
        f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}" stroke="{stroke}" '
        f'stroke-width="{sw}"/>'
    )


def gene(x, y, w, h=22, fill=GENE, stroke=MUTED, label="", sw=1.2, color=INK):
    tip = 9
    pts = (
        f"{x},{y} {x + w - tip},{y} {x + w},{y + h / 2} "
        f"{x + w - tip},{y + h} {x},{y + h}"
    )
    out.append(
        f'<polygon points="{pts}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>'
    )
    if label:
        text(x + (w - tip / 2) / 2, y + h / 2 + 4, label, 10.5, "bold", "middle", color)


def contig(x, y, width):
    out.append(
        f'<line x1="{x}" y1="{y}" x2="{x + width}" y2="{y}" stroke="{LINE}" '
        f'stroke-width="3"/>'
    )


def panel(x, y, w, h, letter, title):
    rect(x, y, w, h, PANEL, LINE, 10)
    circle(x + 24, y + 26, 13, INK)
    text(x + 24, y + 31, letter, 14, "bold", "middle", "#ffffff")
    text(x + 46, y + 31, title, 15, "bold")


def chip(x, y, w, label, color, h=20, size=11):
    rect(x, y, w, h, color, color, h / 2, 0)
    text(x + w / 2, y + h / 2 + 4, label, size, "bold", "middle", "#ffffff")


def pill(x, y, w, label, h=24, stroke=MUTED, size=11):
    rect(x, y, w, h, "#ffffff", stroke, h / 2)
    text(x + w / 2, y + h / 2 + 4, label, size, "bold", "middle")


def box(x, y, w, h, head, detail=None, stroke=LINE, head_size=12.5):
    rect(x, y, w, h, "#ffffff", stroke)
    if detail is None:
        text(x + w / 2, y + h / 2 + 4, head, head_size, "bold", "middle")
    else:
        text(x + w / 2, y + h / 2 - 3, head, head_size, "bold", "middle")
        text(x + w / 2, y + h / 2 + 14, detail, 11, "normal", "middle", MUTED)


def tag(x, y, label, color, w=24):
    rect(x, y, w, 13, color, color, 3, 0)
    text(x + w / 2, y + 10, label, 8.5, "bold", "middle", "#ffffff")


def outputs(files, y=872):
    text(30, y + 32, "Outputs", 14, "bold")
    step = 1260 / len(files)
    for i, (name, what) in enumerate(files):
        x = 120 + i * step
        rect(x, y + 4, step - 12, 46, "#ffffff", LINE)
        text(x + (step - 12) / 2, y + 24, name, 12, "bold", "middle")
        text(x + (step - 12) / 2, y + 41, what, 10.5, "normal", "middle", MUTED)


def decision_table(x, y, columns, colors, rows, rules, col_w=58, row_h=30):
    """Contigs x models presence dots, then one kept/dropped column per rule."""
    name_w = 86
    for j, (label, color) in enumerate(zip(columns, colors)):
        chip(x + name_w + j * col_w + 6, y, col_w - 12, label, color, 18, 10.5)
    rule_x = x + name_w + len(columns) * col_w + 40
    rule_w = 104
    for k, (label, _) in enumerate(rules):
        text(rule_x + k * rule_w + rule_w / 2, y + 13, label, 11.5, "bold", "middle")
    for i, (name, present) in enumerate(rows):
        ry = y + 30 + i * row_h
        rect(x, ry, name_w - 12, row_h - 8, "#ffffff", MUTED, 4)
        text(x + (name_w - 12) / 2, ry + 15, name, 11, "bold", "middle")
        for j, color in enumerate(colors):
            cx = x + name_w + j * col_w + col_w / 2
            if j in present:
                circle(cx, ry + 11, 7, color)
            else:
                circle(cx, ry + 11, 7, "#ffffff", LINE)
        for k, (_, keep) in enumerate(rules):
            kept = keep(present)
            cx = rule_x + k * rule_w + rule_w / 2
            rect(cx - 36, ry + 1, 72, 20, "#e8f5ee" if kept else "#fbeceb",
                 GREEN if kept else RED, 10, 1)
            text(cx, ry + 15, "kept" if kept else "dropped", 10.5, "bold", "middle",
                 GREEN if kept else RED)
    return y + 30 + len(rows) * row_h
