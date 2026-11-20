"""Draw the subsystem table version of Figure 3.

Run python plot_figure_3_table.py beside figure_3_seed_subsystems.tsv.
Python requirements: reportlab, pypdfium2 and Pillow. The PDF has searchable
text and landscape A4 pages. The PNG shows the same hierarchy in two columns.
"""

import argparse
import csv
import io
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from math import ceil
from pathlib import Path

import pypdfium2 as pdfium
import reportlab
from reportlab.lib.colors import HexColor
from reportlab.lib.pagesizes import A4, landscape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen.canvas import Canvas


REGULAR, BOLD = "TableVera", "TableVeraBold"
INK, MUTED, BAR, RULE = "#202C31", "#485B65", "#4682B4", "#CBD4D8"
PAGE_W, PAGE_H = landscape(A4)
MARGIN = 26
TABLE_WIDTH = PAGE_W - 2 * MARGIN
# Positions are in points from the table's left edge.
SUPER, SUPER_N, CLASS, CLASS_N, SUBSYSTEM, BAR_START, BAR_END = 0, 129, 143, 304, 320, 676, 758
BODY_TOP, FOOTER, LEADING = 96, 24, 11


@dataclass(frozen=True)
class Branch:
    superclass: str
    category: str
    leaves: tuple
    total: int


def read_branches(path):
    groups = defaultdict(dict)
    with path.open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            parent, category, leaf = (row[key].strip() for key in ("superclass", "class", "subsystem"))
            count = int(row["essential_genes"])
            leaves = groups[parent, category]
            if min(map(len, (parent, category, leaf))) == 0 or count < 1 or leaf in leaves:
                raise ValueError(f"Invalid or duplicate path: {parent} / {category} / {leaf}")
            leaves[leaf] = count
    if not groups:
        raise ValueError("No SEED assignments found.")
    totals = Counter()
    branches = []
    for (parent, category), leaves in groups.items():
        count = sum(leaves.values())
        totals[parent] += count
        ordered = tuple(sorted(leaves.items(), key=lambda item: (-item[1], item[0].casefold())))
        branches.append(Branch(parent, category, ordered, count))
    branches.sort(key=lambda b: (-totals[b.superclass], b.superclass.casefold(), -b.total, b.category.casefold()))
    return branches, totals


def superclass_name(label):
    return label.capitalize().replace("Dna", "DNA").replace("Rna", "RNA")


@lru_cache(maxsize=None)
def wrap(label, width, font, size):
    lines = []
    for word in label.split():
        candidate = f"{lines[-1]} {word}" if lines else word
        if lines and pdfmetrics.stringWidth(candidate, font, size) <= width:
            lines[-1] = candidate
        else:
            lines.append(word)
    if any(pdfmetrics.stringWidth(line, font, size) > width for line in lines):
        raise ValueError(f"Label exceeds the column width: {label}")
    return tuple(lines)


def group_rows(branches, first, last):
    """Lay out complete classes, with repeated superclass context at a break."""
    placed = []
    top = 0
    for i in range(first, last):
        group = branches[i]
        new_parent = i == first or branches[i - 1].superclass != group.superclass
        continued = i == first and i > 0 and branches[i - 1].superclass == group.superclass
        if i != first:
            top += 8 if new_parent else 5
        start = top
        super_lines = wrap(superclass_name(group.superclass), 103, BOLD, 9.6) if new_parent else ()
        class_lines = wrap(group.category, 130, REGULAR, 9.1)
        leaves = []
        for name, count in group.leaves:
            lines = wrap(name, BAR_START - SUBSYSTEM - 12, REGULAR, 9.15)
            leaves.append((lines, count, top))
            top += LEADING * len(lines) + 4
        parent_lines = len(super_lines) + int(continued)
        top = max(top, start + LEADING * max(parent_lines, len(class_lines)) + 4)
        placed.append((group, start, top, super_lines, continued, class_lines, leaves))
    return placed, top


def balanced_pages(branches, n, maximum=float("inf")):
    """Respect class boundaries and balance vertical space across all pages."""
    target = group_rows(branches, 0, len(branches))[1] / n

    @lru_cache(maxsize=None)
    def solve(first, remaining):
        if remaining == 0:
            return (0, ()) if first == len(branches) else (float("inf"), ())
        best = float("inf"), ()
        for end in range(first + 1, len(branches) - remaining + 2):
            height = group_rows(branches, first, end)[1]
            if height > maximum:
                break
            if remaining == 1 and end < len(branches):
                continue
            cost, rest = solve(end, remaining - 1)
            # Prefer boundaries between superclasses when column heights are close.
            split = end < len(branches) and branches[end - 1].superclass == branches[end].superclass
            proposal = (cost + (height - target) ** 2 + 900 * split, (end,) + rest)
            best = min(best, proposal)
        return best

    ends = solve(0, n)[1]
    return list(zip((0,) + ends[:-1], ends))


def text(canvas, x, y, value, size=9.15, font=REGULAR, color=INK, right=False):
    canvas.setFont(font, size)
    canvas.setFillColor(HexColor(color))
    draw = canvas.drawRightString if right else canvas.drawString
    draw(x, y, str(value))


def rule(canvas, x1, x2, y, width=0.45, color=RULE):
    canvas.setLineWidth(width)
    canvas.setStrokeColor(HexColor(color))
    canvas.line(x1, y, x2, y)


def title(canvas, x, y, branches, totals):
    text(canvas, x, y, "SEED functional assignments of essential genes", 18, BOLD)
    subtitle = (f"{sum(totals.values())} essential genes   ·   {len(totals)} superclasses   ·   "
                f"{len(branches)} classes   ·   {sum(len(g.leaves) for g in branches)} subsystem entries")
    text(canvas, x, y - 18, subtitle, 10.1, color=MUTED)


def table(canvas, x, y, branches, totals, interval, max_count):
    first, last = interval
    placed, height = group_rows(branches, first, last)
    for pos, word, right in [(SUPER, "Superclass", False), (SUPER_N, "Genes", True),
                             (CLASS, "Class", False), (CLASS_N, "Genes", True),
                             (SUBSYSTEM, "Subsystem", False)]:
        text(canvas, x + pos, y + 14, word, 9.3, BOLD, right=right)
    text(canvas, x + (BAR_START + BAR_END) / 2 - 50, y + 26, "Genes per subsystem", 7.3, BOLD)
    text(canvas, x + TABLE_WIDTH, y + 23, "Count", 7.6, BOLD, right=True)
    for fraction in (0, 0.5, 1):
        tick = round(max_count * fraction)
        tick_x = x + BAR_START + tick * (BAR_END - BAR_START) / max_count
        text(canvas, tick_x, y + 10, tick, 7, color=MUTED, right=tick == max_count)
        canvas.setStrokeColor(HexColor("#E4E9EB"))
        canvas.setLineWidth(0.35)
        canvas.line(tick_x, y + 5, tick_x, y - height)
    rule(canvas, x, x + TABLE_WIDTH, y + 4, 0.7, "#8FA2AA")
    for group, top, _, super_lines, continued, class_lines, leaves in placed:
        row_y = y - top - 8
        if top:
            rule(canvas, x if super_lines else x + CLASS, x + TABLE_WIDTH, y - top + 3)
        if super_lines:
            for j, line in enumerate(super_lines):
                text(canvas, x + SUPER, row_y - j * LEADING, line, 9.6, BOLD)
            if continued:
                text(canvas, x + SUPER, row_y - len(super_lines) * LEADING, "(continued)", 8, color=MUTED)
            text(canvas, x + SUPER_N, row_y, totals[group.superclass], 9.6, BOLD, right=True)
        for j, line in enumerate(class_lines):
            text(canvas, x + CLASS, row_y - j * LEADING, line, 9.1)
        text(canvas, x + CLASS_N, row_y, group.total, 9.1, right=True)
        for lines, count, offset in leaves:
            baseline = y - offset - 8
            for j, line in enumerate(lines):
                text(canvas, x + SUBSYSTEM, baseline - j * LEADING, line)
            canvas.setFillColor(HexColor(BAR))
            canvas.rect(x + BAR_START, baseline + 0.5, (BAR_END - BAR_START) * count / max_count,
                        5.7, fill=1, stroke=0)
            text(canvas, x + TABLE_WIDTH, baseline, count, 9.15, right=True)
    rule(canvas, x, x + TABLE_WIDTH, y - height - 2)


NOTE = "Unknown denotes genes without a SEED assignment in Table S1. Superclass and class counts are totals across their subsystems."


def draw_outputs(branches, totals, pdf_path, png_path, dpi):
    max_count = max(count for group in branches for _, count in group.leaves)
    usable = PAGE_H - BODY_TOP - FOOTER - 13
    min_pages = max(1, ceil(group_rows(branches, 0, len(branches))[1] / usable))
    for n in range(min_pages, len(branches) + 1):
        pages = balanced_pages(branches, n, usable)
        if pages:
            break
    else:
        raise ValueError("A complete class is too tall for the print page.")
    canvas = Canvas(str(pdf_path), pagesize=(PAGE_W, PAGE_H), invariant=True)
    canvas.setTitle("Figure 3 | SEED functional assignments of essential genes")
    for i, interval in enumerate(pages, 1):
        title(canvas, MARGIN, PAGE_H - 29, branches, totals)
        table(canvas, MARGIN, PAGE_H - BODY_TOP, branches, totals, interval, max_count)
        text(canvas, MARGIN, FOOTER, NOTE, 8.2, color=MUTED)
        text(canvas, PAGE_W - MARGIN, FOOTER, f"{i} / {len(pages)}", 8.3, color=MUTED, right=True)
        canvas.showPage()
    canvas.save()

    columns = balanced_pages(branches, min(2, len(branches)))
    max_height = max(group_rows(branches, *interval)[1] for interval in columns)
    gutter = 36
    width = 2 * MARGIN + len(columns) * TABLE_WIDTH + (len(columns) - 1) * gutter
    height = BODY_TOP + max_height + 34
    destination = io.BytesIO()
    overview = Canvas(destination, pagesize=(width, height), invariant=True)
    title(overview, MARGIN, height - 29, branches, totals)
    for i, interval in enumerate(columns):
        table(overview, MARGIN + i * (TABLE_WIDTH + gutter), height - BODY_TOP,
              branches, totals, interval, max_count)
    text(overview, MARGIN, 15, NOTE, 8.2, color=MUTED)
    overview.save()
    with pdfium.PdfDocument(destination.getvalue()) as document:
        page = document[0]
        bitmap = page.render(scale=dpi / 72)
        bitmap.to_pil().save(png_path, dpi=(dpi, dpi))
        bitmap.close()
        page.close()
    return pages, columns


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", type=Path, default=Path(__file__).with_name("figure_3_seed_subsystems.tsv"))
    parser.add_argument("--pdf", type=Path, default=Path("Figure_3_table.pdf"))
    parser.add_argument("--png", type=Path, default=Path("Figure_3_table.png"))
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error("--dpi must be positive.")
    fonts = Path(reportlab.__file__).parent / "fonts"
    for name, filename in [(REGULAR, "Vera.ttf"), (BOLD, "VeraBd.ttf")]:
        pdfmetrics.registerFont(TTFont(name, str(fonts / filename)))
    branches, totals = read_branches(args.data)
    draw_outputs(branches, totals, args.pdf, args.png, args.dpi)


if __name__ == "__main__":
    main()
