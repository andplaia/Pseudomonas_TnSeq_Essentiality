"""Plot antiSMASH region 18 in Pc9446 from 1-based, inclusive coordinates.

    python plot_figure_5.py

To regenerate figure_5_region.tsv from the project archive, specify
--sources DIRECTORY containing CP144767.1.region018.gbk, the full RefSeq
P_chlororaphis_ATCC_9446.gff3, Supplementary tables (1-6).xlsx and insertions.tab.
Rebuilding the TSV requires openpyxl; the standalone plot requires Matplotlib.
The source region is numbered 18 by antiSMASH and BGC 17 in the consolidated
Table S2. The functional-role track reproduces the antiSMASH /gene_kind field.
"""

import argparse
import csv
import math
import re
from collections import Counter
from pathlib import Path
from urllib.parse import unquote

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Polygon, Rectangle
from matplotlib.textpath import TextToPath


FIELDS = ("kind", "start", "end", "strand", "refseq_tag", "original_tag", "gene",
          "product", "essentiality", "antismash_gene_kind")
COLORS = {"Essential": "#4682B4", "Nonessential": "#FF6347"}
ROLES = {"biosynthetic": "B", "biosynthetic-additional": "A", "regulatory": "R",
         "transport": "T", "other": "O"}
INK = "#19323D"
MUTED = "#475C65"
ROLE_INK = "#67457C"
ROLE_FILL = "#E2D5EB"
FONT = "DejaVu Sans"
LABEL_SIZE = 9.8
ANGLE = 53
METRICS = TextToPath()
HERE = Path(__file__).resolve().parent
CHROMOSOME = "NZ_CP144767.1"

# Concise labels for unnamed RefSeq CDSs. Both hypothetical proteins use their
# locus-tag suffixes; the accompanying table retains every full product and tag.
PRODUCT_LABELS = {
    "ethanolamine ammonia-lyase subunit EutB": "EutB subunit",
    "N-acetyltransferase": "N-acetyltransferase",
    "DedA family protein": "DedA family",
    "zinc-dependent peptidase": "Zinc peptidase",
    "addiction module antidote protein": "Antitoxin",
    "type II toxin-antitoxin system RelE/ParE family toxin": "RelE/ParE toxin",
    "DUF3828 domain-containing protein": "DUF3828 protein",
    "NAD(P)-dependent alcohol dehydrogenase": "Alcohol dehydrogenase",
    "helix-turn-helix transcriptional regulator": "HTH regulator",
    "exodeoxyribonuclease VII small subunit": "Exonuclease VII subunit",
    "TonB-dependent receptor domain-containing protein": "TonB receptor domain",
    "cobalamin-binding protein": "Cobalamin-binding protein",
    "ABC transporter substrate-binding protein": "ABC substrate-binding protein",
    "phosphatidylglycerophosphatase A": "Phosphatidylglycerophosphatase A",
}


def import_region(directory):
    from openpyxl import load_workbook

    region_file = (directory / "CP144767.1.region018.gbk").read_text()
    offset = int(re.search(r"Orig\. start\s*::\s*(\d+)", region_file)[1])
    right = int(re.search(r"Orig\. end\s*::\s*(\d+)", region_file)[1])
    length = int(re.search(r"^LOCUS\s+\S+\s+(\d+) bp", region_file)[1])
    left = offset + 1
    if right - offset != length:
        raise ValueError("The extraction offsets and GenBank sequence length disagree.")
    features = region_file.split("FEATURES             Location/Qualifiers", 1)[1].split("ORIGIN", 1)[0]
    roles = {}
    core = None
    for feature in re.split(r"(?m)^     (?=\w)", features):
        if feature.startswith("proto_core "):
            first, last = map(int, re.match(r"proto_core\s+(\d+)\.\.(\d+)", feature).groups())
            core = (offset + first, offset + last)
        if feature.startswith("CDS "):
            tag = re.search(r'/locus_tag="([^"]+)"', feature)[1]
            match = re.search(r'/gene_kind="([^"]+)"', feature)
            roles[tag] = match[1] if match else ""

    calls = {}
    book = load_workbook(directory / "Supplementary tables (1-6).xlsx", read_only=True, data_only=True)
    try:
        rows = iter(book["Table S1"].values)
        next(rows)
        cols = next(rows)
        indices = [cols.index(f) for f in ("RefSeq Locus Tag", "Old Locus Tag", "Essentiality")]
        for row in rows:
            tag, old, call = (row[i] for i in indices)
            if tag:
                calls[tag] = old, {"Essential": "Essential", "Non essential": "Nonessential"}[call]
    finally:
        book.close()

    genes = []
    with (directory / "P_chlororaphis_ATCC_9446.gff3").open(encoding="utf-8") as stream:
        for row in csv.reader(stream, delimiter="\t"):
            if len(row) != 9 or row[0] != CHROMOSOME or row[2] != "CDS":
                continue
            start, end = int(row[3]), int(row[4])
            if end < left or start > right:
                continue
            metadata = dict(pair.split("=", 1) for pair in row[8].split(";") if "=" in pair)
            tag = metadata["locus_tag"]
            old, call = calls[tag]
            genes.append(dict(kind="gene", start=start, end=end, strand=row[6],
                              refseq_tag=tag, original_tag=old, gene=metadata.get("gene", ""),
                              product=unquote(metadata["product"]), essentiality=call,
                              antismash_gene_kind=roles.get(old, "")))
    complete = {g["original_tag"] for g in genes if left <= g["start"] <= g["end"] <= right}
    predicted_core = [(g["start"], g["end"]) for g in genes if g["antismash_gene_kind"] == "biosynthetic"]
    if complete != roles.keys() or predicted_core != [core]:
        raise ValueError("Reference features disagree with the antiSMASH region or core annotation.")
    sites = set()
    with (directory / "insertions.tab").open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["scaffold"] == CHROMOSOME and left <= (pos := int(row["pos"])) <= right:
                sites.add(pos)
    return [dict(kind="region", start=left, end=right),
            *sorted(genes, key=lambda g: g["start"]),
            *(dict(kind="site", start=p, end=p) for p in sorted(sites))]


def read_map(rows):
    regions, genes, sites = [], [], []
    for row in rows:
        start, end = int(row["start"]), int(row["end"])
        if not 1 <= start <= end:
            raise ValueError(f"Invalid coordinate in {row}")
        if row["kind"] == "region":
            regions.append((start, end))
        elif row["kind"] == "gene" and row["strand"] in ("+", "-") and row["essentiality"] in COLORS:
            if row["antismash_gene_kind"] and row["antismash_gene_kind"] not in ROLES:
                raise ValueError(f"Unexpected antiSMASH gene kind in {row}")
            genes.append({**row, "start": start, "end": end})
        elif row["kind"] == "site" and start == end:
            sites.append(start)
        else:
            raise ValueError(f"Unknown map record: {row}")
    if len(regions) != 1:
        raise ValueError("Expected a single region.")
    left, right = regions[0]
    if not sites or len(sites) != len(set(sites)) or min(sites) < left or max(sites) > right:
        raise ValueError("Insertion coordinates are empty, repeated or outside the region.")
    if len(genes) != len({g["refseq_tag"] for g in genes}) or not genes:
        raise ValueError("Expected one record per overlapping gene.")
    if any(g["end"] < left or g["start"] > right for g in genes):
        raise ValueError("A gene lies entirely outside the map.")
    return (left, right), sorted(genes, key=lambda g: g["start"]), sorted(sites)


def display_name(g):
    if g["gene"]:
        return g["gene"]
    return PRODUCT_LABELS.get(g["product"], g["refseq_tag"].removeprefix("V4Y05_"))


def label_anchors(centers, minimum_spacing, left, right):
    # Separate anchors of parallel slanted labels without changing gene positions.
    if not centers:
        return []
    if (len(centers) - 1) * minimum_spacing > right - left:
        raise ValueError("Labels cannot fit; increase map width.")
    a = []
    for center in centers:
        a.append(max(center, a[-1] + minimum_spacing if a else left))
    b = [0.] * len(centers)
    for i in range(len(centers) - 1, -1, -1):
        b[i] = min(centers[i], b[i + 1] - minimum_spacing if i + 1 < len(centers) else right)
    return [max(left + i * minimum_spacing,
                min((x + y) / 2, right - (len(centers) - 1 - i) * minimum_spacing))
            for i, (x, y) in enumerate(zip(a, b))]


def render(region, genes, sites):
    left, right = region
    x0, x1, width = 30, 635, 650
    scale = (x1 - x0) / (right - left + 1)

    def x(bp):
        return x0 + (bp - (left - 0.5)) * scale

    def slant_height(g):
        font = FontProperties(family=FONT, size=LABEL_SIZE, style="italic" if g["gene"] else "normal")
        length = METRICS.get_text_width_height_descent(display_name(g), font, False)[0]
        return length * math.sin(math.radians(ANGLE)) + LABEL_SIZE

    plus, minus = ([g for g in genes if g["strand"] == strand] for strand in ("+", "-"))
    upper = 75 + max(map(slant_height, plus), default=0) + 16
    lower = upper + 19
    role_y = lower + 16 + max(map(slant_height, minus), default=0) + 22
    site_y, axis_y = role_y + 29, role_y + 40
    height = axis_y + 89
    fig = plt.figure(figsize=(width / 72, height / 72), dpi=100, facecolor="white")
    ax = fig.add_axes((0, 0, 1, 1), xlim=(0, width), ylim=(height, 0), aspect="equal")
    ax.set_axis_off()

    def text(xpos, ypos, content, size=9, **style):
        return ax.text(xpos, ypos, content, fontsize=size, family=FONT, va=style.pop("va", "center"),
                       color=style.pop("color", INK), **style)

    complete = [g for g in genes if left <= g["start"] <= g["end"] <= right]
    calls = Counter(g["essentiality"] for g in complete)
    text(12, 16, "Terpene-precursor region 18 in Pc9446", 18, weight="bold")
    text(12, 35, f"antiSMASH prediction   ·   CP144767.1   ·   {left:,}–{right:,} bp", 10.2, color=MUTED)
    text(12, 55, f"LB screen, {len(complete)} complete genes", 9.5, weight="bold")
    for xkey, call in [(194, "Essential"), (304, "Nonessential")]:
        ax.add_patch(Rectangle((xkey, 51.5), 11, 7, facecolor=COLORS[call], edgecolor=INK, lw=0.45))
        text(xkey + 15, 55, f"{call} ({calls[call]})", 9.3)
    ax.add_patch(Rectangle((439, 51.5), 11, 7, facecolor="white", hatch="///", edgecolor=INK, lw=0.45))
    text(454, 55, "Boundary-crossing gene", 9.3)

    for strand, group, y, up in [("+", plus, upper, -1), ("−", minus, lower, 1)]:
        text(10, y, strand, 12, weight="bold", color=MUTED)
        ax.plot([x0, x1], [y, y], color="#BBC7CA", lw=0.65, zorder=0)
        centers = []
        for g in group:
            xstart, xend = x(max(g["start"], left) - .5), x(min(g["end"], right) + .5)
            # Both boundary-crossing CDSs have their arrow tips inside this region.
            tail, tip = (xstart, xend) if g["strand"] == "+" else (xend, xstart)
            head = min(4.0, .48 * abs(tip - tail))
            shoulder = tip - math.copysign(head, tip - tail)
            h = 5.6
            pts = [(tail, y-h), (shoulder, y-h), (tip, y), (shoulder, y+h), (tail, y+h)]
            partial = g["start"] < left or g["end"] > right
            ax.add_patch(Polygon(pts, facecolor=COLORS[g["essentiality"]], edgecolor=INK,
                                 linewidth=.5, hatch="///" if partial else None,
                                 gid="gene-" + g["refseq_tag"]))
            centers.append((xstart+xend) / 2)
        anchors = label_anchors(centers, 14.5, x0 + 1, x1 - 1)
        for g, center, anchor in zip(group, centers, anchors):
            ax.plot([center, center, anchor], [y + up * 6, y + up * 9, y + up * 13],
                    color=MUTED, lw=.55)
            text(anchor, y + up * 14, display_name(g), LABEL_SIZE, ha="left" if up < 0 else "right",
                 va="bottom" if up < 0 else "top", rotation=ANGLE, rotation_mode="anchor",
                 style="italic" if g["gene"] else "normal", gid="label-" + g["refseq_tag"])

    text(x0, role_y-12, "antiSMASH gene roles", 9, weight="bold", color=ROLE_INK)
    ax.plot([x0, x1], [role_y, role_y], color="#D6CBDE", lw=.6)
    for g in genes:
        role = g["antismash_gene_kind"]
        if role:
            start, stop = x(g["start"]-.5), x(g["end"]+.5)
            ax.add_patch(Rectangle((start, role_y-4), stop-start, 8, facecolor=ROLE_FILL,
                                   edgecolor=ROLE_INK, lw=.50))
            text((start+stop)/2, role_y, ROLES[role], 7.4, ha="center", weight="bold", color=ROLE_INK)
    text(x0, site_y-12, f"Distinct insertion sites ({len(sites):,})", 9, weight="bold")
    ax.vlines([x(p) for p in sites], site_y-5, site_y+4, color="#304851", linewidth=.35,
              gid="occupied-positions")
    ax.plot([x0, x1], [axis_y, axis_y], color=MUTED, lw=.65)
    step = 2500
    for bp in range(math.ceil(left/step)*step, right+1, step):
        ax.plot([x(bp), x(bp)], [axis_y, axis_y+3.2], color=MUTED, lw=.6)
        text(x(bp), axis_y+11, f"{bp/1000:,.1f}", 8.8, ha="center")
    text((x0+x1)/2, axis_y+27, "Chromosome position (kb)", 10, ha="center")
    text(12, axis_y+47, "antiSMASH roles", 9.2, weight="bold", color=ROLE_INK)
    text(112, axis_y+47, "B  Biosynthetic core     A  Additional biosynthetic     R  Regulatory     T  Transport     O  Other",
         8.9, color=ROLE_INK)
    text(12, axis_y+63, "Blank role positions have no antiSMASH gene-role assignment. "
         "The predicted core is ispA (farnesyl diphosphate synthase).", 8.6, color=MUTED)
    text(12, axis_y+77, "Each insertion mark denotes one occupied chromosome coordinate. "
         "RS labels identify hypothetical proteins by RefSeq locus-tag suffix.", 8.6, color=MUTED)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", type=Path, default=HERE / "figure_5_region.tsv")
    parser.add_argument("--sources", type=Path, help="Rebuild --data from the original project files.")
    parser.add_argument("--output", type=Path, default=Path("Figure_5.png"))
    parser.add_argument("--pdf", type=Path, help="Vector PDF output (default same stem as PNG).")
    args = parser.parse_args()
    if args.sources:
        rows = import_region(args.sources)
        with args.data.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    else:
        with args.data.open(encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
    region, genes, sites = read_map(rows)
    fig = render(region, genes, sites)
    fig.savefig(args.output, dpi=600, facecolor="white")
    with matplotlib.rc_context({"pdf.fonttype": 42}):
        fig.savefig(args.pdf or args.output.with_suffix(".pdf"), facecolor="white",
                    metadata={"Title": "Terpene-precursor region 18 in Pc9446",
                              "CreationDate": None, "ModDate": None})
    plt.close(fig)


if __name__ == "__main__":
    main()
