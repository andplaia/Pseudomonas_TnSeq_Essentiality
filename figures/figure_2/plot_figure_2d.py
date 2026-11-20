"""Draw Figure 2D from the accompanying gene and insertion-coordinate table.

Run with Python 3.10+ and Matplotlib:
    python plot_figure_2d.py
    python plot_figure_2d.py --data figure_2d_regions.tsv --output Figure_2D.png

The table uses 1-based inclusive coordinates on the ATCC 9446 chromosome
(CP144767.1 / NZ_CP144767.1). Gene intervals and strands come from RefSeq
annotation, protein-coding gene classifications from Table S1, and genomic
insertion positions from the archived insertions.tab. As in Figure 1, RNA
genes are gold and pseudogenes are charcoal. Each distinct genomic position contributes
one mark within a neighborhood, including positions outside annotated genes.
The three coordinate axes are independently scaled to the full flanking genes
of the original maps: ATP synthase 204,942–263,859 bp, cell wall
5,838,656–5,869,322 bp, and ribosome 6,319,901–6,360,897 bp. All overlapping
genes are drawn at their annotated coordinates and clipped at the window edge.
Selected gene names are printed.
"""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon
from matplotlib.ticker import FuncFormatter, MultipleLocator


REGIONS = {
    "atp": "ATP synthase and DNA replication",
    "wall": "Peptidoglycan synthesis and cell division",
    "ribosome": "Ribosomal proteins and RNA polymerase",
}
# Only these genes receive text labels. The arrows always include all genes.
LABELS = {
    "atp": {"glmS", "atpD", "atpA", "atpB", "dnaA", "dnaN", "gyrB", "glyS", "hemF"},
    "wall": {"lpxC", "ftsZ", "ftsA", "murC", "murG", "murD", "mraY", "ftsL", "mraZ"},
    "ribosome": {"uvrA", "rplQ", "secY", "rpsE", "rplB", "fusA", "rpoC", "rpoB", "nusG"},
}
INK = "#202C31"
MUTED = "#485B65"
SITE = "#364952"
COLORS = {"Essential": "#4682B4", "Nonessential": "#FF6347",
          "RNA genes": "#B98200", "Pseudogenes": "#41464B"}
PROTEIN_LABELS = {"Essential": "Essential", "Non essential": "Nonessential"}
RNA_BIOTYPES = {"tRNA", "rRNA", "tmRNA", "SRP_RNA", "RNase_P_RNA", "ncRNA"}


def feature_label(row):
    biotype = row["biotype"]
    if biotype == "pseudogene":
        return "Pseudogenes"
    if biotype in RNA_BIOTYPES:
        return "RNA genes"
    if biotype == "protein_coding" and row["essentiality"] in PROTEIN_LABELS:
        return PROTEIN_LABELS[row["essentiality"]]
    raise ValueError(f"Unknown biotype or essentiality for {row['locus_tag']}: {row}")


def read_regions(path):
    regions = {name: {"genes": [], "sites": []} for name in REGIONS}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            region = regions[row["region"]]
            start, end = int(row["start"]), int(row["end"])
            if not 1 <= start <= end:
                raise ValueError(f"Invalid coordinates: {row}")
            if row["kind"] == "region" and "bounds" not in region:
                region["bounds"] = (start, end)
            elif row["kind"] == "site" and start == end:
                region["sites"].append(start)
            elif row["kind"] == "gene" and row["strand"] in ("+", "-"):
                region["genes"].append((start, end, row["strand"], row["gene"], feature_label(row)))
            else:
                raise ValueError(f"Invalid record type, strand or Table S1 assignment: {row}")
    for key, region in regions.items():
        left, right = region["bounds"]
        sites = region["sites"]
        genes = sorted(region["genes"])
        if not sites or len(set(sites)) != len(sites) or min(sites) < left or max(sites) > right:
            raise ValueError(f"Duplicate, missing or out-of-range insertion positions in {key}.")
        if not genes or any(end < left or start > right for start, end, *_ in genes):
            raise ValueError(f"Missing genes or genes outside the window in {key}.")
        names = [name for _, _, _, name, _ in genes if name in LABELS[key]]
        if len(names) != len(LABELS[key]) or set(names) != LABELS[key]:
            raise ValueError(f"Selected labels do not uniquely match the genes in {key}.")
        region["genes"] = genes
    return regions


def gene_arrow(ax, start, end, strand, call, head_width):
    # Extend the inclusive bases to their edges; the arrowhead stays inside the gene.
    tail, tip = (start - 0.5, end + 0.5) if strand == "+" else (end + 0.5, start - 0.5)
    head = min(head_width, abs(tip - tail) * 0.48)
    shoulder = tip - head if strand == "+" else tip + head
    y, h = 0.52, 0.10
    points = [(tail, y - h), (shoulder, y - h), (tip, y), (shoulder, y + h), (tail, y + h)]
    ax.add_patch(Polygon(points, closed=True, facecolor=COLORS[call], edgecolor=INK,
                         linewidth=0.32, joinstyle="round", zorder=3))


def label_genes(ax, genes, selected, bp_per_pt):
    # Keep labels in genomic order while fanning out very close gene centers.
    # Text widths are measured at the final figure size so labels remain separated.
    entries = []
    for start, end, _, name, _ in genes:
        if name in selected:
            center = (start + end) / 2
            artist = ax.text(center, 0.82, name, ha="center", va="bottom", fontsize=8.8,
                             fontstyle="italic", color=INK, zorder=5)
            width = artist.get_window_extent().width / ax.get_window_extent().width
            entries.append((center, width * (ax.get_xlim()[1] - ax.get_xlim()[0]), artist))
    # A forward and a reverse pass give each dense cluster room on both sides.
    gap = 3.0 * bp_per_pt
    centers = [center for center, _, _ in entries]
    forward, backward = centers.copy(), centers.copy()
    for i in range(1, len(entries)):
        space = (entries[i - 1][1] + entries[i][1]) / 2 + gap
        forward[i] = max(centers[i], forward[i - 1] + space)
    for i in range(len(entries) - 2, -1, -1):
        space = (entries[i][1] + entries[i + 1][1]) / 2 + gap
        backward[i] = min(centers[i], backward[i + 1] - space)
    for i, (center, _, artist) in enumerate(entries):
        label_x = (forward[i] + backward[i]) / 2
        artist.set_x(label_x)
        ax.plot([center, center, label_x], [0.63, 0.69, 0.80], lw=0.48, color=MUTED, zorder=2)


def draw_region(ax, key, region):
    left, right = region["bounds"]
    ax.set(xlim=(left - 0.5, right + 0.5), ylim=(0, 1.0), yticks=[])
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.spines["bottom"].set_color("#A4B0B4")
    ax.spines["bottom"].set_linewidth(0.55)
    ax.tick_params(axis="x", labelsize=8.4, pad=3, length=2.5, width=0.55, color=MUTED)
    ax.xaxis.set_major_locator(MultipleLocator(10_000 if right - left > 50_000 else 5_000))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x / 1000:,.0f}"))
    axis_width_pt = ax.get_window_extent().width / ax.figure.dpi * 72
    bp_per_pt = (right - left + 1) / axis_width_pt
    for start, end, strand, _, call in region["genes"]:
        gene_arrow(ax, start, end, strand, call, head_width=3.5 * bp_per_pt)
    label_genes(ax, region["genes"], LABELS[key], bp_per_pt)
    ax.vlines(region["sites"], 0.08, 0.29, colors=SITE, linewidth=0.29, zorder=2)
    ax.text(0, 1.10, REGIONS[key], transform=ax.transAxes, color=INK, weight="bold",
            ha="left", va="bottom", fontsize=11.2)
    ax.text(1, 1.10, f"{len(region['sites']):,} sites", transform=ax.transAxes, color=MUTED,
            ha="right", va="bottom", fontsize=9.2)


def render(regions, output):
    fig = plt.figure(figsize=(7.5, 6.1), dpi=120, facecolor="white")
    fig.text(0.027, 0.965, "D", weight="bold", fontsize=16, va="top", color=INK)
    fig.text(0.085, 0.964, "Insertion sites in three genomic neighborhoods",
             weight="bold", fontsize=14, va="top", color=INK)
    handles = [Patch(facecolor=color, edgecolor=INK, linewidth=0.35, label=label)
               for label, color in COLORS.items()]
    handles.append(Line2D([], [], ls="", marker="|", markersize=9, color=SITE,
                          label="Distinct insertion site"))
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.08, 0.92), ncol=5,
               frameon=False, handlelength=1.25, columnspacing=1.4, handletextpad=0.45,
               fontsize=8.7, borderaxespad=0)
    for key, bottom in zip(REGIONS, (0.635, 0.380, 0.125)):
        ax = fig.add_axes([0.085, bottom, 0.88, 0.144])
        draw_region(ax, key, regions[key])
    fig.text(0.525, 0.047, "Chromosome position (kb)", ha="center", fontsize=10.8, color=INK)
    fig.savefig(output, dpi=600, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", type=Path, default=Path(__file__).with_name("figure_2d_regions.tsv"))
    parser.add_argument("--output", type=Path, default=Path("Figure_2D.png"))
    args = parser.parse_args()
    render(read_regions(args.data), args.output)
