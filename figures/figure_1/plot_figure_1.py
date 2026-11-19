"""Render the Pc9446 genome map with the gene-level track from Mapping_tnseq.ipynb.

Requires Python 3.10+, NumPy, Matplotlib and openpyxl. Example:
    python plot_figure_1.py --sites sites.tsv --gff Pc9446.gff3 --table-s1 supplement.xlsx --output Figure_1.png

The black trace uses the plotting export's insertion indices. The colored
essentiality bands use the assignments recorded in the Table S1 worksheet.
Records absent from that worksheet are drawn in a separate gray band.
"""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import numpy as np
from openpyxl import load_workbook


# Track radii are in plot units, from the center outward.
PERCENTILE = 95
CDS_PLUS = (0.61, 0.70, "#32CD32")
CDS_MINUS = (0.72, 0.81, "#800080")
CLASS_BANDS = {
    "Essential": (0.98, 1.05, "#4682B4"),
    "Non essential": (0.83, 0.92, "#FF6347"),
    "No S1 entry": (0.935, 0.965, "#686E73"),
}
DENSITY_RADII = (1.09, 1.39)
GENOME_RADIUS = 1.44


def read_assignments(path):
    workbook = load_workbook(path, read_only=True, data_only=True)
    assignments = {}
    try:
        rows = iter(workbook["Table S1"].values)
        next(rows)  # The first row is the table title.
        headers = next(rows)
        for values in rows:
            if not any(value is not None for value in values):
                continue
            row = dict(zip(headers, values))
            tag = row["Old Locus Tag"]
            label = row["Essentiality"]
            if not tag or tag in assignments or label not in ("Essential", "Non essential"):
                raise ValueError(f"Invalid or duplicate Table S1 assignment for {tag!r}.")
            coords = tuple(int(float(row[k])) for k in ("Gene start", "Gene end", "Strand"))
            assignments[tag] = (coords, label)
    finally:
        workbook.close()
    if not assignments:
        raise ValueError("Table S1 contains no assignments.")
    return assignments


def read_sites(path, assignments):
    seen = set()
    features = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        first = handle.readline()
        handle.seek(0)
        rows = csv.DictReader(handle, delimiter="\t" if "\t" in first else ",")
        for row in rows:
            tag = row["locus_tag"]
            if not tag or tag in seen:
                raise ValueError(f"Invalid or duplicate locus tag in plotting export: {tag!r}.")
            seen.add(tag)
            coords = tuple(int(row[k]) for k in ("start", "end", "strand"))
            label = "No S1 entry"
            if tag in assignments:
                expected, label = assignments[tag]
                if coords != expected:
                    raise ValueError(f"Genomic coordinates or strand disagree for {tag}.")
            features.append((coords[0], float(row["ins_index"]), label))
    if not features:
        raise ValueError("The sites file contains no features.")
    missing = assignments.keys() - seen
    if missing:
        raise ValueError(f"{len(missing)} Table S1 loci are missing from the plotting export.")
    features.sort(key=lambda entry: entry[0])
    starts, indices, labels = map(np.asarray, zip(*features))
    if not np.all(np.isfinite(indices)) or np.any(indices < 0):
        raise ValueError("Insertion indices must be finite and nonnegative.")
    return starts, indices, labels


def read_annotation(path):
    cds = {"+": [], "-": []}
    sequence_ids = set()
    length = 0
    with path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("##sequence-region "):
                _, name, _, end = line.split()
                sequence_ids.add(name)
                length = max(length, int(end))
            if line.startswith("#") or not line.strip():
                continue
            name, _, feature, start, end, _, strand, _, _ = line.rstrip("\r\n").split("\t")
            sequence_ids.add(name)
            length = max(length, int(end))
            if feature == "CDS":
                if strand not in cds:
                    raise ValueError(f"Unknown strand for CDS at {start}–{end}.")
                cds[strand].append((int(start), int(end)))
    if len(sequence_ids) != 1 or length <= 0 or not all(cds.values()):
        raise ValueError("Expected one chromosome with CDS annotations on both strands.")
    return length, cds


def angle(position, length):
    return 2 * np.pi * np.asarray(position, dtype=float) / length


def xy(theta, radius):
    return np.column_stack((radius * np.cos(theta), radius * np.sin(theta)))


def draw_band(ax, start_angle, end_angle, band):
    inner, outer, color = band
    segments = np.stack((xy(start_angle, inner), xy(end_angle, outer)), axis=1)
    ax.add_collection(LineCollection(segments, linewidths=0.32, colors=color, zorder=2))


def plot_map(sites, gff, table_s1, output):
    starts, indices, labels = read_sites(sites, read_assignments(table_s1))
    length, cds = read_annotation(gff)
    if np.any((starts < 1) | (starts > length)):
        raise ValueError("Feature coordinates fall outside the chromosome.")
    cap = np.percentile(indices, PERCENTILE)
    if cap <= 0:
        raise ValueError("The normalization percentile must be positive.")
    normalized = np.minimum(indices, cap) / cap
    theta = angle(starts, length)
    fig, ax = plt.subplots(figsize=(7.4, 7.5))
    fig.subplots_adjust(left=0.035, right=0.965, top=0.97, bottom=0.13)
    ax.set(aspect="equal", xlim=(-1.68, 1.74), ylim=(-1.60, 1.64))
    ax.axis("off")
    outline = np.linspace(angle(0, length), angle(length, length), 1000)
    low, high = DENSITY_RADII
    # Extend each guide horizontally from its top tangent. The keys are
    # physically aligned with the three radii, outside the chromosome.
    key_x = 1.29
    for level in (0, 0.5, 1):
        r = low + level * (high - low)
        ax.plot(*xy(outline, r).T, color="#B8BEC3", lw=0.5, ls=(0, (2, 3)), zorder=0)
        ax.plot([0, key_x], [r, r], color="#B8BEC3", lw=0.48, ls=(0, (2, 3)), zorder=0)
        ax.plot([key_x, key_x + 0.032], [r, r], color="#23282B", lw=0.75)
        ax.text(key_x + 0.052, r, f"{level:g}", va="center", fontsize=8.5)
    ax.plot([key_x, key_x], [low, high], color="#23282B", lw=0.75)
    ax.text(key_x, high + 0.075, "Normalized\ninsertion index", ha="center", va="bottom",
            fontsize=9, linespacing=1.05)
    ax.plot(*xy(theta, low + normalized * (high-low)).T, lw=0.43, color="#171A1B", zorder=3)
    for strand, band in [("+", CDS_PLUS), ("-", CDS_MINUS)]:
        begin, end = np.asarray(cds[strand]).T
        draw_band(ax, angle(begin, length), angle(end, length), band)
    for label, band in CLASS_BANDS.items():
        positions = theta[labels == label]
        draw_band(ax, positions, positions, band)
    ax.plot(*xy(outline, GENOME_RADIUS).T, color="#6D7478", lw=0.65)
    for pos in range(0, length, 1_000_000):
        a = angle(pos, length)
        ax.plot(*xy(a, np.array([GENOME_RADIUS, GENOME_RADIUS+0.04])).T, color="#555F65", lw=0.8)
        x, y = xy(a, GENOME_RADIUS+0.108)[0]
        ax.text(x, y, f"{pos//1_000_000} Mb", ha="left" if pos == 0 else "center",
                va="center", fontsize=8.5, bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.4})
    ax.text(0, 0.105, "Pseudomonas\nchlororaphis", ha="center", va="center",
            fontsize=15, fontstyle="italic", linespacing=1.06)
    ax.text(0, -0.12, "ATCC 9446", ha="center", va="center", fontsize=12)
    ax.text(0, -0.255, f"{length:,} bp", ha="center", va="center", fontsize=9, color="#525D62")
    class_handles = []
    for label, band in CLASS_BANDS.items():
        count = np.count_nonzero(labels == label)
        display = {"Essential": "Essential", "Non essential": "Nonessential",
                   "No S1 entry": "Absent from Table S1"}[label]
        class_handles.append(Line2D([], [], label=f"{display} ({count:,})", color=band[2], lw=3.5))
    strand_handles = [Line2D([], [], label=name, color=band[2], lw=3.5)
                      for name, band in [("CDS, + strand", CDS_PLUS), ("CDS, − strand", CDS_MINUS)]]
    for handles, height in [(class_handles, 0.070), (strand_handles, 0.025)]:
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, height),
                   ncol=len(handles), frameon=False, fontsize=9.5, columnspacing=1.6, handlelength=1.55)
    fig.savefig(output, dpi=400, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the Pc9446 genome map.")
    parser.add_argument("--sites", type=Path, required=True, help="BioTraDIS gene-level export")
    parser.add_argument("--gff", type=Path, required=True, help="Pc9446 GFF3 annotation")
    parser.add_argument("--table-s1", type=Path, required=True, help="Workbook with the Table S1 assignments")
    parser.add_argument("--output", type=Path, default=Path("Figure_1.png"))
    args = parser.parse_args()
    plot_map(args.sites, args.gff, args.table_s1, args.output)

