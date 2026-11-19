"""Render the Pc9446 genome map with the gene-level track from Mapping_tnseq.ipynb.

Requires Python 3.10+, NumPy and Matplotlib. Example:
    python plot_figure_1.py --sites <BioTraDIS gene export> --gff <Pc9446.gff3> --output Figure_1.png

The sites file may be tab- or comma-delimited and must contain start and
ins_index. The two index bands reproduce the original plot's 0.007 cutoff.
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


# Figure settings. Radii are in plot units. The origin gap holds the radial
# tick labels; it is purely graphical and omits no genomic positions.
GAP_DEGREES = 20
INDEX_CUTOFF = 0.007
PERCENTILE = 95
CDS_PLUS = (0.61, 0.70, "#32CD32")
CDS_MINUS = (0.72, 0.81, "#800080")
INDEX_HIGH = (0.83, 0.93, "#FF6347")
INDEX_LOW = (0.95, 1.05, "#4682B4")
DENSITY_RADII = (1.09, 1.39)
GENOME_RADIUS = 1.44


def read_sites(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        first = handle.readline()
        handle.seek(0)
        rows = csv.DictReader(handle, delimiter="\t" if "\t" in first else ",")
        features = sorted((int(row["start"]), float(row["ins_index"])) for row in rows)
    if not features:
        raise ValueError("The sites file contains no features.")
    starts, indices = np.asarray(features, dtype=float).T
    if not np.all(np.isfinite(indices)) or np.any(indices < 0):
        raise ValueError("Insertion indices must be finite and nonnegative.")
    return starts, indices


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
    fraction = np.asarray(position, dtype=float) / length
    return np.deg2rad(90 + GAP_DEGREES / 2 + (360 - GAP_DEGREES) * fraction)


def xy(theta, radius):
    return np.column_stack((radius * np.cos(theta), radius * np.sin(theta)))


def draw_band(ax, start_angle, end_angle, band):
    inner, outer, color = band
    segments = np.stack((xy(start_angle, inner), xy(end_angle, outer)), axis=1)
    ax.add_collection(LineCollection(segments, linewidths=0.32, colors=color, zorder=2))


def plot_map(sites, gff, output):
    starts, indices = read_sites(sites)
    length, cds = read_annotation(gff)
    if np.any((starts < 1) | (starts > length)):
        raise ValueError("Feature coordinates fall outside the chromosome.")
    cap = np.percentile(indices, PERCENTILE)
    if cap <= 0:
        raise ValueError("The normalization percentile must be positive.")
    normalized = np.minimum(indices, cap) / cap
    theta = angle(starts, length)
    fig, ax = plt.subplots(figsize=(7, 7.4))
    fig.subplots_adjust(left=0.045, right=0.955, top=0.965, bottom=0.115)
    ax.set(aspect="equal", xlim=(-1.64, 1.64), ylim=(-1.62, 1.66))
    ax.axis("off")
    outline = np.linspace(angle(0, length), angle(length, length), 1000)
    low, high = DENSITY_RADII
    for level in (0, 0.5, 1):
        r = low + level * (high - low)
        ax.plot(*xy(outline, r).T, color="#B8BEC3", lw=0.5, ls=(0, (2, 3)), zorder=0)
        ax.plot([0, 0.027], [r, r], color="#23282B", lw=0.8)
        ax.text(0.05, r, f"{level:g}", va="center", fontsize=8)
    ax.plot([0, 0], [low, high], color="#23282B", lw=0.8)
    ax.text(0.015, high + 0.13, "Normalized\ninsertion index", ha="center", va="bottom",
            fontsize=8.5, linespacing=1.05)
    ax.plot(*xy(theta, low + normalized * (high-low)).T, lw=0.43, color="#171A1B", zorder=3)
    for strand, band in [("+", CDS_PLUS), ("-", CDS_MINUS)]:
        begin, end = np.asarray(cds[strand]).T
        draw_band(ax, angle(begin, length), angle(end, length), band)
    below = indices < INDEX_CUTOFF
    draw_band(ax, theta[~below], theta[~below], INDEX_HIGH)
    draw_band(ax, theta[below], theta[below], INDEX_LOW)
    ax.plot(*xy(outline, GENOME_RADIUS).T, color="#6D7478", lw=0.65)
    for pos in range(0, length, 1_000_000):
        a = angle(pos, length)
        ax.plot(*xy(a, np.array([GENOME_RADIUS, GENOME_RADIUS+0.04])).T, color="#555F65", lw=0.8)
        x, y = xy(a, GENOME_RADIUS+0.108)[0]
        ax.text(x, y, f"{pos//1_000_000} Mb", ha="right" if pos == 0 else "center",
                va="center", fontsize=8.5)
    ax.text(0, 0.105, "Pseudomonas\nchlororaphis", ha="center", va="center",
            fontsize=15, fontstyle="italic", linespacing=1.06)
    ax.text(0, -0.12, "ATCC 9446", ha="center", va="center", fontsize=12)
    ax.text(0, -0.255, f"{length:,} bp", ha="center", va="center", fontsize=9, color="#525D62")
    bands = [
        (INDEX_LOW, f"Insertion index < {INDEX_CUTOFF:g}"),
        (INDEX_HIGH, f"Insertion index ≥ {INDEX_CUTOFF:g}"),
        (CDS_PLUS, "CDS, + strand"),
        (CDS_MINUS, "CDS, − strand"),
    ]
    handles = [Line2D([0], [0], label=label, color=band[2], lw=3.7, solid_capstyle="butt")
               for band, label in bands]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.024),
               ncol=2, frameon=False, fontsize=10, columnspacing=2.4, handlelength=1.8)
    fig.savefig(output, dpi=400, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the Pc9446 genome map.")
    parser.add_argument("--sites", type=Path, required=True, help="BioTraDIS gene-level export")
    parser.add_argument("--gff", type=Path, required=True, help="Pc9446 GFF3 annotation")
    parser.add_argument("--output", type=Path, default=Path("Figure_1.png"))
    args = parser.parse_args()
    plot_map(args.sites, args.gff, args.output)

