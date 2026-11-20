"""Draw a standalone gene-position profile from distinct insertion sites.

Figure 2B uses essential genes and Figure 2C uses nonessential genes, as
classified in Table S1. In poolffff.xlsx, locusId identifies the assigned gene
and f is its 5'-to-3' fractional coordinate. The script recovers the genomic
coordinate from f and the annotated endpoints. Each distinct (gene, coordinate)
pair contributes once, regardless of the number of records or reads. A site
assigned to two overlapping genes can therefore contribute to both profiles
of gene-relative position.

Each bin spans two percentage points of the inclusive coding sequence (first
base = 0%, last base = 100%). Internal bin boundaries belong to the higher
bin; 100% belongs to the final bin. Enrichment is the pooled number of sites
in that bin divided by the pooled mean number per bin.

Reproduce the approved panel using the distributed counts:
    python plot_figure_2_sites.py --panel B --counts figure_2b_sites.csv

Recompute the counts from the workbooks:
    python plot_figure_2_sites.py --panel B --pool poolffff.xlsx \
        --table-s1 "Supplementary tables (1-6).xlsx" --save-counts figure_2b_sites.csv
To draw the standalone C review, use --panel C with the same workbooks.
Requires Python 3.10+, NumPy and Matplotlib. Workbook input also needs openpyxl.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
import numpy as np


BINS = 50
CALL = {"B": "Essential", "C": "Non essential"}
TITLE = {"B": "Insertion sites in essential genes", "C": "Insertion sites in nonessential genes"}
DEFAULT_FILE = {"B": "Figure_2B.png", "C": "Figure_2C_review.png"}
COLUMNS = ("panel", "genes_with_sites", "start_percent", "end_percent", "site_count")
TERMINAL = "#CA752C"
CENTRAL = "#367C99"
INK = "#202C31"


def read_gene_spans(path, panel):
    from openpyxl import load_workbook

    book = load_workbook(path, read_only=True, data_only=True)
    spans = {}
    seen = set()
    try:
        rows = iter(book["Table S1"].values)
        next(rows)  # sheet title
        header = next(rows)
        cols = [header.index(key) for key in
                ("RefSeq Locus Tag", "Gene start", "Gene end", "Strand", "Essentiality")]
        for row in rows:
            if not any(value is not None for value in row):
                continue
            tag, start, end, strand, call = (row[i] for i in cols)
            if not tag or tag in seen or call not in CALL.values():
                raise ValueError(f"Invalid gene or assignment in Table S1: {tag!r}, {call!r}.")
            seen.add(tag)
            if call == CALL[panel]:
                if int(start) != start or int(end) != end or start < 1 or start >= end or strand not in (-1, 1):
                    raise ValueError(f"Invalid coding sequence coordinates for {tag}.")
                spans[tag] = (int(start), int(end), int(strand))
    finally:
        book.close()
    if not spans:
        raise ValueError(f"Table S1 contains no {CALL[panel]} calls.")
    return spans


def count_pool_sites(path, spans):
    from openpyxl import load_workbook

    book = load_workbook(path, read_only=True, data_only=True)
    sites = defaultdict(set)
    records = 0
    try:
        rows = iter(book["Sheet1"].values)
        header = next(rows)
        tag_col, frac_col = (header.index(key) for key in ("locusId", "f"))
        for row in rows:
            tag = row[tag_col]
            if tag not in spans:
                continue
            f = row[frac_col]
            if not isinstance(f, (int, float)) or not 0 <= f <= 1:
                raise ValueError(f"Invalid 5'-to-3' relative position for {tag}: {f!r}.")
            start, end, strand = spans[tag]
            offset = f * (end - start)
            bases = round(offset)
            if abs(offset - bases) > 1e-5:
                raise ValueError(f"The relative position for {tag} does not resolve to a genomic base.")
            position = start + bases if strand == 1 else end - bases
            sites[tag].add(position)
            records += 1
    finally:
        book.close()
    if not sites:
        raise ValueError("No pool records matched the selected Table S1 locus tags.")
    counts = [0] * BINS
    for tag, positions in sites.items():
        start, end, strand = spans[tag]
        for position in positions:
            offset = position - start if strand == 1 else end - position
            bin_index = min(BINS - 1, (BINS * offset) // (end - start))
            counts[bin_index] += 1
    print(f"Matched {records:,} pool records from {len(sites):,} genes.")
    return len(sites), counts


def write_counts(path, panel, genes, counts):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(COLUMNS)
        for i, count in enumerate(counts):
            writer.writerow((panel, genes, 2 * i, 2 * (i + 1), count))


def read_counts(path, panel):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != BINS or not set(COLUMNS).issubset(rows[0]):
        raise ValueError("Expected 50 bins with the columns used by --save-counts.")
    genes = int(rows[0]["genes_with_sites"])
    counts = []
    for i, row in enumerate(rows):
        count = int(row["site_count"])
        if (row["panel"] != panel or int(row["genes_with_sites"]) != genes or count < 0
                or (int(row["start_percent"]), int(row["end_percent"])) != (2*i, 2*(i+1))):
            raise ValueError(f"Invalid panel, cohort, bin edges, or site count in row {i+2}.")
        counts.append(count)
    if genes < 1 or sum(counts) < genes:
        raise ValueError("Each included gene must have at least one insertion site.")
    return genes, counts


def render(panel, genes, counts, output):
    counts = np.asarray(counts)
    enrichment = counts / np.mean(counts)
    fig, ax = plt.subplots(figsize=(6.6, 4.65))
    fig.subplots_adjust(left=0.13, right=0.966, bottom=0.24, top=0.82)
    fig.text(0.035, 0.956, panel, weight="bold", size=17, color=INK, va="top")
    fig.text(0.13, 0.951, TITLE[panel], size=14, weight="semibold", color=INK, va="top")
    fig.text(0.13, 0.883, f"{genes:,} genes with mapped sites   ·   {counts.sum():,} distinct sites",
             size=9, color="#45565D", va="top")
    colors = [TERMINAL if i < 5 or i >= 45 else CENTRAL for i in range(BINS)]
    ax.bar(np.arange(BINS) * 2, enrichment, width=2, align="edge", color=colors,
           linewidth=0.45, edgecolor="white", zorder=2)
    ax.axhline(1, color=INK, lw=0.9, ls=(0, (3, 2)), zorder=3)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, max(2, enrichment.max() * 1.075))
    ax.set_ylabel("Insertion-site enrichment (fold)", fontsize=11)
    ax.set_xlabel("Position along coding sequence (5′ to 3′)", fontsize=11, labelpad=6)
    ax.set_xticks(range(0, 101, 20))
    ax.xaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    ax.tick_params(labelsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#E4E8E8", lw=0.65)
    ax.legend(handles=[Patch(facecolor=TERMINAL, label="Terminal 10% at each end"),
                       Patch(facecolor=CENTRAL, label="Central 80%")],
              frameon=False, fontsize=9, loc="upper left")
    fig.text(0.55, 0.055, "Each bin spans 2% of gene length. The dashed line marks the mean sites per bin (1×).",
             ha="center", fontsize=8.2, color=INK)
    fig.savefig(output, dpi=600, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--panel", required=True, choices=("B", "C"))
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--counts", type=Path, help="Previously exported 50-bin site counts")
    inputs.add_argument("--pool", type=Path, help="Positional workbook with locusId and f columns")
    parser.add_argument("--table-s1", type=Path, help="Supplemental workbook required with --pool")
    parser.add_argument("--save-counts", type=Path, help="Write a shareable bin-count CSV")
    parser.add_argument("--output", type=Path, help="Path to the standalone PNG")
    args = parser.parse_args()
    if bool(args.pool) != bool(args.table_s1):
        parser.error("Use --pool and --table-s1 together, or provide --counts alone.")
    if args.pool:
        genes, counts = count_pool_sites(args.pool, read_gene_spans(args.table_s1, args.panel))
    else:
        genes, counts = read_counts(args.counts, args.panel)
    if args.save_counts:
        write_counts(args.save_counts, args.panel, genes, counts)
    render(args.panel, genes, counts, args.output or DEFAULT_FILE[args.panel])
    print(f"{sum(counts):,} distinct sites; first 10%: {sum(counts[:5]):,}; "
          f"last 10%: {sum(counts[-5:]):,}; final 2%: {counts[-1]:,}.")


if __name__ == "__main__":
    main()
