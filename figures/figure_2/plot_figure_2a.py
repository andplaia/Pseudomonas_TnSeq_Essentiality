"""Plot standalone panel A of Figure 2 for P. chlororaphis ATCC 9446.

Inputs are the original 6,206-locus TraDIS export and the supplemental
workbook (sheet "Table S1"). Requires Python 3.10+, matplotlib, numpy
and openpyxl. For example:

    python plot_figure_2a.py --sites sites.tsv --table-s1 supplement.xlsx

The original sites file is trimmed_reads_BADH-2_esen.out.CP144767.1.tradis_gene_insert_sites.csv
(SHA-256 408d9ac861e3c981ae9ef5802b6f2d89b920260bd2540adc7567b87d8cfbf8a5).
Use the "Table S1" worksheet in Supplementary tables (1-6).xlsx for the counts.

The upper graph uses the raw indices and original model specification. The
lower graph pairs the index and distinct-site count from each Table S1 row.
Both inputs use the full annotated span, including the ends.
The shaded interval extends from the lower limit of a two-sided 99.9% locus-
bootstrap interval to the reported essentiality cutoff. Pass --confidence 95
or --confidence 99 for the earlier levels. figure_2a_bootstrap.R regenerates the fit CSV.
Use --sensitivity-percent 10 to shade essential calls that would fall outside
the essential set if the cutoff were lowered by 10% of its reported value.
"""
import argparse
import csv
from math import lgamma
from pathlib import Path
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator
import numpy as np
from openpyxl import load_workbook

HIST_FILL = "#C4D8DB"
HIST_EDGE = "#819B9E"
LOW_FIT = "#C96918"
HIGH_FIT = "#176F9C"
CUTOFF = "#A72E34"
BAND = "#B9BABA"
POINT = "#416C78"
INK = "#202C31"


def read_table(path, cutoff):
    book = load_workbook(path, read_only=True, data_only=True)
    genes = {}
    try:
        rows = iter(book["Table S1"].values)
        next(rows)
        header = next(rows)
        cols = [header.index(k) for k in (
            "RefSeq Locus Tag", "Insertion index", "Insertion count", "Essentiality"
        )]
        for row in rows:
            if not any(cell is not None for cell in row):
                continue
            tag, index, count, call = [row[i] for i in cols]
            index, count = float(index), float(count)
            if (not tag or tag in genes or not np.isfinite(index) or index < 0
                    or not np.isfinite(count) or count < 0 or not count.is_integer()):
                raise ValueError(f"Invalid Table S1 values for {tag!r}.")
            if call not in ("Essential", "Non essential") or (index < cutoff) != (call == "Essential"):
                raise ValueError(f"Check the Table S1 classification and cutoff for {tag!r}.")
            genes[tag] = (index, count, call == "Essential")
    finally:
        book.close()
    if not genes:
        raise ValueError("No genes were found in Table S1.")
    indices, counts, essential = zip(*genes.values())
    return np.array(indices), np.array(counts, dtype=int), np.array(essential, dtype=bool), np.array(list(genes))


def read_fit(sites_path, fit_path, confidence=99.9):
    with fit_path.open(encoding="utf-8-sig") as handle:
        row, = csv.DictReader(handle)
    stats = {key: float(value) for key, value in row.items()}
    lower_key, upper_key = f"lower_{confidence:g}", f"upper_{confidence:g}"
    if lower_key not in stats or upper_key not in stats:
        raise ValueError(f"The fit summary lacks the {confidence:g}% interval. "
                         "Regenerate it with figure_2a_bootstrap.R or select an available confidence level.")
    lower, upper = stats[lower_key], stats[upper_key]
    with sites_path.open(encoding="utf-8-sig", newline="") as handle:
        first = handle.readline()
        handle.seek(0)
        rows = csv.DictReader(handle, delimiter="\t" if "\t" in first else ",")
        indices = np.array([float(r["ins_index"]) for r in rows])
    cutoff = round(stats["fitted_cutoff"], 5)
    if (len(indices) != stats["loci"] or not np.isfinite(indices).all()
            or np.any(indices < 0) or stats["successful"] != stats["replicates"]
            or cutoff != 0.03209 or not 0 < lower < cutoff < upper):
        raise ValueError("The input indices and the verified model summary are inconsistent.")
    bins = int(stats["histogram_bins"])
    width = stats["histogram_bin_width"]
    if not np.isclose(width * (bins - 1), indices.max()):
        raise ValueError("The fit histogram does not cover the observed range.")
    edges = (np.arange(bins + 1) - 0.5) * width
    density, _ = np.histogram(indices, bins=edges, density=True)
    x = np.linspace(0, indices.max(), int(stats["model_grid_points"]))
    rate = stats["exp_rate"]
    low = stats["exp_weight"] * rate * np.exp(-rate * x)
    shape, rate = stats["gamma_shape"], stats["gamma_rate"]
    if shape <= 1 or rate <= 0:
        raise ValueError("Invalid high-index gamma parameters.")
    high = np.zeros_like(x)
    high[1:] = stats["gamma_weight"] * np.exp(
        shape * np.log(rate) - lgamma(shape) + (shape - 1) * np.log(x[1:]) - rate * x[1:])
    return SimpleNamespace(bin_edges=edges, bin_density=density, insertion_index=x,
                           low_density=low, high_density=high, essential_cutoff=cutoff,
                           lower=lower, upper=upper, lower_95=stats["lower_95"],
                           confidence=confidence, loci=len(indices))


def render(model, indices, counts, essential, lower, output, sensitivity_percent=None):
    cutoff = model.essential_cutoff
    selected = essential & (indices >= lower) & (indices < cutoff)
    background = ~selected
    fig, (fit, points) = plt.subplots(2, 1, figsize=(7.8, 6.95), sharex=True,
                                      gridspec_kw={"height_ratios": [1.08, 1.0]})
    fig.subplots_adjust(left=0.145, right=0.97, bottom=0.18, top=0.954, hspace=0.18)
    for ax in (fit, points):
        ax.spines[["right", "top"]].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color("#65757B")
            ax.spines[side].set_linewidth(0.75)
        ax.tick_params(labelsize=9, color="#65757B", length=3.5)
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#E7EAEB", lw=0.65)
        ax.axvspan(lower, cutoff, facecolor=BAND, alpha=0.65, zorder=0.6)
        ax.axvline(cutoff, color=CUTOFF, lw=1.4, ls=(0, (4, 2.4)), zorder=5)

    fig.text(0.035, 0.965, "A", va="top", weight="bold", fontsize=17, color=INK)
    fit.bar(model.bin_edges[:-1], model.bin_density, np.diff(model.bin_edges),
            align="edge", color=HIST_FILL, edgecolor=HIST_EDGE, linewidth=0.44, zorder=1)
    low, = fit.plot(model.insertion_index, model.low_density, color=LOW_FIT, lw=2,
                    label="Low-index exponential fit", zorder=3)
    high, = fit.plot(model.insertion_index, model.high_density, color=HIGH_FIT, lw=2,
                     label="High-index gamma fit", zorder=3)
    fit.set(ylabel="Density", ylim=(0, max(model.low_density.max(), model.bin_density.max())*1.10))
    fit.yaxis.set_major_locator(FixedLocator([0, 5, 10, 15]))
    fit.tick_params(axis="x", bottom=False)
    fit.spines["bottom"].set_visible(False)
    fit.legend(handles=[Patch(facecolor=HIST_FILL, edgecolor=HIST_EDGE,
                              label=f"Observed indices ({model.loci:,} loci)"), low, high],
               loc="upper right", frameon=False, fontsize=9, labelcolor=INK)
    fit.annotate(f"Essentiality cutoff  {cutoff:.5f}",
                 xy=(cutoff, 0.89), xycoords=fit.get_xaxis_transform(),
                 xytext=(0.052, 0.89), textcoords=fit.get_xaxis_transform(),
                 va="center", fontsize=9, color=CUTOFF,
                 arrowprops={"arrowstyle": "-", "color": CUTOFF, "lw": 0.8})

    points.scatter(indices[background], counts[background], s=7, color=POINT, alpha=0.30,
                   linewidths=0, zorder=2, rasterized=True)
    points.scatter(indices[selected], counts[selected], s=13, color=INK, alpha=0.84,
                   linewidths=0.28, edgecolors="white", zorder=4, rasterized=True)
    # A square-root axis keeps zeros and large counts visible in the same view.
    points.set_yscale("function", functions=(lambda v: np.sqrt(np.maximum(v, 0)), np.square))
    points.set_ylim(0, max(100, counts.max()*1.10))
    points.set_yticks([v for v in (0, 25, 100, 250, 500, 1000, 1800) if v <= points.get_ylim()[1]])
    points.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:,.0f}"))
    points.set_xlabel("Insertion index", fontsize=11, labelpad=7)
    points.set_ylabel("Unique insertion sites per gene\n(square-root scale)", fontsize=10)
    points.set_xlim(min(model.bin_edges[0]-0.002, 0), max(model.bin_edges[-1], indices.max())+0.003)
    points.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.2f}"))
    # A magnified view makes the narrow validation band readable at print size.
    zoom = points.inset_axes([0.655, 0.62, 0.325, 0.34], zorder=6)
    # Use the earlier inset scale as a common reference for the validation bands.
    span = cutoff - model.lower_95
    left = min(cutoff - 2.3 * span, lower - 0.25 * span)
    right = cutoff + 1.25 * span
    nearby = (indices >= left) & (indices <= right)
    zoom.axvspan(lower, cutoff, color=BAND, alpha=0.65, lw=0)
    zoom.scatter(indices[nearby & background], counts[nearby & background],
                 s=10, color=POINT, alpha=0.59, linewidths=0)
    zoom.scatter(indices[selected], counts[selected], s=20, color=INK,
                 edgecolor="white", linewidths=0.4, zorder=3)
    zoom.axvline(cutoff, color=CUTOFF, lw=1.1, ls=(0, (4, 2.4)))
    zoom.set_xlim(left, right)
    zoom.set_ylim(-5, max(counts[nearby])*1.12)
    zoom.set_title("Near the cutoff", fontsize=8.3, pad=3, color=INK, loc="left")
    zoom.set_ylabel("Sites", fontsize=7.2, labelpad=2)
    zoom.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))
    zoom.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.3f}"))
    zoom.tick_params(labelsize=7, length=2.6, pad=1.6)
    for spine in zoom.spines.values():
        spine.set_color("#74848A")
        spine.set_linewidth(0.65)
    fig.legend(handles=[
        Line2D([], [], ls="", marker="o", mfc=POINT, mec="none", markersize=5,
               label=f"Table S1 genes (n = {len(indices):,})"),
        Patch(facecolor=BAND, label=f"Essential calls prioritized for validation (n = {selected.sum():,})")
    ], loc="lower center", bbox_to_anchor=(0.55, 0.060), frameon=False,
        ncol=2, columnspacing=1.7, fontsize=9, labelcolor=INK, handlelength=1.5)
    if sensitivity_percent is None:
        band_note = (f"Gray shading extends from the lower {model.confidence:g}% bootstrap limit ({lower:.6f}) "
                     f"to the essentiality cutoff ({cutoff:.5f}).")
    else:
        band_note = (f"Gray shading marks the follow-up range from {sensitivity_percent:g}% below the cutoff "
                     f"({lower:.6f}) to the cutoff ({cutoff:.5f}).")
    fig.text(0.55, 0.030,
             "Counts include all distinct insertion sites across each annotated gene, including both ends.\n"
             + band_note,
             ha="center", va="center", fontsize=8.1, color=INK)
    fig.align_ylabels()
    fig.savefig(output, dpi=600, facecolor="white")
    plt.close(fig)
    return int(selected.sum())


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--table-s1", type=Path, required=True)
    parser.add_argument("--fit", type=Path, default=Path(__file__).with_name("figure_2a_fit.csv"),
                        help="Fit and bootstrap summary from figure_2a_bootstrap.R")
    band = parser.add_mutually_exclusive_group()
    band.add_argument("--confidence", type=float, choices=(95, 99, 99.9), default=99.9,
                      help="Two-sided bootstrap confidence level for the validation band (default: 99.9)")
    band.add_argument("--sensitivity-percent", type=float, metavar="PERCENT",
                      help="Flag essential calls within this percentage below the reported cutoff")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    pct = args.sensitivity_percent
    if pct is not None and not 0 < pct < 100:
        parser.error("--sensitivity-percent must be greater than 0 and less than 100.")
    m = read_fit(args.sites, args.fit, args.confidence if pct is None else 95)
    x, y, ess, tags = read_table(args.table_s1, m.essential_cutoff)
    lower = m.lower if pct is None else m.essential_cutoff * (1 - pct / 100)
    default_name = "Figure_2A.png" if pct is None else f"Figure_2A_{pct:g}pct_sensitivity.png"
    n = render(m, x, y, ess, lower, args.output or Path(default_name), pct)
    print(f"{len(x):,} genes; {ess.sum()} essential; {n} flagged "
          f"({lower:.8f} <= insertion index < {m.essential_cutoff:.5f}).")
    print("Flagged loci:", ", ".join(tags[ess & (x >= lower)]))


if __name__ == "__main__":
    main()
