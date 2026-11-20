"""Draw Figure 4B from the distributed numerical results (Matplotlib required).

    python plot_figure_4b.py

Rebuild the numerical results from Figure 4A's aggregate counts and directional
Fisher tests (SciPy 1.11+ required), and draw the figure:

    python plot_figure_4b.py --from-4a figure_4a_cog_enrichment.tsv \
        --save-statistics figure_4b_cog_enrichment.tsv

Each essential or nonessential gene counts once per distinct COG letter and
once in the denominator for its group if annotated. For each category, the
two-by-two table has rows essential and nonessential, and columns assigned and
not assigned that letter. The point is the conditional maximum-likelihood odds
ratio. The interval is its two-sided, equal-tailed exact 95% confidence
interval conditional on the margins, without a pseudocount. The p and q values
come from the two directional Fisher tests per category in Figure 4A, with a
single Benjamini-Hochberg adjustment across all categories and directions.
The categories follow the order of the prevalence ratios in Figure 4A.
"""

import argparse
import csv
from math import ceil, floor, isfinite, log10
from pathlib import Path
from textwrap import wrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


INK = "#19323D"
MUTED = "#425861"
ESSENTIAL = "#147488"
NONESSENTIAL = "#A6532A"
NEUTRAL = "#64757C"
COUNTS = ("essential_in_category", "essential_annotated",
          "nonessential_in_category", "nonessential_annotated")
TESTS = ("p_essential_enriched", "p_nonessential_enriched",
         "q_essential_enriched", "q_nonessential_enriched")
OR = "conditional_odds_ratio_essential_vs_nonessential"
LOW = "exact_95ci_lower"
HIGH = "exact_95ci_upper"
COLUMNS = ("COG", "category", *COUNTS, OR, LOW, HIGH, *TESTS)


def read_table(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def conditional_results(source):
    from scipy.stats.contingency import odds_ratio

    output = []
    for src in source:
        e, et, n, nt = (int(src[c]) for c in COUNTS)
        result = odds_ratio([[e, et - e], [n, nt - n]], kind="conditional")
        interval = result.confidence_interval(confidence_level=0.95, alternative="two-sided")
        row = {field: src[field] for field in ("COG", "category", *COUNTS, *TESTS)}
        row.update({OR: result.statistic, LOW: interval.low, HIGH: interval.high})
        output.append(row)
    return output


def write_table(path, data):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def prepare(data):
    seen = set()
    denominators = set()
    entries = []
    for row in data:
        letter = row["COG"]
        e, et, n, nt = (int(row[c]) for c in COUNTS)
        value, low, high = (float(row[c]) for c in (OR, LOW, HIGH))
        p_e, p_n, q_e, q_n = (float(row[c]) for c in TESTS)
        if letter in seen or not (et > 0 and nt > 0 and 0 <= e <= et and 0 <= n <= nt):
            raise ValueError(f"Duplicate category or invalid gene count for {letter}")
        if not (0 <= low <= value <= high) or not all(0 <= p <= 1 for p in (p_e, p_n, q_e, q_n)):
            raise ValueError(f"Invalid odds ratio, interval or test result for {letter}")
        seen.add(letter)
        denominators.add((et, nt))
        prevalence_ratio = e * nt / (n * et) if n else float("inf")
        color, marker = NEUTRAL, "o"
        if q_e < 0.05:
            color, marker = ESSENTIAL, "o"
        elif q_n < 0.05:
            color, marker = NONESSENTIAL, "s"
        entries.append((letter, row["category"], value, low, high, color, marker, prevalence_ratio))
    if len(denominators) != 1:
        raise ValueError("All categories must have the same two annotated-gene denominators.")
    return sorted(entries, key=lambda v: (-v[-1], v[0])), denominators.pop()


def draw(data):
    entries, (essential_total, nonessential_total) = prepare(data)
    # Layout coordinates are typographical points, with a dedicated lane for
    # zero. Positive odds ratios have a true logarithmic coordinate.
    width = 582
    left, names_x = 12, 30
    zero_x, log_left, log_right = 254, 278, 436
    or_right, interval_right = 479, width - 12
    rendered = []
    y = 111
    for record in entries:
        label = "\n".join(wrap(record[1], 46, break_long_words=False))
        row_height = 16.4 if "\n" not in label else 24.5
        rendered.append((record, label, y + row_height / 2))
        y += row_height
    baseline, height = y, y + 69
    fig = plt.figure(figsize=(width / 72, height / 72))
    ax = fig.add_axes([0, 0, 1, 1], xlim=(0, width), ylim=(height, 0))
    ax.axis("off")

    def label(x, y, s, size, **kw):
        return ax.text(x, y, s, color=kw.pop("color", INK), va="center", family="DejaVu Sans",
                       fontsize=size, linespacing=1.06, **kw)

    def point(x, y, color, marker, size=4.8):
        ax.plot(x, y, marker=marker, ms=size, mfc="white" if color == NEUTRAL else color,
                mec=color, mew=1, linestyle="none", zorder=4)

    positives = [v for r in entries for v in r[2:5] if v > 0 and isfinite(v)]
    xmin = min(min(positives), 1) / 1.35
    xmax = max(max(positives), 1) * 1.15

    def x_pos(v):
        return log_left + (log10(v) - log10(xmin)) / (log10(xmax) - log10(xmin)) * (log_right - log_left)

    label(left, 16, "B", 15, weight="bold")
    label(names_x + 2, 16, "Odds ratios for COG category membership", 14, weight="bold")
    label(left, 37, f"COG-annotated genes    {essential_total:,} essential   ·   "
          f"{nonessential_total:,} nonessential", 9.7, color=MUTED)
    for x, color, marker, caption in [(left + 2, ESSENTIAL, "o", "Enriched in essential genes"),
                                      (173, NONESSENTIAL, "s", "Enriched in nonessential genes"),
                                      (361, NEUTRAL, "o", "q ≥ 0.05 in both directions")]:
        point(x, 57, color, marker)
        label(x + 8, 57, caption, 8.6)
    label(left, 84, "COG category", 10, weight="bold")
    label((zero_x + log_right) / 2, 80, "Odds ratio, essential / nonessential", 9.7,
          weight="bold", ha="center")
    label(or_right, 84, "OR", 9.6, weight="bold", ha="right")
    label(interval_right, 84, "95% exact CI", 9.6, weight="bold", ha="right")
    label(zero_x, 99, "0", 8, color=MUTED, ha="center")
    for power in range(ceil(log10(xmin)), floor(log10(xmax)) + 1):
        tick = 10.0 ** power
        x = x_pos(tick)
        label(x, 99, f"{tick:g}", 8.1, ha="center", color=INK if tick == 1 else MUTED,
              weight="bold" if tick == 1 else "normal")
        ax.plot([x, x], [107, baseline], color="#415660" if tick == 1 else "#D7E0E3",
                lw=1 if tick == 1 else 0.55, linestyle=(0, (3, 2)) if tick == 1 else "-",
                zorder=1)
    ax.plot([zero_x, zero_x], [107, baseline], lw=0.45, color="#DFE6E8", zorder=0)
    ax.plot([left, interval_right], [108, 108], lw=0.75, color="#9EAEB5")

    def format_or(v):
        return f"{v:.2f}" if isfinite(v) else "∞"

    for idx, (record, name, mid) in enumerate(rendered):
        code, _, estimate, low, high, color, marker, _ = record
        if idx and idx % 5 == 0:
            top = mid - (12.25 if "\n" in name else 8.2)
            ax.plot([left, interval_right], [top, top], lw=0.45, color="#DFE6E8", zorder=0)
        label(left, mid, code, 9.1, weight="bold", color=MUTED)
        label(names_x, mid, name, 9.1)
        lo, hi = x_pos(low) if low else log_left, x_pos(high) if isfinite(high) else log_right
        if low == 0:
            # The gap marks the discontinuity between zero and the positive log scale.
            ax.plot([zero_x, zero_x + 7], [mid, mid], lw=1.1, color=color)
        ax.plot([lo, hi], [mid, mid], lw=1.1, color=color, solid_capstyle="butt")
        if low:
            ax.plot([lo, lo], [mid - 2.4, mid + 2.4], lw=0.85, color=color)
        if isfinite(high):
            ax.plot([hi, hi], [mid - 2.4, mid + 2.4], lw=0.85, color=color)
        else:
            ax.plot(hi, mid, marker=">", ms=4, color=color)
        where = zero_x if estimate == 0 else x_pos(estimate) if isfinite(estimate) else log_right
        point(where, mid, color, marker)
        label(or_right, mid, format_or(estimate), 9.2, weight="bold", ha="right")
        label(interval_right, mid, f"{format_or(low)}–{format_or(high)}", 9, ha="right")

    ax.plot([left, interval_right], [baseline + 2, baseline + 2], lw=0.7, color="#9EAEB5")
    label(left, baseline + 15, "Points show conditional odds ratios with exact 95% confidence intervals. "
          "The dashed line marks equal odds (OR = 1).", 8, color=MUTED)
    label(left, baseline + 32, "Positive odds ratios use a logarithmic scale. The separate zero column "
          "shows the estimates for categories with no essential genes.", 8, color=MUTED)
    label(left, baseline + 49, f"Filled markers indicate q < 0.05 after Benjamini–Hochberg adjustment "
          f"across {2 * len(entries)} one-sided Fisher exact tests.", 8, color=MUTED)
    return fig


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", type=Path, default=here / "figure_4b_cog_enrichment.tsv")
    parser.add_argument("--from-4a", type=Path, help="Recalculate odds ratios and intervals using Figure 4A's results")
    parser.add_argument("--save-statistics", type=Path, help="Write the results calculated with --from-4a")
    parser.add_argument("--output", type=Path, default=Path("Figure_4B.png"))
    parser.add_argument("--pdf", type=Path, help="PDF destination (default: same stem as PNG)")
    args = parser.parse_args()
    if args.save_statistics and not args.from_4a:
        parser.error("--save-statistics requires --from-4a")
    data = conditional_results(read_table(args.from_4a)) if args.from_4a else read_table(args.data)
    fig = draw(data)
    if args.save_statistics:
        write_table(args.save_statistics, data)
    fig.savefig(args.output, dpi=600, facecolor="white")
    with matplotlib.rc_context({"pdf.fonttype": 42}):
        fig.savefig(args.pdf or args.output.with_suffix(".pdf"), facecolor="white",
                    metadata={"Title": "Odds ratios for COG category membership in Pc9446",
                              "CreationDate": None, "ModDate": None})
    plt.close(fig)


if __name__ == "__main__":
    main()
