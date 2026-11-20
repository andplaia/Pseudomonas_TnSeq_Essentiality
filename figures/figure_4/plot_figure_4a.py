"""Plot COG category enrichment using essentiality calls and eggNOG-mapper assignments.

Run beside the distributed TSV:
    python plot_figure_4a.py

To export numerical enrichments and adjusted Fisher exact test results
(requires SciPy 1.11 or later):
    python plot_figure_4a.py --save-statistics figure_4a_cog_enrichment.tsv

To reconstruct the TSV, run with the original workbooks:
    python plot_figure_4a.py --cohort "Supplementary tables (1-6).xlsx" \
        --eggnog anot_eggNog_Pchloro.xlsx --save-data figure_4a_cog_patterns.tsv

The cohort worksheet identifies each essential or nonessential locus by its
RefSeq and old locus tags. The latter is the query identifier in the original
eggNOG-mapper export. The distributed TSV contains the number of genes with
each exact annotation pattern in each group, including missing assignments.
One gene contributes once to every distinct annotated category. Category
prevalence is the number of genes in a category divided by the number with
any COG assignment in the same essentiality group. Fold enrichment is the
essential prevalence divided by the nonessential prevalence. Both groups'
missing assignments are counted separately in the figure note.
"""

import argparse
import csv
from collections import Counter
from math import ceil
from pathlib import Path
from textwrap import wrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Functional category names follow the NCBI COG classification.
# https://www.ncbi.nlm.nih.gov/research/cog/
NAMES = {
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "M": "Cell wall, membrane and envelope biogenesis",
    "N": "Cell motility",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
    "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "X": "Mobilome: prophages, transposons",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
}

FIELDS = ("essentiality", "cog_pattern", "annotation_status", "gene_count")
ASSIGNED = "assigned"
NO_CATEGORY = "no_category_in_export"
NO_QUERY = "query_absent_from_export"
INK = "#19323D"
MUTED = "#425861"
E_COLOR = "#147488"
NE_COLOR = "#A6532A"
GROUPS = ("essential", "nonessential")


def from_workbooks(cohort_path, eggnog_path):
    from openpyxl import load_workbook

    book = load_workbook(eggnog_path, data_only=True, read_only=True)
    try:
        rows = iter(book["Sheet1"].values)
        header = next(rows)
        query_col, cog_col = (header.index(k) for k in ("query", "COG_category"))
        annotations = {}
        for row in rows:
            query = row[query_col]
            if query is None and all(value is None for value in row):
                continue
            if not query or query in annotations:
                raise ValueError(f"Missing or duplicate eggNOG query: {query!r}")
            code = str(row[cog_col] or "").strip()
            annotations[query] = "" if code == "-" else code
    finally:
        book.close()

    patterns = Counter()
    mismatches = Counter()
    seen_refseq, seen_old = set(), set()
    book = load_workbook(cohort_path, data_only=True, read_only=True)
    try:
        rows = iter(book["Table S1"].values)
        next(rows)  # worksheet title
        header = next(rows)
        ref_col, old_col, call_col, compiled_col = (header.index(k) for k in
            ("RefSeq Locus Tag", "Old Locus Tag", "Essentiality", "COG"))
        for row in rows:
            if all(value is None for value in row):
                continue
            call = {"Essential": "essential", "Non essential": "nonessential"}[row[call_col]]
            ref, old = row[ref_col], row[old_col]
            if not ref or not old or ref in seen_refseq or old in seen_old:
                raise ValueError(f"Missing or duplicate locus tag: {ref!r}, {old!r}")
            seen_refseq.add(ref)
            seen_old.add(old)
            code = annotations.get(old, "")
            if old not in annotations:
                status = NO_QUERY
            else:
                status = ASSIGNED if code else NO_CATEGORY
            compiled = str(row[compiled_col] or "").strip()
            if code != ("" if compiled == "-" else compiled):
                mismatches[call] += 1
            patterns[call, code, status] += 1
    finally:
        book.close()
    print("Differences from the compiled COG column by group:", dict(mismatches))
    return [dict(zip(FIELDS, (*pattern, n))) for pattern, n in sorted(patterns.items())]


def summarize(rows):
    counts = {group: Counter() for group in GROUPS}
    statuses = {group: Counter() for group in GROUPS}
    seen = set()
    multi = Counter()
    for row in rows:
        group, code, status, n = (row[key] for key in FIELDS)
        n = int(n)
        key = group, code, status
        if n < 1 or key in seen:
            raise ValueError(f"Invalid count or duplicate input pattern: {row!r}")
        seen.add(key)
        if group not in GROUPS or status not in {ASSIGNED, NO_CATEGORY, NO_QUERY} or bool(code) != (status == ASSIGNED):
            raise ValueError(f"Invalid group or annotation status: {row!r}")
        categories = set(code)
        if not categories.issubset(NAMES):
            raise ValueError(f"Unknown COG category: {code!r}")
        for category in categories:
            counts[group][category] += n
        statuses[group][status] += n
        if len(categories) > 1:
            multi[group] += n
    if any(not statuses[group][ASSIGNED] for group in GROUPS):
        raise ValueError("Each essentiality group must include at least one COG-annotated gene.")
    return counts, statuses, multi


def category_rows(counts, statuses):
    entries = []
    for letter in set(counts["essential"]) | set(counts["nonessential"]):
        e, ne = (counts[group][letter] for group in GROUPS)
        e_pct = 100 * e / statuses["essential"][ASSIGNED]
        ne_pct = 100 * ne / statuses["nonessential"][ASSIGNED]
        fold = e_pct / ne_pct if ne_pct else float("inf")
        entries.append((letter, e, ne, e_pct, ne_pct, fold))
    return sorted(entries, key=lambda item: (-item[-1], item[0]))


def write_statistics(path, counts, statuses):
    """Apply BH correction across both one-sided Fisher tests for every observed category."""
    from scipy.stats import false_discovery_control, fisher_exact

    e_total = statuses["essential"][ASSIGNED]
    ne_total = statuses["nonessential"][ASSIGNED]
    records, pvalues = [], []
    for code, e, ne, _, _, fold in sorted(category_rows(counts, statuses)):
        contingency = [[e, e_total - e], [ne, ne_total - ne]]
        p_e = fisher_exact(contingency, alternative="greater").pvalue
        p_ne = fisher_exact(contingency, alternative="less").pvalue
        records.append([code, NAMES[code], e, e_total, ne, ne_total, fold, p_e, p_ne])
        pvalues.extend([p_e, p_ne])
    adjusted = false_discovery_control(pvalues, method="bh")
    columns = ("COG", "category", "essential_in_category", "essential_annotated",
               "nonessential_in_category", "nonessential_annotated", "fold_essential_vs_nonessential",
               "p_essential_enriched", "p_nonessential_enriched", "q_essential_enriched",
               "q_nonessential_enriched")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for i, record in enumerate(records):
            writer.writerow([*record, adjusted[2 * i], adjusted[2 * i + 1]])


def draw(counts, statuses, multi):
    # The layout uses points so type and marker sizes have physical dimensions.
    width = 554
    left, label_x, plot_x, plot_end = 12, 30, 256, 405
    e_count_x, ne_count_x, fold_x = 443, 483, 541
    entries, cursor = [], 92
    tested = category_rows(counts, statuses)
    for row in tested:
        if row[0] in ("A", "B"):
            continue
        words = wrap(NAMES[row[0]], width=48, break_long_words=False)
        row_height = 16.4 if len(words) == 1 else 24.5
        entries.append((row, "\n".join(words), cursor + row_height / 2))
        cursor += row_height
    baseline, height = cursor, cursor + 94
    fig = plt.figure(figsize=(width / 72, height / 72))
    ax = fig.add_axes([0, 0, 1, 1], xlim=(0, width), ylim=(height, 0))
    ax.axis("off")

    def label(x, y, s, size, **kw):
        ax.text(x, y, s, va="center", color=kw.pop("color", INK), fontsize=size,
                family="DejaVu Sans", linespacing=1.06, **kw)

    label(left, 16, "A", 15, weight="bold")
    label(label_x + 2, 16, "COG category enrichment in essential genes", 14, weight="bold")
    label(left, 39, "COG-annotated genes", 9.7, color=MUTED)
    for x, group, name, marker, color in [(131, "essential", "Essential (E)", "o", E_COLOR),
                                           (308, "nonessential", "Nonessential (NE)", "s", NE_COLOR)]:
        ax.plot(x, 39, marker=marker, ms=5.5, mec=color, mfc=color, ls="none")
        label(x + 8, 39, f"{name}, n = {statuses[group][ASSIGNED]:,}", 9.6)
    label(left, 66, "COG category", 10.2, weight="bold")
    label(plot_x, 63, "Genes in category (%)", 10, weight="bold")
    label((e_count_x + ne_count_x) / 2, 63, "Gene count", 9.7, weight="bold", ha="center")
    label(e_count_x, 80, "E", 8.8, weight="bold", ha="right", color=E_COLOR)
    label(ne_count_x, 80, "NE", 8.8, weight="bold", ha="right", color=NE_COLOR)
    label(fold_x, 63, "Fold", 9.7, weight="bold", ha="right")
    label(fold_x, 80, "E / NE", 8.6, ha="right", color=MUTED)
    max_pct = max(max(row[3:5]) for row, _, _ in entries)
    tick_step = 5 if max_pct <= 30 else 10
    axis_max = max(tick_step, ceil(max_pct / tick_step) * tick_step)
    scale = (plot_end - plot_x) / axis_max
    for tick in range(0, axis_max + 1, tick_step):
        x = plot_x + tick * scale
        label(x, 80, str(tick), 8, ha="center", color=MUTED)
        ax.plot([x, x], [90, baseline], lw=0.5, color="#CDD9DD", zorder=1)
    ax.plot([left, fold_x], [90, 90], lw=0.7, color="#9EAEB5")
    for index, ((code, e, ne, e_pct, ne_pct, fold), name, center) in enumerate(entries):
        if index and index % 5 == 0:
            top = center - (12.25 if "\n" in name else 8.2)
            ax.plot([left, fold_x], [top, top], lw=0.45, color="#DFE6E8", zorder=0)
        label(left, center, code, 9.2, weight="bold", color=MUTED)
        label(label_x, center, name, 9.2)
        e_x, ne_x = plot_x + e_pct * scale, plot_x + ne_pct * scale
        # The vertical separation keeps both marks visible even at zero.
        ax.plot([e_x, ne_x], [center - 2.3, center + 2.3], lw=0.7, color="#9AA9AF", zorder=2)
        ax.plot(e_x, center - 2.3, marker="o", ms=4.9, mec=E_COLOR, mfc=E_COLOR, ls="none", zorder=3)
        ax.plot(ne_x, center + 2.3, marker="s", ms=4.2, mec=NE_COLOR, mfc=NE_COLOR, ls="none", zorder=3)
        label(e_count_x, center, str(e), 9.1, ha="right")
        label(ne_count_x, center, f"{ne:,}", 9.1, ha="right")
        label(fold_x, center, f"{fold:.2f}×" if ne else "∞", 9.3, ha="right", weight="bold")
    ax.plot([left, fold_x], [baseline + 2, baseline + 2], lw=0.7, color="#9EAEB5")
    label(left, baseline + 14, "Fold enrichment is the essential percentage divided by the nonessential percentage. "
          "1× indicates equal prevalence.", 8.0, color=MUTED)
    label(left, baseline + 31, f"Each gene counts once per assigned category. {multi['essential']} essential and "
          f"{multi['nonessential']} nonessential genes have multiple categories.", 8.0, color=MUTED)
    missing_e, missing_ne = (sum(statuses[g].values()) - statuses[g][ASSIGNED] for g in GROUPS)
    label(left, baseline + 48, f"The eggNOG-mapper annotation provides no COG assignment for {missing_e} essential "
          f"and {missing_ne} nonessential genes.", 8.0, color=MUTED)
    label(left, baseline + 65, f"All {len(tested)} categories were tested. Two are listed here with their E and NE gene counts.", 8.0, color=MUTED)
    legend_counts = [f"{code}, {NAMES[code]} ({counts['essential'][code]} E, "
                     f"{counts['nonessential'][code]} NE)" for code in ("A", "B")]
    label(left, baseline + 79, "  ·  ".join(legend_counts), 8.0, color=MUTED)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", type=Path, default=Path(__file__).with_name("figure_4a_cog_patterns.tsv"))
    parser.add_argument("--cohort", type=Path, help="Essentiality workbook with the Table S1 worksheet")
    parser.add_argument("--eggnog", type=Path, help="Original eggNOG-mapper export workbook")
    parser.add_argument("--save-data", type=Path, help="Export grouped annotation pattern counts from the workbooks")
    parser.add_argument("--save-statistics", type=Path, help="Export the full enrichment results and BH-adjusted tests")
    parser.add_argument("--output", type=Path, default=Path("Figure_4A.png"))
    parser.add_argument("--pdf", type=Path, help="PDF destination (default: same stem as PNG)")
    args = parser.parse_args()
    if bool(args.cohort) != bool(args.eggnog) or (args.save_data and not args.eggnog):
        parser.error("--cohort and --eggnog are required together when rebuilding assignments")
    if args.eggnog:
        rows = from_workbooks(args.cohort, args.eggnog)
    else:
        with args.data.open(encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
    counts, statuses, multi = summarize(rows)
    if args.save_statistics:
        write_statistics(args.save_statistics, counts, statuses)
    if args.save_data:
        with args.save_data.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    fig = draw(counts, statuses, multi)
    fig.savefig(args.output, dpi=600, facecolor="white")
    with matplotlib.rc_context({"pdf.fonttype": 42}):
        fig.savefig(args.pdf or args.output.with_suffix(".pdf"), facecolor="white",
                    metadata={"Title": "COG category enrichment in essential genes of Pc9446",
                              "CreationDate": None, "ModDate": None})
    plt.close(fig)
    for group in GROUPS:
        print(f"{group}: {sum(counts[group].values())} gene-category memberships, "
              f"{multi[group]} genes with multiple COG categories. {dict(statuses[group])}")


if __name__ == "__main__":
    main()
