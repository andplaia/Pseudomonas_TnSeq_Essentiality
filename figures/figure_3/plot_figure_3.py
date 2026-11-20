"""Draw Figure 3 from the SEED assignments for the 479 Pc9446 essential genes.

Run with Python 3.10+ and Matplotlib:
    python plot_figure_3.py

The companion TSV lists every subsystem, its class and superclass, and its
essential-gene count. It preserves the 131 terminal paths in the submitted
Cytoscape node and edge exports. The plot sums counts up each path and shows
the two upper levels in decreasing order of superclass size. Equal counts
are ordered alphabetically. Unknown is an explicit category in the source.

To recreate the companion TSV and plot from the original exports:
    python plot_figure_3.py --nodes "cytoscape_nodes(1).tsv" \
        --edges "cytoscape_edges(1).tsv" --save-data figure_3_seed_subsystems.tsv
"""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.textpath import TextToPath


FIELDS = ("superclass", "class", "subsystem", "essential_genes")
INK = "#202C31"
MUTED = "#485B65"
RULE = "#CBD4D8"
BAR = "#4682B4"


def read_table(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def paths_from_exports(nodes_path, edges_path):
    node_rows = read_table(nodes_path)
    nodes = {r["id"]: r for r in node_rows}
    if len(nodes) != len(node_rows):
        raise ValueError("Duplicate SEED node identifier.")
    parents, child_totals = {}, Counter()
    for edge in read_table(edges_path):
        parent, child = edge["source"], edge["target"]
        relationship = (nodes[parent]["level"], nodes[child]["level"])
        if (edge["interaction"] != "is_a" or relationship not in {("SC", "CL"), ("CL", "SS")}
                or child in parents or int(edge["weight"]) != int(nodes[child]["count"])):
            raise ValueError(f"Invalid SEED hierarchy relationship: {parent} -> {child}")
        parents[child] = parent
        child_totals[parent] += int(edge["weight"])
    rows = []
    for key, node in nodes.items():
        if (node["level"] != "SC") != (key in parents):
            raise ValueError(f"Invalid parent for node {key}")
        if node["level"] != "SS":
            if node["level"] not in {"SC", "CL"} or child_totals[key] != int(node["count"]):
                raise ValueError(f"Child counts do not match the total for {key}")
            continue
        class_id = parents[key]
        path = (nodes[parents[class_id]]["label"], nodes[class_id]["label"], node["label"])
        if path[:2] != (node["Superclass"], node["Class"]):
            raise ValueError(f"Conflicting SEED node and edge assignments for {key}")
        rows.append(dict(zip(FIELDS, (*path, int(node["count"])))))
    return rows


def class_counts(rows):
    groups = defaultdict(Counter)
    seen = set()
    for row in rows:
        path = tuple(row[field] for field in FIELDS[:3])
        count = int(row["essential_genes"])
        if not all(path) or path in seen or count < 1:
            raise ValueError(f"Invalid or duplicate SEED path: {path}")
        seen.add(path)
        groups[path[0]][path[1]] += count
    if not groups:
        raise ValueError("The hierarchy contains no subsystem assignments.")
    return groups


def wrap_label(label, width, size, weight="normal"):
    font = FontProperties(family="DejaVu Sans", size=size, weight=weight)
    measure = TextToPath()
    lines = []
    for word in label.split():
        if lines:
            extended = f"{lines[-1]} {word}"
            if measure.get_text_width_height_descent(extended, font, False)[0] <= width:
                lines[-1] = extended
                continue
        lines.append(word)
    return "\n".join(lines)


def superclass_label(label):
    # Only adjust typesetting. The companion retains the exported capitalization.
    return label.capitalize().replace("Dna", "DNA").replace("Rna", "RNA")


def render(groups, output):
    width = 534  # all coordinates and font sizes are in points
    parent_x, total_x, branch_x, class_x = 12, 133, 144, 151
    bar_x, bar_end, count_x = 363, 491, 521
    total_genes = sum(sum(counts.values()) for counts in groups.values())
    nclasses = sum(len(counts) for counts in groups.values())
    max_class = max(max(counts.values()) for counts in groups.values())
    ordered = sorted(groups.items(), key=lambda item: (-sum(item[1].values()), item[0].casefold()))

    top = 70
    layout = []
    for parent, counts in ordered:
        first_y = top
        parent_text = wrap_label(superclass_label(parent), total_x - parent_x - 21, 9, "bold")
        entries = []
        for name, count in sorted(counts.items(), key=lambda item: (-item[1], item[0].casefold())):
            label = wrap_label(name, bar_x - class_x - 8, 8.6)
            row_height = max(14, 10 * len(label.splitlines()) + 4)
            if len(counts) == 1:
                row_height = max(row_height, 10 * len(parent_text.splitlines()) + 4)
            entries.append((label, count, top + row_height / 2))
            top += row_height
        layout.append((parent_text, sum(counts.values()), (first_y + top) / 2, entries, top + 3))
        top += 6

    height = top + 19
    fig = plt.figure(figsize=(width / 72, height / 72), dpi=120, facecolor="white")
    ax = fig.add_axes((0, 0, 1, 1), xlim=(0, width), ylim=(height, 0))
    ax.set_axis_off()

    def label(x, y, value, size, **kwargs):
        ax.text(x, y, value, fontsize=size, va="center", color=kwargs.pop("color", INK),
                family="DejaVu Sans", **kwargs)

    label(parent_x, 14, "SEED functional assignments of essential genes", 14.4, weight="bold")
    label(parent_x, 31, f"{total_genes} essential genes   ·   {len(groups)} superclasses   ·   "
          f"{nclasses} classes", 9.5, color=MUTED)
    label(parent_x, 53, "Superclass", 9.4, weight="bold")
    label(total_x, 53, "Genes", 8.7, ha="right", weight="bold")
    label(class_x, 53, "Class", 9.4, weight="bold")
    label(bar_x, 48, "Genes per class", 9.1, weight="bold")

    scale = (bar_end - bar_x) / max_class
    # Direct counts appear alongside the bars. Ticks provide a shared scale.
    tick_step = 50 if max_class >= 100 else max(1, round(max_class / 3))
    for tick in range(0, max_class + 1, tick_step):
        x = bar_x + tick * scale
        label(x, 61, str(tick), 7.6, ha="center", color=MUTED)
        ax.plot([x, x], [67, layout[-1][-1]], lw=0.4, color="#E6EBED", zorder=0)
    ax.plot([parent_x, count_x], [67, 67], lw=0.75, color="#9DAEB5")
    for parent, total, middle, entries, divider_y in layout:
        label(parent_x, middle, parent, 9, weight="bold", linespacing=1.05)
        label(total_x, middle, str(total), 9, ha="right", weight="bold")
        ax.plot([branch_x, branch_x], [entries[0][2], entries[-1][2]], lw=0.7, color="#A1B2BA")
        for text, count, center in entries:
            ax.plot([branch_x, class_x - 3.5], [center, center], lw=0.7, color="#A1B2BA")
            label(class_x, center, text, 8.6, linespacing=1.05)
            ax.add_patch(Rectangle((bar_x, center - 3.2), count * scale, 6.4,
                                   facecolor=BAR, edgecolor="none"))
            label(count_x, center, str(count), 8.8, ha="right")
        ax.plot([parent_x, count_x], [divider_y, divider_y], lw=0.42, color=RULE)
    label(parent_x, height - 9, "Unknown denotes genes without a SEED assignment in Table S1.",
          8.4, color=MUTED)
    fig.savefig(output, dpi=600, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", type=Path, default=Path(__file__).with_name("figure_3_seed_subsystems.tsv"))
    parser.add_argument("--nodes", type=Path)
    parser.add_argument("--edges", type=Path)
    parser.add_argument("--save-data", type=Path)
    parser.add_argument("--output", type=Path, default=Path("Figure_3.png"))
    args = parser.parse_args()
    if bool(args.nodes) != bool(args.edges) or (args.save_data and not args.nodes):
        parser.error("--nodes and --edges are both required when rebuilding the hierarchy")
    rows = paths_from_exports(args.nodes, args.edges) if args.nodes else read_table(args.data)
    groups = class_counts(rows)
    if args.save_data:
        with args.save_data.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    render(groups, args.output)


if __name__ == "__main__":
    main()
