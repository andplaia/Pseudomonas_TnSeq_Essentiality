"""Draw the three-level SEED network for the 479 essential genes in Pc9446.

Run with Python 3.10+ and Matplotlib:
    python plot_figure_3.py

Both the high-resolution PNG and the searchable vector PDF use the same portrait layout.
Each terminal path in figure_3_seed_subsystems.tsv identifies a subsystem,
its class and its superclass in the BV-BRC functional annotation. Counts
are summed along the paths. Each node's
circle area is proportional to its count, with one scale for all levels.
The group layout reserves space for wrapped names before positioning connected
nodes. A repeated name at different levels or under a different parent
identifies a separate node.

To rebuild the input table from the original Cytoscape exports:
    python plot_figure_3.py --nodes "cytoscape_nodes(1).tsv" \\
        --edges "cytoscape_edges(1).tsv" --save-data figure_3_seed_subsystems.tsv
"""

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from math import cos, pi, sqrt
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Circle, Rectangle
from matplotlib.textpath import TextToPath


FIELDS = ("superclass", "class", "subsystem", "essential_genes")
INK = "#192A35"
MUTED = "#405460"
NODE_COLORS = ("#AE017E", "#66C2A4", "#FFCC00")
NODE_OUTLINES = ("#76115D", "#307360", "#906A13")
EDGE_COLORS = ("#A96579", "#CB8648")
BACKGROUND = {
    "METABOLISM": "#F4F1E5", "PROTEIN PROCESSING": "#EAF0F8",
    "CELLULAR PROCESSES": "#F6EBED", "ENERGY": "#EBF2E6",
    "STRESS RESPONSE, DEFENSE, VIRULENCE": "#F7ECE3", "RNA PROCESSING": "#EEEAF6",
    "DNA PROCESSING": "#E6F2F1", "MEMBRANE TRANSPORT": "#E7EEF1",
    "CELL ENVELOPE": "#F3EAF3", "MISCELLANEOUS": "#ECEDE2",
    "REGULATION AND CELL SIGNALING": "#EAF1EC", "Unknown": "#EEEEF0",
}
FONT = "DejaVu Sans"
LEAF_SIZE = 9.6
CLASS_SIZE = 10.2
ROOT_SIZE = 12.6
RADIUS_UNIT = 1.35
COLUMN_WIDTH = 385
REGION_COLUMNS = 3
GAP = 5
MEASURE = TextToPath()

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



def hierarchy(rows):
    counts = Counter()
    children = defaultdict(list)
    for row in rows:
        path = tuple(row[key] for key in FIELDS[:3])
        n = int(row["essential_genes"])
        if not all(path) or path in counts or n < 1:
            raise ValueError(f"Invalid or duplicate SEED path: {path}")
        for level in range(1, 4):
            child = path[:level]
            if child not in counts:
                children[child[:-1]].append(child)
            counts[child] += n
    if not counts:
        raise ValueError("The hierarchy has no assignments.")
    return counts, children


def radius(n):
    return RADIUS_UNIT * sqrt(n)


@lru_cache(maxsize=None)
def text_width(value, size, weight):
    font = FontProperties(family=FONT, size=size, weight=weight)
    return MEASURE.get_text_width_height_descent(value, font, False)[0]


def wrap(value, limit, size, weight="normal"):
    lines = []
    for word in value.split():
        if lines and text_width(lines[-1] + " " + word, size, weight) <= limit:
            lines[-1] += " " + word
        else:
            lines.append(word)
    return "\n".join(lines)


def text_height(value, size):
    return size * (1.28 * (value.count("\n") + 1))


@dataclass
class Leaf:
    path: tuple
    text: str
    height: float


@dataclass
class Fan:
    path: tuple
    text: str
    wings: list
    height: float
    half_lane: float
    class_x: float
    radius: float
    compact: bool


def plan_fan(path, leaves, counts, width):
    half_lane = 45
    class_text = wrap(path[-1], 2 * half_lane, CLASS_SIZE, "bold")
    class_r = radius(counts[path])
    max_leaf_r = max(radius(counts[p]) for p in leaves)
    longest_word = max(text_width(word, LEAF_SIZE, "normal") for leaf in leaves
                       for word in leaf[-1].split())
    two_wing_space = width / 2 - half_lane - 30 - max_leaf_r - 5
    # A narrow compartment must not cut off a systematic chemical name.
    # Use one wing for that class to keep even its longest word intact.
    compact = len(leaves) <= 3 or longest_word > two_wing_space
    class_x = 57 if compact else width / 2
    outer = width - class_x if compact else width / 2
    text_space = outer - half_lane - 30 - max_leaf_r - 5
    # An alphabetically ordered wing makes related names easy to scan. Find
    # the split minimizing the difference in the two wings' measured heights.
    items = []
    for child in sorted(leaves, key=lambda p: p[-1].casefold()):
        words = wrap(child[-1], text_space, LEAF_SIZE)
        items.append(Leaf(child, words, max(text_height(words, LEAF_SIZE),
                                            2 * radius(counts[child])) + 4.1))
    if compact:
        wings = [(1, items)]
    else:
        total = sum(item.height for item in items)
        split = min(range(1, len(items)), key=lambda i: abs(
            total - 2 * sum(item.height for item in items[:i])))
        wings = [(-1, items[:split]), (1, items[split:])]
    wing_height = max(sum(item.height for item in wing) for _, wing in wings)
    center_height = 2 * (class_r + 3 + text_height(class_text, CLASS_SIZE) + 2)
    return Fan(path, class_text, wings, max(wing_height, center_height) + 8,
               half_lane, class_x, class_r, compact)


@dataclass
class Card:
    path: tuple
    columns: list
    width: float
    height: float


def plan_card(root, counts, children):
    ncols = 2 if len(children[root]) >= 6 else 1
    width = COLUMN_WIDTH * ncols + GAP * (ncols - 1)
    col_width = (width - 12) / ncols
    fans = [plan_fan(p, children[p], counts, col_width - 10) for p in children[root]]
    cols = [[] for _ in range(ncols)]
    heights = [0.0] * ncols
    for fan in sorted(fans, key=lambda f: (-f.height, f.path)):
        i = heights.index(min(heights))
        cols[i].append(fan)
        heights[i] += fan.height + GAP
    return Card(root, cols, width, 44 + max(heights) + 5)


def pack_cards(cards):
    # Put the broadest trees down first; the shorter trees fill the spaces
    # beneath them. All columns share a baseline to form the outer rectangle.
    heights = [0.0] * REGION_COLUMNS
    placed = []
    for card in sorted(cards, key=lambda c: (-len(c.columns), -c.height, c.path)):
        span = len(card.columns)
        i = min(range(len(heights) - span + 1), key=lambda j: max(heights[j:j + span]))
        y = max(heights[i:i + span])
        placed.append((i, y, card))
        heights[i:i + span] = [y + card.height + GAP] * span
    max_y = max(heights) - GAP
    output = []
    for i, y, card in placed:
        extend = min(max_y - y - card.height,
                     min((yy - y - card.height - GAP for j, yy, other in placed
                          if yy > y and j < i + len(card.columns)
                          and i < j + len(other.columns)), default=float("inf")))
        output.append((i * (COLUMN_WIDTH + GAP), y, card, card.height + max(0, extend)))
    return output, REGION_COLUMNS * COLUMN_WIDTH + (REGION_COLUMNS - 1) * GAP, max_y


def root_name(value):
    return value.capitalize().replace("Dna", "DNA").replace("Rna", "RNA")


def build_network(rows):
    counts, children = hierarchy(rows)
    cards = [plan_card(p, counts, children) for p in children[()]]
    places, network_width, network_height = pack_cards(cards)
    margin, top = 9, 56
    width, height = network_width + 2 * margin, network_height + top + 25
    fig = plt.figure(figsize=(width / 72, height / 72), dpi=100, facecolor="white")
    ax = fig.add_axes((0, 0, 1, 1), xlim=(0, width), ylim=(height, 0), aspect="equal")
    ax.set_axis_off()
    ids = {p: i for i, p in enumerate(counts)}
    ax.node_paths = {i: p for p, i in ids.items()}

    def text(x, y, s, size, **style):
        return ax.text(x, y, s, family=FONT, fontsize=size, color=style.pop("color", INK),
                       va=style.pop("va", "center"), linespacing=1.09, zorder=4, **style)

    def node(path, x, y, label_x, label_y, label, size, **style):
        depth, i = len(path) - 1, ids[path]
        dot = Circle((x, y), radius(counts[path]), facecolor=NODE_COLORS[depth],
                     edgecolor=NODE_OUTLINES[depth], linewidth=0.36, zorder=3, gid=f"node-{i}")
        ax.add_patch(dot)
        text(label_x, label_y, label, size, gid=f"label-{i}", **style)

    def edge(parent, child, points):
        xs, ys = zip(*points)
        ax.plot(xs, ys, lw=0.67 if len(parent) == 1 else 0.54,
                color=EDGE_COLORS[len(parent) - 1], zorder=2, solid_capstyle="round",
                solid_joinstyle="round", gid=f"edge-{ids[parent]}-{ids[child]}")

    for offset_x, offset_y, card, card_height in places:
        x, y = margin + offset_x, top + offset_y
        root = card.path
        ax.add_patch(Rectangle((x, y), card.width, card_height,
                               facecolor=BACKGROUND.get(root[0], "#ECF0F0"), edgecolor="none",
                               zorder=0, gid=f"region-{ids[root]}"))
        root_r = radius(counts[root])
        root_x, root_y = x + 20, y + 20
        node(root, root_x, root_y, root_x + root_r + 5, root_y, root_name(root[0]),
             ROOT_SIZE, weight="bold")
        col_width = (card.width - 12) / len(card.columns)
        for c, fans in enumerate(card.columns):
            left = x + 6 + c * col_width
            # Extra height is distributed between and around fans within the region.
            slack = max(0, card_height - 49 - sum(f.height + GAP for f in fans))
            pad = slack / (len(fans) + 1)
            cursor = y + 41 + pad
            for fan in fans:
                hub_x, hub_y = left + 8 + fan.class_x, cursor + fan.height / 2
                origin_y = y + 37
                # Enter each class from above through its central lane. The path
                # follows the margin of its column, outside all subsystem labels.
                gutter = left + 1.5
                edge(root, fan.path, [(root_x, root_y), (root_x, origin_y),
                     (gutter, origin_y), (gutter, cursor + 1), (hub_x, cursor + 1),
                     (hub_x, hub_y)])
                node(fan.path, hub_x, hub_y, hub_x, hub_y + fan.radius + 3,
                     fan.text, CLASS_SIZE, weight="bold", ha="center", va="top")
                label_half = max(text_width(line, CLASS_SIZE, "bold")
                                 for line in fan.text.splitlines()) / 2
                for sign, wing in fan.wings:
                    h = sum(item.height for item in wing)
                    row_top = hub_y - h / 2
                    for item in wing:
                        leaf_y = row_top + item.height / 2
                        row_top += item.height
                        position = (leaf_y - hub_y) / max(1, h / 2)
                        arc = 12 + 11 * cos(position * pi / 2)
                        leaf_x = hub_x + sign * (fan.half_lane + arc)
                        leaf_r = radius(counts[item.path])
                        # Direct spokes connect nearby nodes. Route a spoke around the
                        # class name only when it would otherwise pass through its text.
                        points = [(hub_x, hub_y)]
                        dy = leaf_y - hub_y
                        if dy > fan.radius + 2:
                            dx_at_name = abs(leaf_x - hub_x) * (fan.radius + 2) / dy
                            if dx_at_name < label_half + 2:
                                points.append((hub_x + sign * (label_half + 3),
                                               hub_y + fan.radius + 1))
                        edge(fan.path, item.path, points + [(leaf_x, leaf_y)])
                        node(item.path, leaf_x, leaf_y, leaf_x + sign * (leaf_r + 2.4),
                             leaf_y, item.text, LEAF_SIZE, ha="left" if sign > 0 else "right")
                cursor += fan.height + GAP + pad

    total_genes = sum(counts[root] for root in children[()])
    levels = Counter(map(len, counts))
    text(margin, 15, "SEED functional hierarchy of essential genes", 20, weight="bold")
    text(margin, 38, f"{total_genes} essential genes   ·   {levels[1]} superclasses   ·   "
         f"{levels[2]} classes   ·   {levels[3]} subsystem entries", 11.4, color=MUTED)
    # Use the same size function in the key as in the network.
    legend_x = width - 600
    for name, face, stroke in zip(("Superclass", "Class", "Subsystem"), NODE_COLORS, NODE_OUTLINES):
        ax.add_patch(Circle((legend_x, 25), 4, facecolor=face, edgecolor=stroke,
                            linewidth=0.45, zorder=3))
        text(legend_x + 8, 25, name, 11.3)
        legend_x += 87
    text(legend_x + 3, 25, "Area ∝ gene count", 11.3, color=MUTED)
    legend_x += 143
    for n in (1, 16, 128):
        r = radius(n)
        ax.add_patch(Circle((legend_x + r, 25), r, facecolor="#D6DFE4",
                            edgecolor="#52616A", linewidth=0.45, zorder=3))
        text(legend_x + 2 * r + 5, 25, str(n), 11.2)
        legend_x += 2 * r + 38
    text(margin, height - 12, "Unknown denotes genes without a SEED subsystem assignment "
         "in the BV-BRC functional annotation.", 10.2, color=MUTED)
    return fig, ax


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", type=Path, default=Path(__file__).with_name("figure_3_seed_subsystems.tsv"))
    parser.add_argument("--nodes", type=Path)
    parser.add_argument("--edges", type=Path)
    parser.add_argument("--save-data", type=Path)
    parser.add_argument("--output", type=Path, default=Path("Figure_3.png"))
    parser.add_argument("--pdf", type=Path, help="Vector PDF destination (default: same stem as PNG)")
    parser.add_argument("--dpi", type=int, default=450, help="PNG resolution (default: 450 dpi)")
    args = parser.parse_args()
    if bool(args.nodes) != bool(args.edges) or (args.save_data and not args.nodes):
        parser.error("--nodes and --edges are both required when rebuilding the hierarchy")
    rows = paths_from_exports(args.nodes, args.edges) if args.nodes else read_table(args.data)
    if args.save_data:
        with args.save_data.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    fig, _ = build_network(rows)
    fig.savefig(args.output, dpi=args.dpi, facecolor="white")
    # Type 42 embeds the TrueType face and retains searchable Unicode text.
    with matplotlib.rc_context({"pdf.fonttype": 42}):
        fig.savefig(args.pdf or args.output.with_suffix(".pdf"), facecolor="white",
                    metadata={"Title": "SEED functional hierarchy of essential genes in Pc9446",
                              "CreationDate": None, "ModDate": None})
    plt.close(fig)


if __name__ == "__main__":
    main()
