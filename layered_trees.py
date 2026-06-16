"""
layered_trees.py
────────────────────────────────────────────────────────────────
Build layered (Sugiyama / hierarchical) layouts for each human
population gene tree and export:

  outputs/
    <GENE>_layered.png        — high-res plot (300 dpi)
    <GENE>_layered_data.json  — full layer data (positions, edges, AFs)
    all_layered_data.json     — combined dataset for all genes
    layer_summary.csv         — tabular summary of all nodes × genes

Layout algorithm
────────────────
1. BFS from root → assign layer index (depth)
2. Reingold-Tilford x-positioning within each layer
   (leaves spaced evenly; internals centred over their children)
3. Branch lengths encoded as VERTICAL spacing between layers
   (longer branches → more vertical distance)

Requirements:
    pip install networkx matplotlib
"""

# /// script
# requires-python = ">=3.10"
# dependencies = ["networkx", "matplotlib"]
# ///

import csv
import json
import math
import os
import pickle
from collections import defaultdict, deque

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
import networkx as nx

# ── Paths ────────────────────────────────────────────────────────
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
NETWORKS_DIR = os.path.join(SCRIPT_DIR, "networks")
OUT_DIR      = os.path.join(SCRIPT_DIR, "outputs")
os.makedirs(OUT_DIR, exist_ok=True)

GENES = ["LCT", "EPAS1", "HLA-A", "GAPDH", "MAPT"]

# ── Colour palette (one colour per population) ───────────────────
POP_COLORS = {
    "African":                "#f87171",   # red-400
    "American":               "#fb923c",   # orange-400
    "Ashkenazi Jewish":       "#facc15",   # yellow-400
    "East Asian":             "#4ade80",   # green-400
    "Finnish":                "#34d399",   # emerald-400
    "Middle Eastern":         "#22d3ee",   # cyan-400
    "Non-Finnish European":   "#818cf8",   # indigo-400
    "South Asian":            "#e879f9",   # fuchsia-400
}
INTERNAL_COLOR = "#334155"   # slate-700
EDGE_COLOR     = "#475569"   # slate-600
BG_COLOR       = "#0b1120"
TEXT_COLOR      = "#e2e8f8"


# ════════════════════════════════════════════════════════════════
#  STEP 1 – Load NetworkX graphs
# ════════════════════════════════════════════════════════════════

def load_trees() -> dict[str, nx.DiGraph]:
    pkl = os.path.join(NETWORKS_DIR, "all_trees.pkl")
    if os.path.exists(pkl):
        with open(pkl, "rb") as f:
            return pickle.load(f)
    # fallback: load per-gene GraphML
    trees = {}
    for gene in GENES:
        safe = gene.replace("-", "_")
        path = os.path.join(NETWORKS_DIR, f"{safe}.graphml")
        if os.path.exists(path):
            trees[gene] = nx.read_graphml(path)
    return trees


# ════════════════════════════════════════════════════════════════
#  STEP 2 – Layered layout computation
# ════════════════════════════════════════════════════════════════

def assign_layers(G: nx.DiGraph) -> dict[str, int]:
    """BFS from root → layer index (depth from root)."""
    root = G.graph.get("root") or next(n for n in G if G.in_degree(n) == 0)
    layers = {root: 0}
    queue  = deque([root])
    while queue:
        node = queue.popleft()
        for child in G.successors(node):
            if child not in layers:
                layers[child] = layers[node] + 1
                queue.append(child)
    return layers


def assign_x_positions(G: nx.DiGraph, layers: dict[str, int]) -> dict[str, float]:
    """
    Reingold–Tilford-inspired x assignment:
      - leaves at each layer are spread evenly
      - internal nodes centred over their children
    Bottom-up traversal (leaves first).
    """
    root = G.graph.get("root") or next(n for n in G if G.in_degree(n) == 0)

    # topological order → reverse for bottom-up
    topo = list(nx.topological_sort(G))
    leaves = [n for n in topo if G.out_degree(n) == 0]

    # assign leaf x-positions left to right
    x = {}
    leaf_counter = [0]

    def place(node):
        children = list(G.successors(node))
        if not children:
            x[node] = float(leaf_counter[0])
            leaf_counter[0] += 1
        else:
            for c in children:
                place(c)
            x[node] = (x[children[0]] + x[children[-1]]) / 2

    place(root)

    # Normalise to [0, 1]
    mn, mx = min(x.values()), max(x.values())
    span = mx - mn if mx != mn else 1.0
    return {n: (v - mn) / span for n, v in x.items()}


def assign_y_positions(
    G: nx.DiGraph,
    layers: dict[str, int],
    use_branch_lengths: bool = False,
) -> dict[str, float]:
    """
    y position based on layer depth (uniform spacing).
    Returned as values in [0, 1] where 0 = root (top), 1 = deepest leaf.
    """
    if not use_branch_lengths:
        max_layer = max(layers.values()) or 1
        return {n: lyr / max_layer for n, lyr in layers.items()}

    # Cumulative branch length from root
    root = G.graph.get("root") or next(n for n in G if G.in_degree(n) == 0)
    depth = {root: 0.0}
    for node in nx.topological_sort(G):
        bl = G.nodes[node].get("branch_length", 0.0)
        if isinstance(bl, str):
            bl = float(bl)
        for child in G.successors(node):
            cbl = G.nodes[child].get("branch_length", 0.0)
            if isinstance(cbl, str):
                cbl = float(cbl)
            depth[child] = depth[node] + cbl

    mn, mx = min(depth.values()), max(depth.values())
    span = mx - mn if mx != mn else 1.0
    return {n: (v - mn) / span for n, v in depth.items()}


def compute_layout(G: nx.DiGraph) -> dict[str, dict]:
    """
    Returns a dict: node_id → {x, y, layer, label, is_leaf, af, color}
    x in [0,1], y in [0,1] (0=root at top, 1=leaves at bottom)
    """
    layers  = assign_layers(G)
    x_pos   = assign_x_positions(G, layers)
    y_pos   = assign_y_positions(G, layers, use_branch_lengths=True)

    layout = {}
    for node, data in G.nodes(data=True):
        label   = str(data.get("label", ""))
        is_leaf = bool(data.get("is_leaf", False))
        af_raw  = data.get("af", -1.0)
        af      = float(af_raw) if af_raw is not None else -1.0

        color = POP_COLORS.get(label, INTERNAL_COLOR) if is_leaf else INTERNAL_COLOR

        layout[node] = {
            "x":       round(x_pos.get(node, 0.0), 6),
            "y":       round(y_pos.get(node, 0.0), 6),
            "layer":   layers.get(node, 0),
            "label":   label,
            "is_leaf": is_leaf,
            "af":      round(af, 7) if af >= 0 else None,
            "color":   color,
        }
    return layout


# ════════════════════════════════════════════════════════════════
#  STEP 3 – Plot
# ════════════════════════════════════════════════════════════════

def plot_layered_tree(
    G: nx.DiGraph,
    layout: dict,
    gene: str,
    description: str,
    ax: plt.Axes | None = None,
    standalone: bool = True,
) -> plt.Figure | None:
    """
    Draw a layered tree on ax.
    If standalone=True, creates its own figure and saves to disk.
    """
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 7), facecolor=BG_COLOR)
        ax.set_facecolor(BG_COLOR)

    ax.set_facecolor(BG_COLOR)

    # ── Draw edges (L-shaped elbows) ────────────────────────────
    for u, v in G.edges():
        xu, yu = layout[u]["x"], layout[u]["y"]
        xv, yv = layout[v]["x"], layout[v]["y"]
        ax.plot([xu, xu], [yu, yv], color=EDGE_COLOR, linewidth=1.6,
                solid_capstyle="round", zorder=1)
        ax.plot([xu, xv], [yv, yv], color=EDGE_COLOR, linewidth=1.6,
                solid_capstyle="round", zorder=1)

    # ── Draw nodes ──────────────────────────────────────────────
    max_layer = max(v["layer"] for v in layout.values()) or 1
    for node_id, info in layout.items():
        x, y = info["x"], info["y"]
        color   = info["color"]
        is_leaf = info["is_leaf"]
        label   = info["label"]

        if is_leaf:
            ax.scatter(x, y, s=140, color=color, zorder=4,
                       edgecolors="#0f172a", linewidths=1.8)
            # Population label — horizontal, below dot
            ax.text(x, y + 0.06, label,
                    ha="center", va="bottom",
                    fontsize=8, color=color,
                    fontweight="bold", zorder=5)
            # AF value
            if info["af"] is not None:
                ax.text(x, y + 0.015, f"{info['af']:.4f}",
                        ha="center", va="bottom",
                        fontsize=6, color="#94a3b8", alpha=0.9,
                        zorder=5)
        else:
            ax.scatter(x, y, s=60, color="#1e3a5f", zorder=3,
                       edgecolors="#475569", linewidths=1.2)

    # ── Layer labels on left axis ────────────────────────────────
    layers_present = sorted(set(v["layer"] for v in layout.values()))
    for lyr in layers_present:
        y_vals = [v["y"] for v in layout.values() if v["layer"] == lyr]
        y_lyr  = sum(y_vals) / len(y_vals)
        if lyr > 0:
            ax.axhline(y_lyr, color="#1e293b", linewidth=0.6,
                       linestyle=":", alpha=0.4, zorder=0)
        ax.text(-0.12, y_lyr, f"L{lyr}",
                ha="right", va="center",
                fontsize=7, color="#475569")

    # ── Axes cosmetics ───────────────────────────────────────────
    ax.set_xlim(-0.18, 1.12)
    ax.set_ylim(-0.08, 1.22)   # root at bottom=0, leaves at top=1
    ax.invert_yaxis()          # root at top visually
    ax.axis("off")

    # Title
    ax.set_title(
        f"{gene}\n{description}",
        color=TEXT_COLOR, fontsize=10.5, fontweight="bold",
        pad=14, loc="left",
    )

    # Population legend (standalone only)
    if standalone:
        present_pops = {info["label"] for info in layout.values() if info["is_leaf"]}
        handles = [
            mpatches.Patch(color=POP_COLORS[lbl], label=lbl)
            for lbl in POP_COLORS if lbl in present_pops
        ]
        leg = ax.legend(
            handles=handles,
            loc="lower right",
            fontsize=7.5,
            framealpha=0.18,
            facecolor="#0f172a",
            edgecolor="#334155",
            labelcolor=TEXT_COLOR,
            title="Population",
            title_fontsize=8,
        )
        leg.get_title().set_color(TEXT_COLOR)

    if standalone and fig is not None:
        out_path = os.path.join(OUT_DIR, f"{gene.replace('-','_')}_layered.png")
        fig.savefig(out_path, dpi=300, bbox_inches="tight",
                    facecolor=BG_COLOR, edgecolor="none")
        plt.close(fig)
        print(f"  → {out_path}")
        return None
    return fig


# ════════════════════════════════════════════════════════════════
#  STEP 4 – Panel plot (all genes side-by-side)
# ════════════════════════════════════════════════════════════════

def plot_panel(all_trees, all_layouts, all_descs):
    n = len(all_trees)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 8), facecolor=BG_COLOR)
    fig.suptitle(
        "Human Population Gene Trees  ·  Layered Hierarchical Layout",
        color=TEXT_COLOR, fontsize=12, fontweight="bold", y=1.01,
    )

    for ax, (gene, G) in zip(axes, all_trees.items()):
        ax.set_facecolor(BG_COLOR)
        plot_layered_tree(
            G, all_layouts[gene], gene, all_descs.get(gene, ""),
            ax=ax, standalone=False,
        )

    # Shared legend
    handles = [mpatches.Patch(color=c, label=lbl) for lbl, c in POP_COLORS.items()]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=4,
        fontsize=7,
        framealpha=0.15,
        facecolor=BG_COLOR,
        edgecolor="#334155",
        labelcolor=TEXT_COLOR,
        bbox_to_anchor=(0.5, -0.04),
    )

    fig.tight_layout(pad=2)
    out_path = os.path.join(OUT_DIR, "all_genes_layered_panel.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight",
                facecolor=BG_COLOR, edgecolor="none")
    plt.close(fig)
    print(f"  → {out_path}  (panel)")


# ════════════════════════════════════════════════════════════════
#  STEP 5 – Export layer data
# ════════════════════════════════════════════════════════════════

def export_json(all_layouts, all_trees, all_descs):
    """Export per-gene and combined layer data as JSON."""
    combined = {}

    for gene, layout in all_layouts.items():
        G    = all_trees[gene]
        desc = all_descs.get(gene, "")

        # Build edge list with branch length
        edges = []
        for u, v, edata in G.edges(data=True):
            edges.append({
                "source":        u,
                "target":        v,
                "branch_length": round(float(edata.get("weight", 0.0)), 7),
            })

        gene_data = {
            "gene":        gene,
            "description": desc,
            "root":        G.graph.get("root"),
            "n_nodes":     G.number_of_nodes(),
            "n_edges":     G.number_of_edges(),
            "nodes":       {nid: info for nid, info in layout.items()},
            "edges":       edges,
            "layers": {
                str(lyr): [
                    nid for nid, info in layout.items()
                    if info["layer"] == lyr
                ]
                for lyr in sorted(set(v["layer"] for v in layout.values()))
            },
        }
        combined[gene] = gene_data

        # Per-gene file
        safe = gene.replace("-", "_")
        path = os.path.join(OUT_DIR, f"{safe}_layered_data.json")
        with open(path, "w") as f:
            json.dump(gene_data, f, indent=2)
        print(f"  → {path}")

    # Combined file
    combined_path = os.path.join(OUT_DIR, "all_layered_data.json")
    with open(combined_path, "w") as f:
        json.dump(combined, f, indent=2)
    print(f"  → {combined_path}  (combined)")

    return combined


def export_csv(all_layouts, all_trees):
    """Export a flat CSV: gene, node_id, layer, x, y, label, is_leaf, af."""
    csv_path = os.path.join(OUT_DIR, "layer_summary.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "gene", "node_id", "layer", "x", "y",
            "label", "is_leaf", "af",
        ])
        for gene, layout in all_layouts.items():
            for nid, info in layout.items():
                writer.writerow([
                    gene, nid,
                    info["layer"], info["x"], info["y"],
                    info["label"], info["is_leaf"],
                    info["af"] if info["af"] is not None else "",
                ])
    print(f"  → {csv_path}  (CSV summary)")


# ════════════════════════════════════════════════════════════════
#  MAIN
# ════════════════════════════════════════════════════════════════

def main():
    print("Loading NetworkX trees …")
    trees = load_trees()
    if not trees:
        print("ERROR: no tree files found in networks/. Run pangenome_nx.py first.")
        return

    all_layouts = {}
    all_descs   = {}

    print(f"\nComputing layered layouts for {len(trees)} genes …")
    for gene, G in trees.items():
        print(f"\n── {gene} ─────────────────────────────────────")
        layout = compute_layout(G)
        all_layouts[gene] = layout
        all_descs[gene]   = G.graph.get("description", "")

        # Layer summary to stdout
        by_layer = defaultdict(list)
        for nid, info in layout.items():
            by_layer[info["layer"]].append(info)
        for lyr in sorted(by_layer):
            nodes_in_layer = by_layer[lyr]
            names = [n["label"] or n.get("node_id", "?") for n in nodes_in_layer]
            print(f"  Layer {lyr}: {len(nodes_in_layer)} node(s)  →  {names}")

    print("\n── Plotting individual trees ────────────────────────────")
    for gene, G in trees.items():
        plot_layered_tree(
            G, all_layouts[gene], gene, all_descs.get(gene, ""),
            standalone=True,
        )

    print("\n── Plotting panel ──────────────────────────────────────")
    plot_panel(trees, all_layouts, all_descs)

    print("\n── Exporting JSON data ─────────────────────────────────")
    export_json(all_layouts, trees, all_descs)

    print("\n── Exporting CSV summary ───────────────────────────────")
    export_csv(all_layouts, trees)

    # ── Print a human-readable data table ───────────────────────
    print("\n══════ LAYER DATA PREVIEW ══════════════════════════════")
    for gene, layout in all_layouts.items():
        print(f"\n{gene}  ({all_descs.get(gene,'')})")
        print(f"  {'Node':30s}  {'Layer':5s}  {'x':6s}  {'y':6s}  {'AF':8s}")
        print(f"  {'-'*30}  {'-'*5}  {'-'*6}  {'-'*6}  {'-'*8}")
        for nid, info in sorted(layout.items(), key=lambda kv: (kv[1]['layer'], kv[1]['x'])):
            name = info['label'] if info['label'] else f"[{nid}]"
            af_s = f"{info['af']:.5f}" if info['af'] is not None else "   —"
            print(f"  {name:30s}  {info['layer']:5d}  {info['x']:6.3f}  {info['y']:6.3f}  {af_s:8s}")

    print(f"\n✓ All outputs written to {OUT_DIR}/")


if __name__ == "__main__":
    main()
