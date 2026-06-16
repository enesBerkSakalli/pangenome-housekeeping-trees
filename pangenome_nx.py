"""
pangenome_nx.py
────────────────────────────────────────────────────────────────
Build human population gene trees as NetworkX DiGraph objects.

For each gene (LCT, EPAS1, HLA-A, GAPDH, MAPT) we:
  1. Load gnomAD population allele-frequency data
  2. Compute a pairwise AF-distance matrix across 8 populations
  3. Run UPGMA clustering → Newick string
  4. Parse the Newick back into a NetworkX DiGraph

Each node carries:
  - label        : human-readable population name (leaf) or clade id (internal)
  - is_leaf      : bool
  - branch_length: edge weight to parent
  - af           : allele frequency in gnomAD (leaves only, else None)

Exports (per gene, into ./networks/):
  - <GENE>.graphml   – for Gephi / Cytoscape / networkx.read_graphml()
  - <GENE>.json      – node-link JSON  (networkx.node_link_data)
  - <GENE>.pkl       – Python pickle   (fastest for re-use in notebooks)

A summary pickle `all_trees.pkl` is also written:
  { gene: nx.DiGraph, ... }

Usage:
    python pangenome_nx.py

Requirements:
    pip install networkx
"""

# /// script
# requires-python = ">=3.10"
# dependencies = ["networkx"]
# ///

import copy
import json
import math
import os
import pickle
import re

import networkx as nx
from networkx.readwrite import json_graph

# ── Paths ────────────────────────────────────────────────────────
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
SCRATCH_DIR = os.path.join(SCRIPT_DIR, "scratch")
OUT_DIR     = os.path.join(SCRIPT_DIR, "networks")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Population mapping ───────────────────────────────────────────
POP_ORDER = ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas"]
POP_LABELS = {
    "afr": "African",
    "amr": "American",
    "asj": "Ashkenazi Jewish",
    "eas": "East Asian",
    "fin": "Finnish",
    "mid": "Middle Eastern",
    "nfe": "Non-Finnish European",
    "sas": "South Asian",
}

GENES = {
    "LCT":   {"file": "LCT_variants.json",   "description": "Lactase (dairy tolerance adaptation)"},
    "EPAS1": {"file": "EPAS1_variants.json",  "description": "HIF2α (high-altitude adaptation)"},
    "HLA-A": {"file": "HLA-A_variants.json",  "description": "HLA Class I (immune recognition)"},
    "GAPDH": {"file": "GAPDH_variants.json",  "description": "GAPDH (housekeeping control)"},
    "MAPT":  {"file": "MAPT_variants.json",   "description": "Microtubule-associated protein tau"},
}


# ════════════════════════════════════════════════════════════════
#  STEP 1 – Load gnomAD data
# ════════════════════════════════════════════════════════════════

def load_variants(filepath: str) -> list:
    with open(filepath) as f:
        data = json.load(f)
    try:
        return data["data"]["gene"]["variants"]
    except (KeyError, TypeError):
        try:
            return data["data"]["region"]["variants"]
        except (KeyError, TypeError):
            return []


def compute_pop_afs(variants: list) -> dict[str, float]:
    """Return {pop_id: allele_frequency} for POP_ORDER populations."""
    ac = {p: 0 for p in POP_ORDER}
    an = {p: 0 for p in POP_ORDER}
    for v in variants:
        for src_key in ("exome", "genome"):
            src = v.get(src_key)
            if not src:
                continue
            for pop in src.get("populations", []):
                pid = pop.get("id")
                if pid in POP_ORDER:
                    ac[pid] += pop.get("ac", 0)
                    an[pid] += pop.get("an", 0)
    return {p: (ac[p] / an[p] if an[p] > 0 else 0.0) for p in POP_ORDER}


# ════════════════════════════════════════════════════════════════
#  STEP 2 – Distance matrix & UPGMA
# ════════════════════════════════════════════════════════════════

def distance_matrix(pop_afs: dict[str, float]) -> list[list[float]]:
    n = len(POP_ORDER)
    mat = [[0.0] * n for _ in range(n)]
    for i, p1 in enumerate(POP_ORDER):
        for j, p2 in enumerate(POP_ORDER):
            if i < j:
                d = abs(pop_afs[p1] - pop_afs[p2])
                mat[i][j] = mat[j][i] = d
    return mat


def upgma(matrix: list[list[float]], labels: list[str]) -> str:
    """Run UPGMA; returns a Newick string."""
    dists   = copy.deepcopy(matrix)
    # (newick_label, cluster_size)
    clusters = [(lbl, 1) for lbl in labels]
    heights  = [0.0] * len(labels)

    while len(clusters) > 1:
        # find closest pair
        best_d, bi, bj = math.inf, 0, 1
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                if dists[i][j] < best_d:
                    best_d, bi, bj = dists[i][j], i, j

        new_h   = best_d / 2.0
        bl_i    = new_h - heights[bi]
        bl_j    = new_h - heights[bj]
        ni, si  = clusters[bi]
        nj, sj  = clusters[bj]
        new_nwk = f"({ni}:{bl_i:.6f},{nj}:{bl_j:.6f})"
        new_sz  = si + sj

        # compute distances to new cluster
        keep = [k for k in range(len(clusters)) if k not in (bi, bj)]
        new_dists_to_others = [
            (dists[bi][k] * si + dists[bj][k] * sj) / new_sz
            for k in keep
        ]

        # rebuild
        new_clusters = [clusters[k] for k in keep] + [(new_nwk, new_sz)]
        new_heights  = [heights[k]  for k in keep] + [new_h]
        m = len(new_clusters)
        new_mat = [[0.0] * m for _ in range(m)]
        for a in range(len(keep)):
            for b in range(len(keep)):
                new_mat[a][b] = dists[keep[a]][keep[b]]
            new_mat[a][m - 1] = new_mat[m - 1][a] = new_dists_to_others[a]

        clusters, heights, dists = new_clusters, new_heights, new_mat

    return clusters[0][0] + ";"


# ════════════════════════════════════════════════════════════════
#  STEP 3 – Newick → NetworkX DiGraph
# ════════════════════════════════════════════════════════════════

def newick_to_nx(newick: str, pop_afs: dict[str, float]) -> nx.DiGraph:
    """
    Parse Newick string and return a NetworkX DiGraph.

    Nodes
    -----
    Each node has:
      node_id       : str  — unique id ("root", "internal_0", …, or pop name)
      label         : str  — display name
      is_leaf       : bool
      af            : float | None  — allele frequency (leaves only)
      branch_length : float         — length of edge *from parent*

    Edges
    -----
    parent → child, weight = branch_length
    """
    G  = nx.DiGraph()
    _counter = {"n": 0}

    def _uid(prefix="internal"):
        _counter["n"] += 1
        return f"{prefix}_{_counter['n']}"

    # Tokenise
    tokens = re.findall(r"[^:;,()\s]+|[:;,()]", newick)
    pos    = [0]

    def parse():
        node_id = None

        if pos[0] < len(tokens) and tokens[pos[0]] == "(":
            pos[0] += 1                           # consume '('
            node_id = _uid()
            G.add_node(node_id, label="", is_leaf=False, af=-1.0, branch_length=0.0)

            # first child
            child_id = parse()
            G.add_edge(node_id, child_id, weight=G.nodes[child_id]["branch_length"])

            while pos[0] < len(tokens) and tokens[pos[0]] == ",":
                pos[0] += 1                       # consume ','
                child_id = parse()
                G.add_edge(node_id, child_id, weight=G.nodes[child_id]["branch_length"])

            pos[0] += 1                           # consume ')'

        # optional name
        if pos[0] < len(tokens) and tokens[pos[0]] not in (":", ",", ")", ";"):
            raw_name = tokens[pos[0]]; pos[0] += 1
            if node_id is None:
                # leaf: use population name directly as node_id if unique
                leaf_id = raw_name.replace(" ", "_")
                node_id = leaf_id
                G.add_node(
                    node_id,
                    label=raw_name,
                    is_leaf=True,
                    af=float(pop_afs.get(raw_name, -1.0)),
                    branch_length=0.0,
                )
            else:
                G.nodes[node_id]["label"] = raw_name
        elif node_id is None:
            node_id = _uid("leaf")
            G.add_node(node_id, label="", is_leaf=True, af=-1.0, branch_length=0.0)

        # optional branch length
        if pos[0] < len(tokens) and tokens[pos[0]] == ":":
            pos[0] += 1
            bl = float(tokens[pos[0]]); pos[0] += 1
            G.nodes[node_id]["branch_length"] = bl

        return node_id

    root_id = parse()
    G.graph["root"] = root_id
    return G


# ════════════════════════════════════════════════════════════════
#  STEP 4 – Build & export all trees
# ════════════════════════════════════════════════════════════════

def build_and_export():
    all_trees: dict[str, nx.DiGraph] = {}
    pop_labels = [POP_LABELS[p] for p in POP_ORDER]

    for gene, meta in GENES.items():
        filepath = os.path.join(SCRATCH_DIR, meta["file"])
        if not os.path.exists(filepath):
            print(f"  ⚠  {meta['file']} not found — skipping {gene}")
            continue

        print(f"\n── {gene} ─────────────────────────────────────")
        variants = load_variants(filepath)
        pop_afs  = compute_pop_afs(variants)
        mat      = distance_matrix(pop_afs)
        newick   = upgma(mat, pop_labels)

        # AF dict keyed by LABEL (for nx node lookup)
        label_to_af = {POP_LABELS[p]: pop_afs[p] for p in POP_ORDER}

        G = newick_to_nx(newick, label_to_af)
        # Annotate graph-level metadata
        G.graph.update({
            "gene":        gene,
            "description": meta["description"],
            "n_variants":  len(variants),
            "source":      "gnomAD v4",
            "method":      "UPGMA on pairwise allele-frequency distance",
        })

        print(f"  Nodes : {G.number_of_nodes()} "
              f"({sum(1 for _,d in G.nodes(data=True) if d['is_leaf'])} leaves)")
        print(f"  Edges : {G.number_of_edges()}")
        print(f"  Root  : {G.graph['root']}")

        # ── Export ────────────────────────────────────────────
        safe = gene.replace("-", "_")

        # GraphML
        gml_path = os.path.join(OUT_DIR, f"{safe}.graphml")
        nx.write_graphml(G, gml_path)
        print(f"  → {gml_path}")

        # Node-link JSON  (networkx.node_link_data)
        json_path = os.path.join(OUT_DIR, f"{safe}.json")
        with open(json_path, "w") as f:
            json.dump(json_graph.node_link_data(G), f, indent=2)
        print(f"  → {json_path}")

        # Pickle
        pkl_path = os.path.join(OUT_DIR, f"{safe}.pkl")
        with open(pkl_path, "wb") as f:
            pickle.dump(G, f)
        print(f"  → {pkl_path}")

        all_trees[gene] = G

    # ── Combined pickle ───────────────────────────────────────
    combined = os.path.join(OUT_DIR, "all_trees.pkl")
    with open(combined, "wb") as f:
        pickle.dump(all_trees, f)
    print(f"\n✓ All trees written to {OUT_DIR}/")
    print(f"✓ Combined pickle: {combined}")

    return all_trees


# ════════════════════════════════════════════════════════════════
#  Quick demo
# ════════════════════════════════════════════════════════════════

def demo(all_trees: dict[str, nx.DiGraph]):
    print("\n── Quick demo ──────────────────────────────────────────")
    for gene, G in all_trees.items():
        root = G.graph["root"]
        leaves = [n for n, d in G.nodes(data=True) if d["is_leaf"]]
        max_af_leaf = max(leaves, key=lambda n: G.nodes[n].get("af") or 0)
        print(
            f"  {gene:8s} | {G.number_of_nodes()} nodes | "
            f"highest AF leaf: {G.nodes[max_af_leaf]['label']!r:25s} "
            f"(af={G.nodes[max_af_leaf]['af']:.5f})"
        )

    # Example: load from GraphML
    print("\n── Example: reload LCT from GraphML ────────────────────")
    G2 = nx.read_graphml(os.path.join(OUT_DIR, "LCT.graphml"))
    print(f"  Reloaded LCT: {G2.number_of_nodes()} nodes, "
          f"{G2.number_of_edges()} edges")

    # Example: shortest path between two leaf populations
    print("\n── Example: path African → Finnish in LCT tree ─────────")
    G_lct = all_trees["LCT"]
    # Convert to undirected for path queries
    G_und = G_lct.to_undirected()
    src = next(n for n, d in G_lct.nodes(data=True) if d.get("label") == "African")
    tgt = next(n for n, d in G_lct.nodes(data=True) if d.get("label") == "Finnish")
    path = nx.shortest_path(G_und, src, tgt)
    path_labels = [G_lct.nodes[n].get("label") or n for n in path]
    print(f"  Path: {' → '.join(path_labels)}")


if __name__ == "__main__":
    all_trees = build_and_export()
    demo(all_trees)
