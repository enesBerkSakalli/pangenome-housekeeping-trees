"""
Build population-expanded variant trees as NetworkX objects.

Each gene tree uses the selected variant leaves from ``variant_3d_analysis.py``.
Every leaf is assigned to the population where its allele frequency is highest,
variants are clustered within that population, and the eight population subtree
roots are connected by a top-level UPGMA tree of population centroids.
"""

from __future__ import annotations

import csv
import json
import os
import pickle
import re
from dataclasses import dataclass, field
from statistics import mean

import networkx as nx
import numpy as np
from networkx.readwrite import json_graph
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import pdist

ROOT = os.path.dirname(os.path.abspath(__file__))
INPUT_BUNDLE = os.path.join(ROOT, "outputs_3d", "networks_3d", "all_variant_trees.pkl")
OUT_DIR = os.path.join(ROOT, "outputs_3d", "networks_population_expanded")
SUMMARY_PATH = os.path.join(ROOT, "outputs_3d", "population_expanded_tree_summary.csv")

GENE_ORDER = ["GAPDH", "ACTB", "RPLP0"]
POP_ORDER = [
    "African",
    "American",
    "Ashkenazi Jewish",
    "East Asian",
    "Finnish",
    "Middle Eastern",
    "Non-Finnish European",
    "South Asian",
]

POP_COLORS = {
    "African": "#FF6B6B",
    "American": "#F7B267",
    "Ashkenazi Jewish": "#F4E285",
    "East Asian": "#5CC8FF",
    "Finnish": "#7BD389",
    "Middle Eastern": "#B388EB",
    "Non-Finnish European": "#4D96FF",
    "South Asian": "#FF9F9F",
}

AF_KEYS = [f"af_{pop.lower().replace(' ', '_').replace('-', '_')}" for pop in POP_ORDER]


@dataclass
class VariantLeaf:
    source_node: str
    variant_id: str
    rsid: str
    consequence: str
    consequence_group: str
    af_vector: list[float]
    assigned_population: str
    assigned_af: float
    score: float
    mean_af: float
    max_af: float
    spread: float


@dataclass
class TreeNode:
    node_id: str
    label: str
    children: list["TreeNode"] = field(default_factory=list)
    leaves: list[VariantLeaf] = field(default_factory=list)
    branch_length: float = 0.0
    height: float = 0.0
    population_label: str | None = None
    is_population_root: bool = False
    payload: VariantLeaf | None = None

    @property
    def is_leaf(self) -> bool:
        return self.payload is not None

    @property
    def taxon_count(self) -> int:
        return len(self.leaves)


def safe_id(value: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9_]+", "_", value)
    return clean.strip("_")


def load_source_graphs() -> dict[str, nx.DiGraph]:
    with open(INPUT_BUNDLE, "rb") as fh:
        return pickle.load(fh)


def extract_variant_leaves(graph: nx.DiGraph) -> list[VariantLeaf]:
    leaves = []
    for node, attrs in graph.nodes(data=True):
        if not attrs.get("is_leaf"):
            continue
        af_vector = [float(attrs.get(key, 0.0)) for key in AF_KEYS]
        pop_index = int(np.argmax(af_vector))
        leaves.append(
            VariantLeaf(
                source_node=str(node),
                variant_id=str(attrs.get("variant_id", attrs.get("label", node))),
                rsid=str(attrs.get("rsid", "")),
                consequence=str(attrs.get("consequence", "")),
                consequence_group=str(attrs.get("consequence_group", "Background")),
                af_vector=af_vector,
                assigned_population=POP_ORDER[pop_index],
                assigned_af=af_vector[pop_index],
                score=float(attrs.get("score", 0.0)),
                mean_af=float(attrs.get("mean_af", mean(af_vector) if af_vector else 0.0)),
                max_af=float(attrs.get("max_af", max(af_vector) if af_vector else 0.0)),
                spread=float(attrs.get("spread", max(af_vector) - min(af_vector) if af_vector else 0.0)),
            )
        )
    leaves.sort(key=lambda item: (item.assigned_population, -item.assigned_af, item.variant_id))
    return leaves


def cluster_variants(gene: str, population: str, leaves: list[VariantLeaf]) -> TreeNode:
    prefix = f"{gene}_{safe_id(population)}"
    ordered = sorted(leaves, key=lambda item: (-item.assigned_af, item.variant_id))
    if len(ordered) == 1:
        leaf = ordered[0]
        return TreeNode(
            node_id=f"{prefix}_leaf_{safe_id(leaf.source_node)}",
            label=leaf.variant_id,
            leaves=[leaf],
            population_label=population,
            is_population_root=True,
            payload=leaf,
        )

    matrix = np.array([leaf.af_vector for leaf in ordered], dtype=float)
    condensed = pdist(matrix, metric="euclidean")
    if not np.any(condensed):
        condensed = condensed + 1e-9
    linkage_matrix = linkage(condensed, method="average")
    scipy_root = to_tree(linkage_matrix, rd=False)

    leaf_nodes = {
        idx: TreeNode(
            node_id=f"{prefix}_leaf_{safe_id(leaf.source_node)}",
            label=leaf.variant_id,
            leaves=[leaf],
            population_label=population,
            payload=leaf,
        )
        for idx, leaf in enumerate(ordered)
    }

    def build(node) -> TreeNode:
        if node.is_leaf():
            return leaf_nodes[node.id]
        left = build(node.get_left())
        right = build(node.get_right())
        node_height = float(node.dist) / 2.0
        left.branch_length = max(1e-9, node_height - left.height)
        right.branch_length = max(1e-9, node_height - right.height)
        return TreeNode(
            node_id=f"{prefix}_internal_{node.id}",
            label=f"{population} clade {node.id}",
            children=[left, right],
            leaves=left.leaves + right.leaves,
            height=node_height,
            population_label=population,
        )

    root = build(scipy_root)
    root.is_population_root = True
    root.label = f"{population} expanded subtree"
    root.branch_length = 0.0
    return root


def population_centroid(leaves: list[VariantLeaf]) -> list[float]:
    matrix = np.array([leaf.af_vector for leaf in leaves], dtype=float)
    return list(np.mean(matrix, axis=0))


def cluster_population_roots(gene: str, population_roots: dict[str, TreeNode]) -> TreeNode:
    ordered_pops = [pop for pop in POP_ORDER if pop in population_roots]
    matrix = np.array([population_centroid(population_roots[pop].leaves) for pop in ordered_pops], dtype=float)
    condensed = pdist(matrix, metric="euclidean")
    if not np.any(condensed):
        condensed = condensed + 1e-9
    linkage_matrix = linkage(condensed, method="average")
    scipy_root = to_tree(linkage_matrix, rd=False)

    def build(node) -> TreeNode:
        if node.is_leaf():
            return population_roots[ordered_pops[node.id]]
        left = build(node.get_left())
        right = build(node.get_right())
        node_height = max(left.height, right.height) + max(float(node.dist) / 2.0, 1e-9)
        left.branch_length = max(1e-9, node_height - left.height)
        right.branch_length = max(1e-9, node_height - right.height)
        return TreeNode(
            node_id=f"{gene}_top_internal_{node.id}",
            label=f"{gene} population split {node.id}",
            children=[left, right],
            leaves=left.leaves + right.leaves,
            height=node_height,
        )

    root = build(scipy_root)
    root.node_id = f"{gene}_root"
    root.label = f"{gene} population-expanded tree"
    root.branch_length = 0.0
    return root


def add_tree_to_graph(graph: nx.DiGraph, node: TreeNode, parent: str | None = None) -> None:
    population = node.population_label or ""
    payload = {
        "label": node.label,
        "is_leaf": node.is_leaf,
        "taxon_count": int(node.taxon_count),
        "branch_length": float(node.branch_length),
        "cluster_height": float(node.height),
        "population_label": population,
        "population_color": POP_COLORS.get(population, ""),
        "is_population_root": bool(node.is_population_root),
    }
    if node.is_leaf and node.payload is not None:
        leaf = node.payload
        payload.update(
            {
                "variant_id": leaf.variant_id,
                "rsid": leaf.rsid,
                "consequence": leaf.consequence,
                "consequence_group": leaf.consequence_group,
                "assigned_population": leaf.assigned_population,
                "assigned_af": float(leaf.assigned_af),
                "score": float(leaf.score),
                "mean_af": float(leaf.mean_af),
                "max_af": float(leaf.max_af),
                "spread": float(leaf.spread),
            }
        )
        for pop, value in zip(POP_ORDER, leaf.af_vector):
            payload[f"af_{pop.lower().replace(' ', '_').replace('-', '_')}"] = float(value)
    graph.add_node(node.node_id, **payload)
    if parent is not None:
        graph.add_edge(parent, node.node_id, weight=float(node.branch_length))
    for child in node.children:
        add_tree_to_graph(graph, child, node.node_id)


def build_gene_tree(gene: str, source_graph: nx.DiGraph) -> nx.DiGraph:
    leaves = extract_variant_leaves(source_graph)
    by_population = {pop: [] for pop in POP_ORDER}
    for leaf in leaves:
        by_population[leaf.assigned_population].append(leaf)

    missing = [pop for pop, pop_leaves in by_population.items() if not pop_leaves]
    if missing:
        raise RuntimeError(f"{gene} has no leaves assigned to: {', '.join(missing)}")

    population_roots = {
        pop: cluster_variants(gene, pop, pop_leaves)
        for pop, pop_leaves in by_population.items()
    }
    root = cluster_population_roots(gene, population_roots)

    graph = nx.DiGraph()
    graph.graph["gene"] = gene
    graph.graph["description"] = source_graph.graph.get("description", "")
    graph.graph["root"] = root.node_id
    graph.graph["taxon_definition"] = "variant"
    graph.graph["tree_model"] = "population-expanded constrained UPGMA"
    graph.graph["population_assignment"] = "variant assigned to population with maximum AF"
    graph.graph["selected_taxa"] = len(leaves)
    graph.graph["total_variants_in_source"] = source_graph.graph.get("total_variants_in_source", "")
    graph.graph["source_selection"] = source_graph.graph.get("selection", "")
    graph.graph["source_distance"] = source_graph.graph.get("distance", "")
    graph.graph["source_clustering"] = source_graph.graph.get("clustering", "")
    graph.graph["populations"] = ",".join(POP_ORDER)
    add_tree_to_graph(graph, root)
    return graph


def export_graphs(graphs: dict[str, nx.DiGraph]) -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "all_population_expanded_trees.pkl"), "wb") as fh:
        pickle.dump(graphs, fh)

    with open(SUMMARY_PATH, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["gene", "nodes", "edges", "leaves", *POP_ORDER])
        for gene, graph in graphs.items():
            counts = {
                pop: sum(
                    1
                    for _, attrs in graph.nodes(data=True)
                    if attrs.get("is_leaf") and attrs.get("assigned_population") == pop
                )
                for pop in POP_ORDER
            }
            writer.writerow(
                [
                    gene,
                    graph.number_of_nodes(),
                    graph.number_of_edges(),
                    sum(1 for _, attrs in graph.nodes(data=True) if attrs.get("is_leaf")),
                    *[counts[pop] for pop in POP_ORDER],
                ]
            )

    for gene, graph in graphs.items():
        safe_gene = safe_id(gene)
        with open(os.path.join(OUT_DIR, f"{safe_gene}.pkl"), "wb") as fh:
            pickle.dump(graph, fh)
        with open(os.path.join(OUT_DIR, f"{safe_gene}.json"), "w") as fh:
            json.dump(json_graph.node_link_data(graph), fh, indent=2)
        nx.write_graphml(graph, os.path.join(OUT_DIR, f"{safe_gene}.graphml"))


def main() -> None:
    source_graphs = load_source_graphs()
    graphs = {}
    for gene in GENE_ORDER:
        graph = build_gene_tree(gene, source_graphs[gene])
        graphs[gene] = graph
        leaves = sum(1 for _, attrs in graph.nodes(data=True) if attrs.get("is_leaf"))
        pop_roots = sum(1 for _, attrs in graph.nodes(data=True) if attrs.get("is_population_root"))
        print(f"{gene}: leaves={leaves} nodes={graph.number_of_nodes()} population_subtrees={pop_roots}")
    export_graphs(graphs)
    print(f"Wrote population-expanded NetworkX trees to {OUT_DIR}")


if __name__ == "__main__":
    main()
