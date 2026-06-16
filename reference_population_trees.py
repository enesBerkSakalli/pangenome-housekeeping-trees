"""
Build shared-taxa population trees for neutral reference genes.

Genes:
- GAPDH
- ACTB
- RPLP0

Taxa are the same across all trees: the 8 gnomAD population groups.
"""

from __future__ import annotations

import json
import math
import os
import pickle
from dataclasses import dataclass

import networkx as nx
import numpy as np
from networkx.readwrite import json_graph
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRATCH_DIR = os.path.join(SCRIPT_DIR, "scratch")
OUT_DIR = os.path.join(SCRIPT_DIR, "outputs_reference")
NETWORKS_DIR = os.path.join(OUT_DIR, "networks_reference")

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
    "GAPDH": {"file": "GAPDH_variants.json", "description": "GAPDH (glycolytic housekeeping reference)"},
    "ACTB": {"file": "ACTB_variants.json", "description": "ACTB (beta-actin structural reference)"},
    "RPLP0": {"file": "RPLP0_variants.json", "description": "RPLP0 (ribosomal housekeeping reference)"},
}


@dataclass
class ClusterNode:
    node_id: str
    label: str
    members: list[str]
    children: list["ClusterNode"]
    branch_length: float
    height: float

    @property
    def size(self) -> int:
        return len(self.members)

    @property
    def is_leaf(self) -> bool:
        return not self.children


def load_variants(path: str) -> list[dict]:
    with open(path) as fh:
        data = json.load(fh)
    return data["data"]["gene"]["variants"]


def compute_pop_afs(variants: list[dict]) -> dict[str, float]:
    ac = {pop: 0 for pop in POP_ORDER}
    an = {pop: 0 for pop in POP_ORDER}
    for variant in variants:
        for source_name in ("exome", "genome"):
            source = variant.get(source_name) or {}
            for pop in source.get("populations", []):
                pop_id = pop.get("id")
                if pop_id in ac:
                    ac[pop_id] += pop.get("ac", 0)
                    an[pop_id] += pop.get("an", 0)
    return {pop: (ac[pop] / an[pop] if an[pop] else 0.0) for pop in POP_ORDER}


def distance_matrix(pop_afs: dict[str, float]) -> np.ndarray:
    matrix = np.zeros((len(POP_ORDER), len(POP_ORDER)), dtype=float)
    for row, left in enumerate(POP_ORDER):
        for col, right in enumerate(POP_ORDER):
            if row < col:
                distance = abs(pop_afs[left] - pop_afs[right])
                matrix[row, col] = distance
                matrix[col, row] = distance
    return matrix


def build_cluster_tree(pop_afs: dict[str, float]) -> ClusterNode:
    matrix = distance_matrix(pop_afs)
    condensed = squareform(matrix)
    linkage_matrix = linkage(condensed, method="average")
    scipy_root = to_tree(linkage_matrix, rd=False)

    leaves = {
        idx: ClusterNode(
            node_id=POP_ORDER[idx],
            label=POP_LABELS[POP_ORDER[idx]],
            members=[POP_ORDER[idx]],
            children=[],
            branch_length=0.0,
            height=0.0,
        )
        for idx in range(len(POP_ORDER))
    }

    def build(node) -> ClusterNode:
        if node.is_leaf():
            return leaves[node.id]
        left = build(node.get_left())
        right = build(node.get_right())
        node_height = float(node.dist) / 2.0
        left.branch_length = max(0.0, node_height - left.height)
        right.branch_length = max(0.0, node_height - right.height)
        return ClusterNode(
            node_id=f"internal_{node.id}",
            label=f"clade_{node.id}",
            members=left.members + right.members,
            children=[left, right],
            branch_length=0.0,
            height=node_height,
        )

    return build(scipy_root)


def cluster_to_graph(root: ClusterNode, gene: str, description: str, pop_afs: dict[str, float]) -> nx.DiGraph:
    graph = nx.DiGraph()
    graph.graph["gene"] = gene
    graph.graph["description"] = description
    graph.graph["root"] = root.node_id
    graph.graph["taxon_definition"] = "population"
    graph.graph["taxa_order"] = [POP_LABELS[pop] for pop in POP_ORDER]
    graph.graph["pop_afs"] = {POP_LABELS[pop]: float(pop_afs[pop]) for pop in POP_ORDER}

    stack = [root]
    while stack:
        node = stack.pop()
        payload = {
            "label": node.label,
            "is_leaf": node.is_leaf,
            "taxon_count": node.size,
            "branch_length": float(node.branch_length),
            "cluster_height": float(node.height),
        }
        if node.is_leaf:
            pop_id = node.members[0]
            payload["population_id"] = pop_id
            payload["population_label"] = POP_LABELS[pop_id]
            payload["af"] = float(pop_afs[pop_id])
        graph.add_node(node.node_id, **payload)
        for child in node.children:
            graph.add_edge(node.node_id, child.node_id, weight=float(child.branch_length))
            stack.append(child)
    return graph


def export_graphs(graphs: dict[str, nx.DiGraph]) -> None:
    os.makedirs(NETWORKS_DIR, exist_ok=True)
    with open(os.path.join(NETWORKS_DIR, "all_reference_trees.pkl"), "wb") as fh:
        pickle.dump(graphs, fh)

    summary = []
    for gene, graph in graphs.items():
        safe = gene.replace("-", "_")
        with open(os.path.join(NETWORKS_DIR, f"{safe}.pkl"), "wb") as fh:
            pickle.dump(graph, fh)
        with open(os.path.join(NETWORKS_DIR, f"{safe}.json"), "w") as fh:
            json.dump(json_graph.node_link_data(graph), fh, indent=2)
        graphml_graph = graph.copy()
        graphml_graph.graph["taxa_order"] = ",".join(graph.graph["taxa_order"])
        graphml_graph.graph["pop_afs"] = json.dumps(graph.graph["pop_afs"], sort_keys=True)
        nx.write_graphml(graphml_graph, os.path.join(NETWORKS_DIR, f"{safe}.graphml"))
        summary.append(
            {
                "gene": gene,
                "nodes": graph.number_of_nodes(),
                "edges": graph.number_of_edges(),
                "taxa": sum(1 for _, attrs in graph.nodes(data=True) if attrs.get("is_leaf")),
            }
        )

    with open(os.path.join(OUT_DIR, "reference_tree_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    graphs = {}
    for gene, meta in GENES.items():
        variants = load_variants(os.path.join(SCRATCH_DIR, meta["file"]))
        pop_afs = compute_pop_afs(variants)
        root = build_cluster_tree(pop_afs)
        graph = cluster_to_graph(root, gene, meta["description"], pop_afs)
        graphs[gene] = graph
        print(f"{gene}: taxa=8 nodes={graph.number_of_nodes()} edges={graph.number_of_edges()}")

    export_graphs(graphs)
    print(f"Wrote shared-taxa trees to {NETWORKS_DIR}")


if __name__ == "__main__":
    main()
