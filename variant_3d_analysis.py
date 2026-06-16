"""
Redo the pangenome analysis at the variant level and export 3D layered trees.

Taxa definition
---------------
Each taxon is a single variant from the raw gnomAD gene dump. For each gene, the
script selects the most population-informative variants and clusters them by
their 8-population allele-frequency profiles.

Outputs
-------
outputs_3d/
  <GENE>_variant_tree.json
  all_variant_trees.json
  variant_tree_summary.csv

variant_tree_data.js
  Browser-ready bundle for the 3D viewer.
"""

from __future__ import annotations

import csv
import json
import math
import os
from collections import Counter, deque
from dataclasses import dataclass
from statistics import mean, pvariance

import networkx as nx
import numpy as np
from networkx.readwrite import json_graph
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import pdist

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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRATCH_DIR = os.path.join(SCRIPT_DIR, "scratch")
OUT_DIR = os.path.join(SCRIPT_DIR, "outputs_3d")
NETWORKS_DIR = os.path.join(OUT_DIR, "networks_3d")
JS_BUNDLE_PATH = os.path.join(SCRIPT_DIR, "variant_tree_data.js")

TARGET_TAXA = 1024
MIN_TAXA = 250
Z_LAYER_COUNT = 3

LAYER_LABELS = {
    0: "Layer 1 · trunk clades",
    1: "Layer 2 · mid branches",
    2: "Layer 3 · taxon frontier",
}

CONSEQUENCE_GROUPS = {
    "Protein-altering": "#ff7a59",
    "Regulatory/splice": "#ffd166",
    "Background": "#4cc9f0",
    "Internal": "#94a3b8",
}


@dataclass
class VariantTaxon:
    taxon_id: str
    variant_id: str
    rsid: str
    consequence: str
    consequence_group: str
    af_vector: list[float]
    mean_af: float
    max_af: float
    spread: float
    variance: float
    score: float
    active_populations: int


@dataclass
class ClusterNode:
    node_id: str
    label: str
    members: list[VariantTaxon]
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
    try:
        return data["data"]["gene"]["variants"]
    except (KeyError, TypeError):
        try:
            return data["data"]["region"]["variants"]
        except (KeyError, TypeError):
            return []


def consequence_group(raw: str) -> str:
    text = (raw or "").lower()
    protein_terms = (
        "missense",
        "frameshift",
        "stop_gained",
        "stop_lost",
        "start_lost",
        "protein_altering",
        "inframe",
        "coding_sequence",
    )
    regulatory_terms = (
        "splice",
        "utr",
        "upstream",
        "downstream",
        "promoter",
        "regulatory",
        "non_coding_transcript",
    )
    if any(term in text for term in protein_terms):
        return "Protein-altering"
    if any(term in text for term in regulatory_terms):
        return "Regulatory/splice"
    return "Background"


def variant_vector(variant: dict) -> list[float]:
    ac = {pop: 0 for pop in POP_ORDER}
    an = {pop: 0 for pop in POP_ORDER}
    for source_name in ("exome", "genome"):
        source = variant.get(source_name) or {}
        for pop in source.get("populations", []):
            pop_id = pop.get("id")
            if pop_id in ac:
                ac[pop_id] += pop.get("ac", 0)
                an[pop_id] += pop.get("an", 0)
    return [ac[pop] / an[pop] if an[pop] else 0.0 for pop in POP_ORDER]


def build_taxon(variant: dict, index: int) -> VariantTaxon | None:
    vec = variant_vector(variant)
    if max(vec) <= 0:
        return None
    active = sum(1 for value in vec if value > 0)
    mean_af = mean(vec)
    spread = max(vec) - min(vec)
    variance = pvariance(vec)
    score = (spread * 0.55) + (math.sqrt(variance) * 0.35) + (mean_af * 0.10)
    consequence = variant.get("consequence", "unknown")
    variant_id = variant.get("variant_id", f"variant_{index}")
    rsids = variant.get("rsids") or []
    rsid = rsids[0] if rsids else ""
    return VariantTaxon(
        taxon_id=f"taxon_{index:04d}",
        variant_id=variant_id,
        rsid=rsid,
        consequence=consequence,
        consequence_group=consequence_group(consequence),
        af_vector=vec,
        mean_af=mean_af,
        max_af=max(vec),
        spread=spread,
        variance=variance,
        score=score,
        active_populations=active,
    )


def select_taxa(variants: list[dict], target_taxa: int = TARGET_TAXA) -> list[VariantTaxon]:
    taxa = []
    for index, variant in enumerate(variants):
        taxon = build_taxon(variant, index)
        if taxon is not None:
            taxa.append(taxon)

    taxa.sort(
        key=lambda item: (
            item.active_populations >= 2,
            item.score,
            item.spread,
            item.max_af,
        ),
        reverse=True,
    )
    return taxa[:target_taxa]


def upgma(taxa: list[VariantTaxon]) -> ClusterNode:
    if len(taxa) == 1:
        taxon = taxa[0]
        return ClusterNode(
            node_id=taxon.taxon_id,
            label=taxon.variant_id,
            members=[taxon],
            children=[],
            branch_length=0.0,
            height=0.0,
        )

    matrix = np.array([taxon.af_vector for taxon in taxa], dtype=float)
    condensed = pdist(matrix, metric="euclidean")
    linkage_matrix = linkage(condensed, method="average")
    scipy_root = to_tree(linkage_matrix, rd=False)

    leaves = {
        idx: ClusterNode(
            node_id=taxon.taxon_id,
            label=taxon.variant_id,
            members=[taxon],
            children=[],
            branch_length=0.0,
            height=0.0,
        )
        for idx, taxon in enumerate(taxa)
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


def leaf_order(node: ClusterNode, acc: list[ClusterNode]) -> list[ClusterNode]:
    if node.is_leaf:
        acc.append(node)
        return acc
    for child in node.children:
        leaf_order(child, acc)
    return acc


def structural_layers(root: ClusterNode) -> tuple[dict[str, int], dict[str, int]]:
    depth_map = {root.node_id: 0}
    queue = deque([root])
    nodes_by_id = {root.node_id: root}
    while queue:
        node = queue.popleft()
        for child in node.children:
            depth_map[child.node_id] = depth_map[node.node_id] + 1
            nodes_by_id[child.node_id] = child
            queue.append(child)
    max_depth = max(depth_map.values()) or 1
    z_layers = {}
    for node_id, depth in depth_map.items():
        tier = min(Z_LAYER_COUNT - 1, int((depth / max_depth) * Z_LAYER_COUNT))
        if tier == Z_LAYER_COUNT:
            tier = Z_LAYER_COUNT - 1
        z_layers[node_id] = tier
    return depth_map, z_layers


def dominant_group(node: ClusterNode) -> str:
    if node.is_leaf:
        return node.members[0].consequence_group
    counts = Counter(member.consequence_group for member in node.members)
    return counts.most_common(1)[0][0]


def cluster_to_networkx(
    root: ClusterNode,
    gene: str,
    description: str,
    selected_taxa: int,
    total_variants_in_source: int,
) -> nx.DiGraph:
    graph = nx.DiGraph()
    graph.graph["gene"] = gene
    graph.graph["description"] = description
    graph.graph["root"] = root.node_id
    graph.graph["taxon_definition"] = "variant"
    graph.graph["selected_taxa"] = selected_taxa
    graph.graph["total_variants_in_source"] = total_variants_in_source
    graph.graph["selection"] = f"top {selected_taxa} population-informative variants by AF spread/variance"
    graph.graph["distance"] = "euclidean distance on 8-population AF vectors"
    graph.graph["clustering"] = "UPGMA"

    queue = deque([root])
    while queue:
        node = queue.popleft()
        group = dominant_group(node)
        payload = {
            "label": node.label,
            "is_leaf": node.is_leaf,
            "taxon_count": node.size,
            "branch_length": float(node.branch_length),
            "cluster_height": float(node.height),
            "dominant_group": group,
        }
        if node.is_leaf:
            taxon = node.members[0]
            payload.update(
                {
                    "variant_id": taxon.variant_id,
                    "rsid": taxon.rsid,
                    "consequence": taxon.consequence,
                    "consequence_group": taxon.consequence_group,
                    "active_populations": taxon.active_populations,
                    "score": float(taxon.score),
                    "mean_af": float(taxon.mean_af),
                    "max_af": float(taxon.max_af),
                    "spread": float(taxon.spread),
                    "variance": float(taxon.variance),
                }
            )
            for pop_label, value in zip((POP_LABELS[p] for p in POP_ORDER), taxon.af_vector):
                payload[f"af_{pop_label.lower().replace(' ', '_').replace('-', '_')}"] = float(value)
        graph.add_node(node.node_id, **payload)
        for child in node.children:
            graph.add_edge(node.node_id, child.node_id, weight=float(child.branch_length))
            queue.append(child)
    return graph


def build_layout(root: ClusterNode) -> tuple[dict[str, dict], list[dict]]:
    leaves = leaf_order(root, [])
    leaf_positions = {leaf.node_id: idx for idx, leaf in enumerate(leaves)}
    depth_map, z_tiers = structural_layers(root)
    max_depth = max(depth_map.values()) or 1
    max_height = max(node.height for node in walk_nodes(root)) or 1.0

    nodes = {}
    edges = []

    def assign(node: ClusterNode) -> tuple[float, float]:
        child_results = []
        if node.is_leaf:
            x = float(leaf_positions[node.node_id])
        else:
            child_results = [(child, assign(child)) for child in node.children]
            x = sum(coord[0] for _, coord in child_results) / len(child_results)
        y = node.height
        z_tier = z_tiers[node.node_id]
        z = (z_tier - 1) * 180.0

        leaf_payload = {}
        if node.is_leaf:
            taxon = node.members[0]
            leaf_payload = {
                "variant_id": taxon.variant_id,
                "rsid": taxon.rsid,
                "consequence": taxon.consequence,
                "consequence_group": taxon.consequence_group,
                "af_vector": dict(zip((POP_LABELS[p] for p in POP_ORDER), [round(v, 8) for v in taxon.af_vector])),
                "score": round(taxon.score, 8),
                "active_populations": taxon.active_populations,
            }

        nodes[node.node_id] = {
            "id": node.node_id,
            "label": node.label,
            "is_leaf": node.is_leaf,
            "taxon_count": node.size,
            "x": round((x / max(1, len(leaves) - 1)) * 860.0 - 430.0, 4),
            "y": round((y / max_height) * 520.0 - 260.0, 4),
            "z": z,
            "depth": depth_map[node.node_id],
            "layer3d": z_tier,
            "layer_label": LAYER_LABELS[z_tier],
            "branch_length": round(node.branch_length, 8),
            "cluster_height": round(node.height, 8),
            "color": CONSEQUENCE_GROUPS["Internal"] if not node.is_leaf else CONSEQUENCE_GROUPS[leaf_payload["consequence_group"]],
            "dominant_group": dominant_group(node),
            **leaf_payload,
        }

        for child, _child_coords in child_results:
            child_node = nodes[child.node_id]
            edges.append(
                {
                    "source": node.node_id,
                    "target": child.node_id,
                    "length": round(child.branch_length, 8),
                    "points": [
                        [nodes[node.node_id]["x"], nodes[node.node_id]["y"], nodes[node.node_id]["z"]],
                        [child_node["x"], child_node["y"], child_node["z"]],
                    ],
                }
            )
        return x, y

    assign(root)
    return nodes, edges


def walk_nodes(root: ClusterNode) -> list[ClusterNode]:
    acc = []
    stack = [root]
    while stack:
        node = stack.pop()
        acc.append(node)
        stack.extend(reversed(node.children))
    return acc


def to_newick(node: ClusterNode) -> str:
    if node.is_leaf:
        safe_label = node.label.replace(":", "_").replace(",", "_")
        return f"{safe_label}:{node.branch_length:.8f}"
    inner = ",".join(to_newick(child) for child in node.children)
    return f"({inner}){node.label}:{node.branch_length:.8f}"


def analyse_gene(gene: str, meta: dict) -> dict:
    variants = load_variants(os.path.join(SCRATCH_DIR, meta["file"]))
    taxa = select_taxa(variants, TARGET_TAXA)
    if len(taxa) < MIN_TAXA:
        raise RuntimeError(f"{gene} yielded only {len(taxa)} taxa; minimum required is {MIN_TAXA}")

    root = upgma(taxa)
    graph = cluster_to_networkx(
        root=root,
        gene=gene,
        description=meta["description"],
        selected_taxa=len(taxa),
        total_variants_in_source=len(variants),
    )
    nodes, edges = build_layout(root)
    consequence_counts = Counter(taxon.consequence_group for taxon in taxa)

    return {
        "gene": gene,
        "description": meta["description"],
        "method": {
            "taxon_definition": "variant",
            "selection": f"top {len(taxa)} population-informative variants by AF spread/variance",
            "distance": "euclidean distance on 8-population AF vectors",
            "clustering": "UPGMA",
            "layers3d": list(LAYER_LABELS.values()),
        },
        "stats": {
            "total_variants_in_source": len(variants),
            "selected_taxa": len(taxa),
            "node_count": len(nodes),
            "edge_count": len(edges),
            "max_branch_height": max(node["cluster_height"] for node in nodes.values()),
            "consequence_groups": dict(consequence_counts),
        },
        "populations": [POP_LABELS[pop] for pop in POP_ORDER],
        "graph": graph,
        "nodes": nodes,
        "edges": edges,
        "newick": to_newick(root) + ";",
    }


def export_outputs(results: dict[str, dict]) -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(NETWORKS_DIR, exist_ok=True)

    summary_path = os.path.join(OUT_DIR, "variant_tree_summary.csv")
    with open(summary_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "gene",
                "selected_taxa",
                "node_count",
                "edge_count",
                "total_variants_in_source",
                "protein_altering_taxa",
                "regulatory_splice_taxa",
                "background_taxa",
            ]
        )
        for gene, payload in results.items():
            counts = payload["stats"]["consequence_groups"]
            writer.writerow(
                [
                    gene,
                    payload["stats"]["selected_taxa"],
                    payload["stats"]["node_count"],
                    payload["stats"]["edge_count"],
                    payload["stats"]["total_variants_in_source"],
                    counts.get("Protein-altering", 0),
                    counts.get("Regulatory/splice", 0),
                    counts.get("Background", 0),
                ]
            )

    serializable_results = {}
    networkx_bundle = {}
    for gene, payload in results.items():
        graph = payload["graph"]
        networkx_bundle[gene] = graph
        payload_copy = dict(payload)
        payload_copy["graph"] = json_graph.node_link_data(graph)
        serializable_results[gene] = payload_copy

    combined_path = os.path.join(OUT_DIR, "all_variant_trees.json")
    with open(combined_path, "w") as fh:
        json.dump(serializable_results, fh, indent=2)

    for gene, payload in serializable_results.items():
        per_gene_path = os.path.join(OUT_DIR, f"{gene.replace('-', '_')}_variant_tree.json")
        with open(per_gene_path, "w") as fh:
            json.dump(payload, fh, indent=2)

    with open(JS_BUNDLE_PATH, "w") as fh:
        fh.write("window.PANGENOME_VARIANT_TREES = ")
        json.dump(serializable_results, fh, indent=2)
        fh.write(";\n")

    for gene, graph in networkx_bundle.items():
        safe_gene = gene.replace("-", "_")
        with open(os.path.join(NETWORKS_DIR, f"{safe_gene}.pkl"), "wb") as fh:
            import pickle

            pickle.dump(graph, fh)
        with open(os.path.join(NETWORKS_DIR, f"{safe_gene}.json"), "w") as fh:
            json.dump(json_graph.node_link_data(graph), fh, indent=2)
        nx.write_graphml(graph, os.path.join(NETWORKS_DIR, f"{safe_gene}.graphml"))

    with open(os.path.join(NETWORKS_DIR, "all_variant_trees.pkl"), "wb") as fh:
        import pickle

        pickle.dump(networkx_bundle, fh)


def main() -> None:
    results = {}
    for gene, meta in GENES.items():
        print(f"Analyzing {gene}...", flush=True)
        payload = analyse_gene(gene, meta)
        results[gene] = payload
        print(
            f"  selected_taxa={payload['stats']['selected_taxa']} "
            f"nodes={payload['stats']['node_count']} edges={payload['stats']['edge_count']}"
            ,
            flush=True,
        )

    export_outputs(results)
    print(f"\nWrote 3D analysis outputs to {OUT_DIR}", flush=True)
    print(f"Wrote browser bundle to {JS_BUNDLE_PATH}", flush=True)


if __name__ == "__main__":
    main()
