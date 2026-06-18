"""
Create a DataDiVR project from the inferred stacked housekeeping-gene trees.

The project is built as DataDiVR scene/layout variants over one stable
NetworkX graph:
- scene 00 shows all three stacked gene trees
- scenes 01-03 keep the same graph and coordinates but focus one gene at a time

No dependency installation is required. This script deliberately avoids the
DataDiVR uploader dependency chain and writes the project files/textures
directly in the format used by DataDiVR_WebApp/static/projects.
"""

from __future__ import annotations

import copy
import json
import math
import os
import pickle
import shutil
import struct
import zlib
from collections import Counter
from typing import Iterable

import networkx as nx
import numpy as np

import plotly_stacked_trees as pst


ROOT = os.path.dirname(os.path.abspath(__file__))
DATADIVR_REPO = os.path.join(ROOT, "external", "DataDiVR_WebApp")
PROJECT_NAME = "Pangenome_Housekeeping_Stacked_Trees"
PROJECT_DIR = os.path.join(DATADIVR_REPO, "static", "projects", PROJECT_NAME)
PORTABLE_PROJECT_DIR = os.path.join(ROOT, "datadivr_project", PROJECT_NAME)
OUTPUT_JSON = os.path.join(ROOT, "outputs_3d", f"{PROJECT_NAME}.json")
OUTPUT_AUDIT = os.path.join(ROOT, "outputs_3d", f"{PROJECT_NAME}_datadivr_audit.json")
OUTPUT_SCENES_PICKLE = os.path.join(ROOT, "outputs_3d", f"{PROJECT_NAME}_networkx_scenes.pkl")
OUTPUT_PATHS_JSON = os.path.join(ROOT, "outputs_3d", f"{PROJECT_NAME}_paths.json")
OUTPUT_PATH_CONNECTIONS_JSON = os.path.join(ROOT, "outputs_3d", f"{PROJECT_NAME}_path_connections.json")
OUTPUT_COORDINATES_JSON = os.path.join(ROOT, "outputs_3d", "datadivr_coordinate_mappings.json")

SCENE_NAMES = [
    "00_all_genes_stacked",
    *[f"{index:02d}_{gene}_focus" for index, gene in enumerate(pst.GENE_ORDER, start=1)],
]

INFO = (
    "Three stacked NetworkX unrooted housekeeping-gene trees inferred from "
    "cached NCBI/Ensembl protein homology alignments. The same evidence-backed "
    "taxa are used in every gene layer; focus scenes fade non-selected genes."
)
PATH_GROUP_DOMINANCE_THRESHOLD = 0.60
ANCESTOR_CONTEXT_LEVELS = (1, 2)
BACKBONE_LINK_COLOR = (202, 220, 250, 70)
AMBIGUOUS_LINK_COLOR = (142, 160, 190, 48)
FADED_LINK_COLOR = (52, 60, 76, 10)
INTERLAYER_HUMAN_CONNECTOR_COLOR = (255, 242, 204, 150)
INTERLAYER_NEUTRAL_CONNECTOR_COLOR = (207, 230, 255, 128)
INTERLAYER_FLOW_NODE_COLOR = (18, 24, 34, 35)
INTERLAYER_PORT_NODE_COLOR = (36, 48, 62, 76)
INTERLAYER_SELECTED_TAXON_NODE_COLOR = (28, 34, 46, 44)
INTERLAYER_CLADE_FLOW_ALPHA = 126
INTERLAYER_SELECTED_TAXON_FLOW_ALPHA = 38
INTERLAYER_SELECTED_TAXON_FOCUS_ALPHA = 92

PATH_ANIMATION_SETTINGS = {
    "enabled": True,
    "mode": "flow_pulses",
    "drawPathCurves": True,
    "pulseEnabled": False,
    "maxVisiblePaths": 260,
    "curveSegments": 36,
    "pulseRadius": 0.072,
    "pulseSpeed": 0.18,
    "pulseStagger": 0.007,
    "curveOpacity": 0.14,
    "pulseOpacity": 0.88,
    "focusSceneSpeedBoost": 1.35,
}

SUBTREE_HIGHLIGHT_ANIMATION = {
    "enabled": True,
    "mode": "subtree_focus",
    "initialSetup": "Ray-finned fish",
    "presentation": {
        "enabled": True,
        "autoPlay": True,
        "stageSeconds": 4.0,
        "showStageLabel": True,
        "sequence": [
            {
                "mode": "overview",
                "title": "All tree layers",
                "layoutIndex": 0,
                "highlight": False,
            },
            *[
                {
                    "mode": "layer",
                    "title": f"{gene} tree layer",
                    "layoutIndex": index,
                    "highlight": False,
                }
                for index, gene in enumerate(pst.GENE_ORDER, start=1)
            ],
            *[
                {
                    "mode": "subtree_layer",
                    "title": f"Ray-finned fish subtree in {gene}",
                    "layoutIndex": index,
                    "setup": "Ray-finned fish",
                    "highlight": True,
                    "linkGroups": ["Ray-finned fish paths"],
                }
                for index, gene in enumerate(pst.GENE_ORDER, start=1)
            ],
            {
                "mode": "subtree_all",
                "title": "Ray-finned fish subtree and inter-layer connections",
                "layoutIndex": 0,
                "setup": "Ray-finned fish",
                "highlight": True,
                "linkGroups": [
                    "Ray-finned fish paths",
                    "Ray-finned fish clade-level flow corridor",
                    "Ray-finned fish selected same-taxon history tracks",
                ],
            },
        ],
    },
    "dimNonFocus": True,
    "dimNodeOpacity": 0.035,
    "dimLinkOpacity": 0.004,
    "highlightNodeOpacity": 1.0,
    "highlightLinkOpacity": 0.88,
    "linkPulseOpacity": 0.12,
    "pulseScale": 0.72,
    "pulseSpeed": 0.85,
    "setups": [
        {
            "name": "Ray-finned fish",
            "nodeGroup": "Ray-finned fish",
            "linkGroups": [
                "Ray-finned fish paths",
                "Ray-finned fish clade-level flow corridor",
                "Ray-finned fish selected same-taxon history tracks",
            ],
            "color": [34, 211, 238],
            "layouts": SCENE_NAMES,
        },
        {
            "name": "Glires",
            "nodeGroup": "Glires",
            "linkGroups": [
                "Glires paths",
                "Glires clade-level flow corridor",
                "Glires selected same-taxon history tracks",
            ],
            "color": [56, 189, 248],
            "layouts": SCENE_NAMES,
        },
    ],
}

DEFAULT_PDATA = {
    "layoutsDD": "0",
    "layoutsRGBDD": "0",
    "linksDD": "0",
    "linksRGBDD": "0",
    "selectionsDD": "0",
    "slider-label_scale": "420",
    "slider-node_size": "288",
    "slider-link_transparency": "360",
    "slider-link_size": "260",
    "activeNode": "0",
    "cbnode": [],
    "prot_size": "",
    "CGlayouts": "0",
    "CGvis": "0",
    "analytics": 0,
    "annotation-1": 0,
    "annotation-2": 0,
    "annotation-Operations": 0,
    "annotationOperationsActive": False,
}


def clamp(value: float, low: float = 0.0, high: float = 1.0) -> float:
    return max(low, min(high, float(value)))


def hex_to_rgba(hex_color: str, alpha: int) -> tuple[int, int, int, int]:
    value = hex_color.lstrip("#")
    return (
        int(value[0:2], 16),
        int(value[2:4], 16),
        int(value[4:6], 16),
        int(clamp(alpha, 0, 255)),
    )


def safe_name(value: object) -> str:
    return str(value).replace("/", "_").replace("\\", "_").replace(" ", "_")


def node_display_name(attrs: dict) -> str:
    if attrs.get("is_interlayer_clade_port"):
        return f"{attrs.get('gene_pair', 'gene pair')} {attrs.get('interlayer_group', 'clade')} flow port"
    if attrs.get("is_interlayer_clade_control"):
        return f"{attrs.get('gene_pair', 'gene pair')} {attrs.get('interlayer_group', 'clade')} flow corridor"
    if attrs.get("is_interlayer_taxon_control"):
        return f"{attrs.get('interlayer_label', 'selected taxon')} history-flow control"
    if attrs.get("species_label"):
        return str(attrs["species_label"])
    if attrs.get("is_human_sequence_anchor") or attrs.get("is_human_variant_context"):
        return "Human subtree"
    if attrs.get("is_human_population_summary"):
        return "Human subtree leaf"
    return str(attrs.get("label") or attrs.get("name") or "internal split")


def node_annotation(gene: str, original_node: str, attrs: dict) -> dict:
    annotation = {
        "gene": gene,
        "node_id": original_node,
        "display": node_display_name(attrs),
        "node_type": "leaf" if attrs.get("is_leaf") else "internal",
        "taxon_clade": attrs.get("mammal_clade", ""),
        "species": attrs.get("species", ""),
        "species_label": attrs.get("species_label", ""),
    }
    optional_keys = [
        "orthology_type",
        "taxonomy_level",
        "ensembl_gene_id",
        "protein_id",
        "taxon_id",
        "percent_identity",
        "percent_positive",
        "alignment_coverage",
        "variant_id",
        "assigned_af",
        "visual_leaf_support",
        "node_radius",
        "supporting_variant_count",
        "supporting_mixed_candidates",
    ]
    for key in optional_keys:
        if attrs.get(key) not in (None, ""):
            annotation[key] = attrs[key]
    if attrs.get("is_human_population_summary"):
        annotation["human_display_policy"] = "neutral human subtree leaf; population labels intentionally hidden"
    return annotation


def base_node_color(attrs: dict) -> tuple[int, int, int, int]:
    if attrs.get("is_interlayer_clade_port"):
        return hex_to_rgba(
            pst.MAMMAL_CLADE_COLORS.get(str(attrs.get("interlayer_group")), "#CFE6FF"),
            92,
        )
    if attrs.get("is_interlayer_clade_control"):
        return INTERLAYER_FLOW_NODE_COLOR
    if attrs.get("is_interlayer_taxon_control"):
        return INTERLAYER_SELECTED_TAXON_NODE_COLOR
    if attrs.get("is_leaf") and attrs.get("is_mammal_ortholog"):
        clade = attrs.get("mammal_clade", "Species ortholog context")
        return hex_to_rgba(pst.MAMMAL_CLADE_COLORS.get(clade, "#E5E7EB"), 230)
    if attrs.get("is_human_sequence_anchor") or attrs.get("is_human_variant_context") or attrs.get("is_human_population_summary"):
        return hex_to_rgba(pst.HUMAN_SUBTREE_COLOR, 220 if attrs.get("is_leaf") else 120)
    if attrs.get("is_leaf"):
        return hex_to_rgba(pst.HUMAN_SUBTREE_COLOR, 220)
    return (202, 220, 250, 86)


def faded_node_color(attrs: dict) -> tuple[int, int, int, int]:
    if attrs.get("is_leaf"):
        return (118, 132, 152, 56)
    return (88, 100, 118, 30)


def edge_color_for_child(attrs: dict, faded: bool = False) -> tuple[int, int, int, int]:
    if faded:
        return FADED_LINK_COLOR
    if attrs.get("is_mammal_ortholog"):
        clade = attrs.get("mammal_clade", "Species ortholog context")
        return hex_to_rgba(pst.MAMMAL_CLADE_COLORS.get(clade, "#E5E7EB"), 196)
    if attrs.get("is_human_sequence_anchor") or attrs.get("is_human_variant_context") or attrs.get("is_human_population_summary"):
        return hex_to_rgba(pst.HUMAN_SUBTREE_COLOR, 188)
    return AMBIGUOUS_LINK_COLOR


def node_radius_for_support(attrs: dict, support: int, max_support: int) -> float:
    if attrs.get("is_interlayer_clade_control") or attrs.get("is_interlayer_taxon_control"):
        return 0.34
    if attrs.get("is_interlayer_clade_port"):
        return 0.58
    if attrs.get("is_leaf"):
        return 1.0
    if max_support <= 1:
        return 0.72
    normalized = math.log1p(max(1, int(support))) / math.log1p(max_support)
    return round(0.58 + (0.92 * (normalized**0.72)), 4)


def leaf_path_group(attrs: dict) -> str:
    if (
        attrs.get("is_human_sequence_anchor")
        or attrs.get("is_human_variant_context")
        or attrs.get("is_human_population_summary")
        or (attrs.get("is_leaf") and not attrs.get("is_mammal_ortholog"))
    ):
        return "Human subtree"
    if attrs.get("is_mammal_ortholog"):
        return str(attrs.get("mammal_clade") or "Species ortholog context")
    return "Species ortholog context"


def subtree_leaf_group_counts(graph: nx.DiGraph) -> dict[str, Counter]:
    cache: dict[str, Counter] = {}

    def count_for(node: str) -> Counter:
        if node in cache:
            return cache[node]
        attrs = graph.nodes[node]
        children = list(graph.successors(node))
        if attrs.get("is_leaf") or not children:
            cache[node] = Counter({leaf_path_group(attrs): 1})
            return cache[node]
        counts: Counter = Counter()
        for child in children:
            counts.update(count_for(child))
        cache[node] = counts
        return counts

    for node in graph.nodes:
        count_for(node)
    return cache


def path_group_from_counts(counts: Counter) -> tuple[str, float, int]:
    total = int(sum(counts.values()))
    if total <= 0:
        return "Backbone", 0.0, 0
    dominant_group, dominant_count = counts.most_common(1)[0]
    dominance = float(dominant_count / total)
    if dominant_group != "Species ortholog context" and dominance >= PATH_GROUP_DOMINANCE_THRESHOLD:
        return dominant_group, dominance, total
    if len(counts) == 1 and dominant_group != "Species ortholog context":
        return dominant_group, 1.0, total
    return "Backbone", dominance, total


def edge_path_group_assignments(graph: nx.DiGraph) -> dict[tuple[str, str], dict]:
    counts_by_node = subtree_leaf_group_counts(graph)
    assignments = {}
    for parent, child in graph.edges:
        counts = counts_by_node[str(child)]
        group, dominance, total = path_group_from_counts(counts)
        assignments[(str(parent), str(child))] = {
            "path_group": group,
            "path_group_dominance": round(dominance, 4),
            "path_group_leaf_count": total,
            "path_group_counts": dict(counts),
        }
    return assignments


def edge_color_for_path_group(
    path_group: str,
    dominance: float,
    faded: bool = False,
    focused: bool = False,
) -> tuple[int, int, int, int]:
    if faded:
        return FADED_LINK_COLOR
    if path_group == "Backbone":
        return BACKBONE_LINK_COLOR
    if path_group == "Human subtree":
        return hex_to_rgba(pst.HUMAN_SUBTREE_COLOR, 142 if focused else 112)
    if path_group in pst.MAMMAL_CLADE_COLORS:
        alpha = int(86 + 44 * clamp(dominance, PATH_GROUP_DOMINANCE_THRESHOLD, 1.0))
        if focused:
            alpha = min(152, alpha + 22)
        return hex_to_rgba(pst.MAMMAL_CLADE_COLORS[path_group], alpha)
    return edge_color_for_child({}, faded=faded)


def interlayer_taxon_key(attrs: dict) -> str | None:
    if attrs.get("is_mammal_ortholog") and attrs.get("species"):
        return f"species::{attrs['species']}"
    if attrs.get("is_human_sequence_anchor"):
        return f"species::{pst.HUMAN_SEQUENCE_LEAF}"
    return None


def interlayer_taxon_index(graph: nx.DiGraph) -> dict[str, str]:
    index = {}
    for node, attrs in graph.nodes(data=True):
        key = interlayer_taxon_key(attrs)
        if key:
            index[key] = str(node)
    return index


def interlayer_taxon_label(attrs: dict) -> str:
    if attrs.get("is_human_sequence_anchor"):
        return "Human"
    return str(attrs.get("species_label") or attrs.get("label") or attrs.get("species") or "species")


def interlayer_taxon_group(attrs: dict) -> str:
    if attrs.get("is_human_sequence_anchor"):
        return "Human"
    return str(attrs.get("mammal_clade") or "Species ortholog context")


def edge_color_for_interlayer_group(
    group: str,
    faded: bool = False,
    focused: bool = False,
) -> tuple[int, int, int, int]:
    if faded:
        return FADED_LINK_COLOR
    if group == "Human":
        color = INTERLAYER_HUMAN_CONNECTOR_COLOR
    else:
        color = hex_to_rgba(
            pst.MAMMAL_CLADE_COLORS.get(group, "#CFE6FF"),
            124,
        )
    if focused:
        return (color[0], color[1], color[2], min(154, color[3] + 20))
    return color


def edge_color_for_interlayer_edge(
    attrs: dict,
    faded: bool = False,
    focused: bool = False,
) -> tuple[int, int, int, int]:
    if faded:
        return FADED_LINK_COLOR
    group = str(attrs.get("interlayer_group") or "Species")
    if group == "Human":
        base = INTERLAYER_HUMAN_CONNECTOR_COLOR
    else:
        base = hex_to_rgba(
            pst.MAMMAL_CLADE_COLORS.get(group, "#CFE6FF"),
            INTERLAYER_CLADE_FLOW_ALPHA,
        )
    if attrs.get("is_interlayer_taxon_flow_edge"):
        if focused:
            return (base[0], base[1], base[2], INTERLAYER_SELECTED_TAXON_FOCUS_ALPHA)
        return (48, 56, 72, INTERLAYER_SELECTED_TAXON_FLOW_ALPHA)
    if attrs.get("is_interlayer_ancestor_context_edge"):
        alpha = 58 if not focused else 92
        return (base[0], base[1], base[2], alpha)
    if focused:
        return (base[0], base[1], base[2], min(154, base[3] + 20))
    return base


def scene_node_group_name(attrs: dict, faded: bool) -> str:
    original = attrs.get("original_attrs", attrs)
    if original.get("is_interlayer_clade_port"):
        return f"{original.get('interlayer_group', 'Clade')} clade flow ports"
    if original.get("is_interlayer_clade_control"):
        return "Clade flow corridor controls"
    if original.get("is_interlayer_taxon_control"):
        return "Selected taxon flow controls"
    if faded:
        return "Faded non-focused leaves" if original.get("is_leaf") else "Faded non-focused backbone"
    if original.get("is_leaf") and original.get("is_mammal_ortholog"):
        return str(original.get("mammal_clade") or "Species ortholog context")
    if (
        original.get("is_human_sequence_anchor")
        or original.get("is_human_variant_context")
        or original.get("is_human_population_summary")
        or (original.get("is_leaf") and not original.get("is_mammal_ortholog"))
    ):
        return "Human subtree leaves" if original.get("is_leaf") else "Human subtree connector"
    return "Backbone / internal splits"


def named_group_selections(scenes: list[nx.Graph], node_order: list[str]) -> list[dict]:
    node_to_id = {node: index for index, node in enumerate(node_order)}
    group_order = {
        group: index
        for index, group in enumerate(
            [
                *pst.MAMMAL_CLADE_ORDER,
                "Human subtree leaves",
                "Human subtree connector",
                *[f"{group} clade flow ports" for group in pst.MAMMAL_CLADE_ORDER],
                "Human clade flow ports",
                "Clade flow corridor controls",
                "Selected taxon flow controls",
                "Backbone / internal splits",
                "Faded non-focused leaves",
                "Faded non-focused backbone",
            ]
        )
    }
    selections = []
    for scene in scenes:
        focus_gene = scene.graph.get("focus_gene")
        groups: dict[tuple[str, tuple[int, int, int, int]], list[str]] = {}
        for node in node_order:
            attrs = scene.nodes[node]
            if attrs.get("is_interlayer_flow_node"):
                faded = bool(
                    focus_gene not in (None, "all")
                    and focus_gene not in attrs.get("flow_genes", [])
                )
            else:
                faded = bool(focus_gene not in (None, "all") and attrs["gene"] != focus_gene)
            group_name = scene_node_group_name(attrs, faded)
            color = tuple(int(value) for value in attrs["nodecolor"])
            groups.setdefault((group_name, color), []).append(str(node_to_id[node]))
        for (group_name, color), node_ids in sorted(
            groups.items(),
            key=lambda item: (group_order.get(item[0][0], 10_000), item[0][0], item[0][1]),
        ):
            selections.append(
                {
                    "name": group_name,
                    "nodes": node_ids,
                    "layoutname": scene.graph["layoutname"],
                    "labelcolor": list(color),
                }
            )
    return selections


def scene_link_group_name(attrs: dict, faded: bool) -> str:
    if faded:
        if attrs.get("is_interlayer_connection"):
            return "Faded non-focused interlayer flows"
        return "Faded non-focused gene paths"
    if attrs.get("is_interlayer_connection"):
        group = str(attrs.get("interlayer_group") or "Species")
        if attrs.get("is_interlayer_taxon_flow_edge"):
            if group == "Human":
                return "Human selected same-taxon history tracks"
            return f"{group} selected same-taxon history tracks"
        if group == "Human":
            return "Human clade-level flow corridor"
        return f"{group} clade-level flow corridor"
    path_group = str(attrs.get("path_group") or "Backbone")
    if path_group == "Backbone":
        return "Backbone / mixed ancestral paths"
    if path_group == "Human subtree":
        return "Human subtree paths"
    return f"{path_group} paths"


def named_link_group_selections(scenes: list[nx.Graph], edge_order: list[tuple[str, str]]) -> list[dict]:
    link_order = {edge: index for index, edge in enumerate(edge_order)}
    link_order.update({(target, source): index for (source, target), index in link_order.items()})
    group_order = {
        group: index
        for index, group in enumerate(
            [
                *[f"{group} paths" for group in pst.MAMMAL_CLADE_ORDER],
                "Human subtree paths",
                "Backbone / mixed ancestral paths",
                *[f"{group} clade-level flow corridor" for group in pst.MAMMAL_CLADE_ORDER],
                "Human clade-level flow corridor",
                *[f"{group} selected same-taxon history tracks" for group in pst.MAMMAL_CLADE_ORDER],
                "Human selected same-taxon history tracks",
                "Faded non-focused gene paths",
                "Faded non-focused interlayer flows",
            ]
        )
    }
    selections = []
    for scene in scenes:
        focus_gene = scene.graph.get("focus_gene")
        groups: dict[tuple[str, tuple[int, int, int, int]], list[str]] = {}
        for source, target in edge_order:
            attrs = scene.edges[source, target]
            if attrs.get("is_interlayer_connection"):
                flow_genes = attrs.get("flow_genes", [])
                touches_focus = bool(
                    focus_gene not in (None, "all")
                    and flow_genes
                    and focus_gene in flow_genes
                )
                faded = bool(focus_gene not in (None, "all") and flow_genes and not touches_focus)
            else:
                faded = bool(focus_gene not in (None, "all") and attrs["gene"] != focus_gene)
            group_name = scene_link_group_name(attrs, faded)
            color = tuple(int(value) for value in attrs["linkcolor"])
            groups.setdefault((group_name, color), []).append(str(link_order[(source, target)]))
        for (group_name, color), link_ids in sorted(
            groups.items(),
            key=lambda item: (group_order.get(item[0][0], 10_000), item[0][0], item[0][1]),
        ):
            selections.append(
                {
                    "name": group_name,
                    "links": link_ids,
                    "layoutname": scene.graph["layoutname"],
                    "labelcolor": list(color),
                }
            )
    return selections


def build_inferred_gene_trees() -> tuple[dict[str, nx.DiGraph], dict[str, dict[str, np.ndarray]]]:
    source_graphs = pst.load_graphs()
    context = pst.load_mammal_context()
    human_subtree_source_gene = "ACTB"
    gene_source_graphs = {}
    for gene in pst.GENE_ORDER:
        if gene in source_graphs:
            gene_source_graphs[gene] = source_graphs[gene]
            continue
        proxy = source_graphs[human_subtree_source_gene].copy()
        proxy.graph.update(source_graphs[human_subtree_source_gene].graph)
        proxy.graph["human_subtree_proxy_source_gene"] = human_subtree_source_gene
        proxy.graph["human_subtree_proxy_note"] = (
            f"{gene} ortholog layer uses the existing {human_subtree_source_gene} neutral human-subtree "
            "variant display because no gene-specific human variant tree is cached for this gene"
        )
        gene_source_graphs[gene] = proxy
    graphs = {
        gene: pst.compose_inferred_taxon_aligned_tree(gene_source_graphs[gene], gene, context)
        for gene in pst.GENE_ORDER
    }
    layouts = {
        gene: pst.compute_unrooted_3d_layout(graphs[gene], gene, geometry="hyperbolic")
        for gene in pst.GENE_ORDER
    }
    return graphs, layouts


def normalize_positions(layouts: dict[str, dict[str, np.ndarray]]) -> dict[tuple[str, str], list[float]]:
    keyed_points = []
    for gene in pst.GENE_ORDER:
        for node, point in layouts[gene].items():
            keyed_points.append(((gene, str(node)), np.array(point, dtype=float)))
    matrix = np.array([point for _, point in keyed_points], dtype=float)
    low = matrix.min(axis=0)
    high = matrix.max(axis=0)
    span = np.maximum(high - low, 1e-9)

    normalized = {}
    for key, point in keyed_points:
        # Keep a margin so the DataDiVR coordinate texture never encodes exact 0 or 1.
        norm = 0.045 + ((point - low) / span) * 0.91
        normalized[key] = [round(float(clamp(value, 0.001, 0.999)), 7) for value in norm]
    return normalized


def combined_node_id(gene: str, node: str) -> str:
    return f"{gene}|{node}"


def add_interlayer_species_connections(
    base: nx.Graph,
    gene_graphs: dict[str, nx.DiGraph],
    layouts: dict[str, dict[str, np.ndarray]],
) -> dict:
    indices = {gene: interlayer_taxon_index(gene_graphs[gene]) for gene in pst.GENE_ORDER}
    shared_keys = sorted(set.intersection(*(set(index) for index in indices.values())))

    taxon_records = []
    grouped_taxa: dict[str, list[dict]] = {}
    for key in shared_keys:
        reference_gene = pst.GENE_ORDER[0]
        reference_node = indices[reference_gene][key]
        reference_attrs = gene_graphs[reference_gene].nodes[reference_node]
        label = interlayer_taxon_label(reference_attrs)
        group = interlayer_taxon_group(reference_attrs)
        record = {"key": key, "label": label, "group": group}
        taxon_records.append(record)
        grouped_taxa.setdefault(group, []).append(record)

    edge_records = []
    node_records = []
    clade_flow_records = []
    selected_taxon_flow_records = []
    path_records = []
    group_order = [group for group in [*pst.MAMMAL_CLADE_ORDER, "Human"] if group in grouped_taxa]
    adjacent_pairs = list(zip(pst.GENE_ORDER[:-1], pst.GENE_ORDER[1:]))
    parent_maps: dict[str, dict[str, str | None]] = {}

    def parent_map_for_gene(gene: str) -> dict[str, str | None]:
        if gene in parent_maps:
            return parent_maps[gene]
        graph = gene_graphs[gene]
        root = pst.root_node(graph)
        parent: dict[str, str | None] = {root: None}
        stack = [root]
        while stack:
            node = stack.pop()
            for child in graph.successors(node):
                child = str(child)
                if child in parent:
                    continue
                parent[child] = str(node)
                stack.append(child)
        parent_maps[gene] = parent
        return parent

    def ancestor_at_distance(gene: str, node: str, distance: int) -> str | None:
        parent = parent_map_for_gene(gene)
        current: str | None = str(node)
        for _ in range(distance):
            if current is None:
                return None
            current = parent.get(current)
        return current

    def rounded_position(point: np.ndarray | list[float]) -> list[float]:
        array = np.array(point, dtype=float)
        return [round(float(clamp(value, 0.001, 0.999)), 7) for value in array]

    def centroid(points: list[list[float]]) -> list[float]:
        return rounded_position(np.array(points, dtype=float).mean(axis=0))

    def endpoint_position(gene: str, key: str) -> list[float]:
        return base.nodes[combined_node_id(gene, indices[gene][key])]["pos"]

    def flow_focus_genes_for_port(gene: str) -> list[str]:
        index = pst.GENE_ORDER.index(gene)
        focus_genes = {gene}
        if index > 0:
            focus_genes.add(pst.GENE_ORDER[index - 1])
        if index < len(pst.GENE_ORDER) - 1:
            focus_genes.add(pst.GENE_ORDER[index + 1])
        return [item for item in pst.GENE_ORDER if item in focus_genes]

    def flow_control_positions(
        start: list[float],
        end: list[float],
        lane_index: int,
        lane_total: int,
        offset_scale: float,
    ) -> list[float]:
        start_array = np.array(start, dtype=float)
        end_array = np.array(end, dtype=float)
        delta = end_array - start_array
        angle = math.tau * ((lane_index + 0.5) / max(lane_total, 1))
        lane_fraction = (lane_index % 5) / 4.0
        radius = offset_scale * (0.52 + 0.46 * lane_fraction)
        lane_xy = np.array(
            [
                clamp(0.5 + (math.cos(angle) * radius), 0.44, 0.56),
                clamp(0.5 + (math.sin(angle) * radius), 0.44, 0.56),
            ],
            dtype=float,
        )
        control = start_array + (delta * 0.50)
        control[:2] = lane_xy
        return rounded_position(control)

    def clade_corridor_control_positions(
        start: list[float],
        end: list[float],
        group_index: int,
        pair_index: int,
        group_total: int,
        context_level: int = 0,
    ) -> list[float]:
        start_array = np.array(start, dtype=float)
        end_array = np.array(end, dtype=float)
        delta = end_array - start_array
        angle = (
            (-math.pi / 2.0)
            + (math.tau * group_index / max(group_total, 1))
            + (pair_index * 0.10)
            + (context_level * 0.055)
        )
        radius = 0.18 + (0.035 * (group_index % 4)) + (0.040 * context_level)
        lane_anchor = np.array(
            [
                clamp(0.5 + (math.cos(angle) * radius), 0.16, 0.84),
                clamp(0.5 + (math.sin(angle) * radius), 0.16, 0.84),
            ],
            dtype=float,
        )
        tangent = np.array([-math.sin(angle), math.cos(angle)], dtype=float)
        pair_splay = (0.018 + 0.006 * context_level) if pair_index == 0 else (-0.018 - 0.006 * context_level)
        control = start_array + (delta * 0.50)
        control[:2] = (control[:2] * 0.25) + (lane_anchor * 0.75) + (tangent * pair_splay)
        z_bow = (0.025 + 0.010 * context_level) * math.sin(angle + pair_index)
        control[2] = clamp(control[2] + z_bow, 0.045, 0.955)
        return rounded_position(control)

    def flow_node_attrs(
        node_id: str,
        label: str,
        pos: list[float],
        group: str,
        node_kind: str,
        flow_genes: list[str],
        gene_pair: str,
        taxon_count: int = 0,
        taxon_label: str | None = None,
    ) -> dict:
        radius = node_radius_for_support(
            {f"is_interlayer_{node_kind}": True},
            taxon_count,
            max(1, taxon_count),
        )
        original_attrs = {
            "label": label,
            "is_leaf": False,
            "taxon_count": taxon_count,
            "interlayer_group": group,
            "is_interlayer_flow_node": True,
            f"is_interlayer_{node_kind}": True,
            "gene_pair": gene_pair,
            "flow_genes": flow_genes,
            "visual_leaf_support": taxon_count,
            "node_radius": radius,
        }
        return {
            "name": label,
            "pos": pos,
            "gene": "interlayer_flow",
            "flow_genes": flow_genes,
            "gene_pair": gene_pair,
            "original_node": node_id,
            "original_attrs": original_attrs,
            "annotation": {
                "display": label,
                "node_type": f"interlayer {node_kind.replace('_', ' ')}",
                "interlayer_group": group,
                "gene_pair": gene_pair,
                "connected_taxa": taxon_count,
                "taxon": taxon_label or "",
                "meaning": (
                    "helper node used to draw a segmented DataDiVR flow path; "
                    "DataDiVR only supports straight, fixed-width links"
                ),
            },
            "cluster": "Interlayer gene-history flows",
            "interlayer_group": group,
            "is_interlayer_flow_node": True,
            f"is_interlayer_{node_kind}": True,
            "is_leaf": False,
            "visual_leaf_support": taxon_count,
            "node_radius": radius,
        }

    port_nodes: dict[tuple[str, str], str] = {}

    def clade_port_node(gene: str, group: str) -> str:
        key = (gene, group)
        if key in port_nodes:
            return port_nodes[key]
        endpoint_positions = [endpoint_position(gene, record["key"]) for record in grouped_taxa[group]]
        port_pos = centroid(endpoint_positions)
        node_id = f"__flow_port__{safe_name(group)}__{gene}"
        base.add_node(
            node_id,
            **flow_node_attrs(
                node_id=node_id,
                label=f"{gene} {group} clade-flow port",
                pos=port_pos,
                group=group,
                node_kind="clade_port",
                flow_genes=flow_focus_genes_for_port(gene),
                gene_pair=f"{gene} clade port",
                taxon_count=len(grouped_taxa[group]),
            ),
        )
        node_records.append({"node": node_id, "group": group, "gene": gene, "kind": "clade_port"})
        port_nodes[key] = node_id
        return node_id

    def add_flow_edge(source: str, target: str, attrs: dict, edge_type: str) -> None:
        base.add_edge(source, target, **attrs)
        edge_records.append(
            {
                "source": source,
                "target": target,
                "group": attrs.get("interlayer_group"),
                "gene_pair": attrs.get("gene_pair"),
                "edge_type": edge_type,
                "flow_kind": attrs.get("flow_kind"),
            }
        )

    def mrca_node_for_group(gene: str, group: str) -> str | None:
        graph = gene_graphs[gene]
        selected_nodes = [
            indices[gene][record["key"]]
            for record in grouped_taxa[group]
            if record["key"] in indices[gene]
        ]
        if not selected_nodes:
            return None
        if len(selected_nodes) == 1:
            return str(selected_nodes[0])

        root = pst.root_node(graph)
        parent: dict[str, str | None] = {root: None}
        stack = [root]
        while stack:
            node = stack.pop()
            for child in graph.successors(node):
                child = str(child)
                if child in parent:
                    continue
                parent[child] = str(node)
                stack.append(child)

        def path_from_root(node: str) -> list[str]:
            path = []
            current: str | None = str(node)
            while current is not None and current in parent:
                path.append(current)
                current = parent[current]
            return list(reversed(path))

        paths = [path_from_root(str(node)) for node in selected_nodes]
        if not all(paths):
            return str(selected_nodes[0])
        mrca = paths[0][0]
        for nodes_at_depth in zip(*paths):
            if len(set(nodes_at_depth)) != 1:
                break
            mrca = nodes_at_depth[0]
        return str(mrca)

    def add_segmented_clade_flow(
        source: str,
        target: str,
        group: str,
        group_index: int,
        pair_index: int,
        source_gene: str,
        target_gene: str,
        flow_kind: str,
        label_suffix: str,
        path_kind: str = "clade_flow",
        context_level: int = 0,
    ) -> None:
        gene_pair = f"{source_gene}->{target_gene}"
        start = base.nodes[source]["pos"]
        end = base.nodes[target]["pos"]
        control_pos = clade_corridor_control_positions(
            start,
            end,
            group_index,
            pair_index,
            len(group_order),
            context_level=context_level,
        )
        control = f"__flow_control__{safe_name(group)}__{source_gene}_{target_gene}__{flow_kind}"
        flow_genes = [source_gene, target_gene]
        base.add_node(
            control,
            **flow_node_attrs(
                node_id=control,
                label=f"{gene_pair} {group} {label_suffix} control",
                pos=control_pos,
                group=group,
                node_kind="clade_control",
                flow_genes=flow_genes,
                gene_pair=gene_pair,
                taxon_count=len(grouped_taxa[group]),
            ),
        )
        node_records.append(
            {
                "node": control,
                "group": group,
                "gene_pair": gene_pair,
                "kind": "clade_control",
                "flow_kind": flow_kind,
            }
        )
        flow_attrs = {
            "gene": "interlayer",
            "flow_genes": flow_genes,
            "gene_pair": gene_pair,
            "is_interlayer_connection": True,
            "is_interlayer_clade_flow_edge": True,
            "is_interlayer_ancestor_context_edge": path_kind == "ancestor_context_flow",
            "flow_kind": flow_kind,
            "context_level": context_level,
            "interlayer_group": group,
            "interlayer_label": f"{group} {label_suffix}",
            "path_group": "Interlayer clade-level flow corridor",
            "path_group_dominance": 1.0,
            "path_group_leaf_count": len(grouped_taxa[group]),
            "path_group_counts": {group: len(grouped_taxa[group])},
        }
        add_flow_edge(source, control, copy.deepcopy(flow_attrs), "clade_flow_segment")
        add_flow_edge(control, target, copy.deepcopy(flow_attrs), "clade_flow_segment")
        clade_flow_records.append(
            {
                "gene_pair": gene_pair,
                "group": group,
                "taxa": len(grouped_taxa[group]),
                "source": source,
                "target": target,
                "flow_kind": flow_kind,
                "path_kind": path_kind,
                "context_level": context_level,
                "path_nodes": [source, control, target],
            }
        )
        path_records.append(
            {
                "name": f"{gene_pair} {group} {flow_kind}",
                "kind": path_kind,
                "gene_pair": gene_pair,
                "group": group,
                "flow_kind": flow_kind,
                "context_level": context_level,
                "taxa": len(grouped_taxa[group]),
                "path_nodes": [source, control, target],
            }
        )

    middle_gene = pst.GENE_ORDER[1]
    first_gene = pst.GENE_ORDER[0]
    third_gene = pst.GENE_ORDER[2]
    for group_index, group in enumerate(group_order):
        first_mrca = mrca_node_for_group(first_gene, group)
        middle_mrca = mrca_node_for_group(middle_gene, group)
        third_mrca = mrca_node_for_group(third_gene, group)
        if first_mrca is None or middle_mrca is None or third_mrca is None:
            continue
        first_mrca_key = combined_node_id(first_gene, first_mrca)
        middle_mrca_key = combined_node_id(middle_gene, middle_mrca)
        third_mrca_key = combined_node_id(third_gene, third_mrca)
        if first_mrca_key not in base or middle_mrca_key not in base or third_mrca_key not in base:
            continue
        add_segmented_clade_flow(
            source=first_mrca_key,
            target=middle_mrca_key,
            group=group,
            group_index=group_index,
            pair_index=0,
            source_gene=first_gene,
            target_gene=middle_gene,
            flow_kind="converge_to_middle_mrca",
            label_suffix=f"MRCA-to-MRCA convergence into {middle_gene}",
        )
        add_segmented_clade_flow(
            source=middle_mrca_key,
            target=third_mrca_key,
            group=group,
            group_index=group_index,
            pair_index=1,
            source_gene=middle_gene,
            target_gene=third_gene,
            flow_kind="diverge_from_middle_mrca",
            label_suffix=f"MRCA-to-MRCA divergence from {middle_gene}",
        )
        for context_level in ANCESTOR_CONTEXT_LEVELS:
            first_ancestor = ancestor_at_distance(first_gene, first_mrca, context_level)
            middle_ancestor = ancestor_at_distance(middle_gene, middle_mrca, context_level)
            third_ancestor = ancestor_at_distance(third_gene, third_mrca, context_level)
            if first_ancestor is None or middle_ancestor is None or third_ancestor is None:
                continue
            first_ancestor_key = combined_node_id(first_gene, first_ancestor)
            middle_ancestor_key = combined_node_id(middle_gene, middle_ancestor)
            third_ancestor_key = combined_node_id(third_gene, third_ancestor)
            if first_ancestor_key not in base or middle_ancestor_key not in base or third_ancestor_key not in base:
                continue
            add_segmented_clade_flow(
                source=first_ancestor_key,
                target=middle_ancestor_key,
                group=group,
                group_index=group_index,
                pair_index=0,
                source_gene=first_gene,
                target_gene=middle_gene,
                flow_kind=f"ancestor_context_level_{context_level}",
                label_suffix=f"ancestor context level {context_level} around {middle_gene}",
                path_kind="ancestor_context_flow",
                context_level=context_level,
            )
            add_segmented_clade_flow(
                source=middle_ancestor_key,
                target=third_ancestor_key,
                group=group,
                group_index=group_index,
                pair_index=1,
                source_gene=middle_gene,
                target_gene=third_gene,
                flow_kind=f"ancestor_context_level_{context_level}",
                label_suffix=f"ancestor context level {context_level} around {middle_gene}",
                path_kind="ancestor_context_flow",
                context_level=context_level,
            )

    history_records = pst.history_flow_records(gene_graphs, layouts)
    taxon_lane_total = max(1, len(history_records))
    for flow_index, record in enumerate(history_records):
        source_gene = record["source_gene"]
        target_gene = record["target_gene"]
        source_node = str(record["source_node"])
        target_node = str(record["target_node"])
        source_key = combined_node_id(source_gene, source_node)
        target_key = combined_node_id(target_gene, target_node)
        if source_key not in base or target_key not in base:
            continue
        source_attrs = record["source_attrs"]
        group = interlayer_taxon_group(source_attrs)
        label = interlayer_taxon_label(source_attrs)
        gene_pair = f"{source_gene}->{target_gene}"
        control_pos = flow_control_positions(
            base.nodes[source_key]["pos"],
            base.nodes[target_key]["pos"],
            flow_index,
            taxon_lane_total,
            offset_scale=0.090,
        )
        safe_key = safe_name(record["key"])
        control = f"__taxon_flow_control__{safe_key}__{source_gene}_{target_gene}"
        flow_genes = [source_gene, target_gene]
        base.add_node(
            control,
            **flow_node_attrs(
                node_id=control,
                label=f"{gene_pair} {label} history-track control",
                pos=control_pos,
                group=group,
                node_kind="taxon_control",
                flow_genes=flow_genes,
                gene_pair=gene_pair,
                taxon_count=1,
                taxon_label=label,
            ),
        )
        node_records.append(
            {
                "node": control,
                "group": group,
                "gene_pair": gene_pair,
                "kind": "taxon_control",
                "taxon": label,
            }
        )
        taxon_flow_attrs = {
            "gene": "interlayer",
            "flow_genes": flow_genes,
            "gene_pair": gene_pair,
            "is_interlayer_connection": True,
            "is_interlayer_taxon_flow_edge": True,
            "flow_kind": "selected_same_taxon_history_track",
            "interlayer_key": record["key"],
            "interlayer_group": group,
            "interlayer_label": label,
            "path_group": "Interlayer selected same-taxon history track",
            "path_group_dominance": 1.0,
            "path_group_leaf_count": 1,
            "path_group_counts": {group: 1},
            "flow_score": round(float(record["score"]), 6),
            "screen_shift": round(float(record["screen_shift"]), 6),
            "human_distance_delta": round(float(record["human_distance_delta"]), 6),
            "identity_delta": round(float(record["identity_delta"]), 6),
            "must_show": bool(record["must_show"]),
        }
        add_flow_edge(source_key, control, copy.deepcopy(taxon_flow_attrs), "taxon_flow_segment")
        add_flow_edge(control, target_key, copy.deepcopy(taxon_flow_attrs), "taxon_flow_segment")
        selected_taxon_flow_records.append(
            {
                "gene_pair": gene_pair,
                "key": record["key"],
                "label": label,
                "group": group,
                "score": round(float(record["score"]), 6),
                "must_show": bool(record["must_show"]),
                "path_nodes": [source_key, control, target_key],
            }
        )
        path_records.append(
            {
                "name": f"{gene_pair} {label}",
                "kind": "selected_taxon_flow",
                "gene_pair": gene_pair,
                "group": group,
                "taxon_key": record["key"],
                "label": label,
                "score": round(float(record["score"]), 6),
                "must_show": bool(record["must_show"]),
                "path_nodes": [source_key, control, target_key],
            }
        )

    return {
        "shared_taxa": taxon_records,
        "flow_nodes": node_records,
        "edges": edge_records,
        "clade_flow_records": clade_flow_records,
        "selected_taxon_flow_records": selected_taxon_flow_records,
        "path_records": path_records,
    }


def build_combined_graphs() -> tuple[list[nx.Graph], dict]:
    gene_graphs, layouts = build_inferred_gene_trees()
    positions = normalize_positions(layouts)
    path_groups_by_gene = {
        gene: edge_path_group_assignments(gene_graphs[gene])
        for gene in pst.GENE_ORDER
    }
    node_support_by_gene = {
        gene: pst.descendant_leaf_supports(gene_graphs[gene])
        for gene in pst.GENE_ORDER
    }
    max_node_support_by_gene = {
        gene: max(node_support_by_gene[gene].values(), default=1)
        for gene in pst.GENE_ORDER
    }
    base = nx.Graph()

    for gene in pst.GENE_ORDER:
        graph = gene_graphs[gene]
        for node, attrs in graph.nodes(data=True):
            node_key = combined_node_id(gene, str(node))
            support = node_support_by_gene[gene].get(str(node), 1)
            node_radius = node_radius_for_support(
                attrs,
                support,
                max_node_support_by_gene[gene],
            )
            base.add_node(
                node_key,
                name=node_display_name(attrs),
                pos=positions[(gene, str(node))],
                gene=gene,
                original_node=str(node),
                original_attrs=copy.deepcopy(attrs),
                annotation={
                    **node_annotation(gene, str(node), attrs),
                    "visual_leaf_support": support,
                    "node_radius": node_radius,
                },
                cluster=str(attrs.get("mammal_clade", "")) or None,
                visual_leaf_support=support,
                node_radius=node_radius,
            )
        for parent, child, edge_attrs in graph.edges(data=True):
            parent_key = combined_node_id(gene, str(parent))
            child_key = combined_node_id(gene, str(child))
            child_attrs = graph.nodes[child]
            path_group = path_groups_by_gene[gene][(str(parent), str(child))]
            base.add_edge(
                parent_key,
                child_key,
                gene=gene,
                original_edge=[str(parent), str(child)],
                weight=float(edge_attrs.get("weight", 1.0)),
                child_attrs=copy.deepcopy(child_attrs),
                **path_group,
            )

    interlayer_records = add_interlayer_species_connections(base, gene_graphs, layouts)

    base.graph.update(
        {
            "projectname": PROJECT_NAME,
            "info": INFO,
            "geometry": "stacked_3d_hyperbolic_unrooted_tree_layouts",
            "source": "plotly_stacked_trees.py inferred NetworkX trees",
            "path_records": interlayer_records["path_records"],
        }
    )

    scenes = []
    for scene_name in SCENE_NAMES:
        focus_gene = scene_name.split("_")[1] if scene_name.endswith("_focus") else None
        scene = nx.Graph()
        scene.graph.update(base.graph)
        scene.graph["layoutname"] = scene_name
        scene.graph["focus_gene"] = focus_gene or "all"
        for node, attrs in base.nodes(data=True):
            copied = copy.deepcopy(attrs)
            if copied.get("is_interlayer_flow_node"):
                faded = bool(focus_gene and focus_gene not in copied.get("flow_genes", []))
                copied["nodecolor"] = faded_node_color(copied["original_attrs"]) if faded else base_node_color(copied)
            else:
                faded = bool(focus_gene and copied["gene"] != focus_gene)
                copied["nodecolor"] = faded_node_color(copied["original_attrs"]) if faded else base_node_color(copied["original_attrs"])
            scene.add_node(node, **copied)
        for source, target, attrs in base.edges(data=True):
            if attrs.get("is_interlayer_connection"):
                flow_genes = attrs.get("flow_genes", [])
                touches_focus = bool(focus_gene and focus_gene in flow_genes)
                faded = bool(focus_gene and flow_genes and not touches_focus)
                focused = bool(focus_gene and touches_focus)
            else:
                faded = bool(focus_gene and attrs["gene"] != focus_gene)
                focused = bool(focus_gene and attrs["gene"] == focus_gene)
            scene.add_edge(
                source,
                target,
                **{
                    key: copy.deepcopy(value)
                    for key, value in attrs.items()
                    if key != "child_attrs"
                },
            )
            if attrs.get("is_interlayer_connection"):
                scene.edges[source, target]["linkcolor"] = edge_color_for_interlayer_edge(
                    attrs,
                    faded=faded,
                    focused=focused,
                )
            else:
                scene.edges[source, target]["linkcolor"] = edge_color_for_path_group(
                    attrs.get("path_group", "Backbone"),
                    float(attrs.get("path_group_dominance", 0.0)),
                    faded=faded,
                    focused=focused,
                )
        scenes.append(scene)

    interlayer_edge_records = interlayer_records["edges"]
    clade_flow_edge_records = [
        record for record in interlayer_edge_records if record.get("edge_type") == "clade_flow_segment"
    ]
    taxon_flow_edge_records = [
        record for record in interlayer_edge_records if record.get("edge_type") == "taxon_flow_segment"
    ]
    interlayer_taxon_records = interlayer_records["shared_taxa"]
    interlayer_by_group = Counter(record["group"] for record in interlayer_taxon_records)
    direct_clade_flow_records = [
        record for record in interlayer_records["clade_flow_records"]
        if record.get("path_kind", "clade_flow") == "clade_flow"
    ]
    ancestor_context_flow_records = [
        record for record in interlayer_records["clade_flow_records"]
        if record.get("path_kind") == "ancestor_context_flow"
    ]
    clade_flows_by_pair = Counter(record["gene_pair"] for record in direct_clade_flow_records)
    ancestor_context_flows_by_pair = Counter(record["gene_pair"] for record in ancestor_context_flow_records)
    selected_taxon_flows_by_pair = Counter(record["gene_pair"] for record in interlayer_records["selected_taxon_flow_records"])
    adjacent_pair_taxa = {
        f"{source_gene}->{target_gene}": len(interlayer_taxon_records)
        for source_gene, target_gene in zip(pst.GENE_ORDER[:-1], pst.GENE_ORDER[1:])
    }
    audit = {
        "project_name": PROJECT_NAME,
        "scene_names": SCENE_NAMES,
        "gene_order": pst.GENE_ORDER,
        "nodes": base.number_of_nodes(),
        "edges": base.number_of_edges(),
        "within_gene_tree_edges": sum(graph.number_of_edges() for graph in gene_graphs.values()),
        "interlayer_species_connections": len(interlayer_records["selected_taxon_flow_records"]),
        "interlayer_routed_edges": len(interlayer_edge_records),
        "interlayer_clade_flow_edges": len(clade_flow_edge_records),
        "interlayer_selected_taxon_flow_edges": len(taxon_flow_edge_records),
        "interlayer_flow_nodes": len(interlayer_records["flow_nodes"]),
        "interlayer_direct_mrca_flow_segments": len(direct_clade_flow_records),
        "interlayer_ancestor_context_flow_segments": len(ancestor_context_flow_records),
        "interlayer_selected_taxon_flow_tracks": len(interlayer_records["selected_taxon_flow_records"]),
        "explicit_path_records": len(interlayer_records["path_records"]),
        "interlayer_shared_taxa": len(interlayer_taxon_records),
        "adjacent_gene_pair_species_connections": adjacent_pair_taxa,
        "visible_selected_taxon_flows_by_pair": dict(sorted(selected_taxon_flows_by_pair.items())),
        "clade_flow_corridors_by_pair": dict(sorted(clade_flows_by_pair.items())),
        "ancestor_context_flows_by_pair": dict(sorted(ancestor_context_flows_by_pair.items())),
        "interlayer_shared_taxa_by_group": dict(sorted(interlayer_by_group.items())),
        "interlayer_routing_policy": (
            f"clade-level corridors connect MRCA_{pst.GENE_ORDER[0]}(S) to "
            f"MRCA_{pst.GENE_ORDER[1]}(S) to MRCA_{pst.GENE_ORDER[2]}(S) "
            "for each shared taxon group S; additional ancestor_context_flow "
            "paths connect parent/grandparent ancestors around those MRCAs and "
            "are not labeled as direct most-recent common ancestors; faint "
            "selected same-taxon tracks remain separate"
        ),
        "per_gene": {},
        "taxon_policy": (
            "same data-backed ortholog species per gene; human subtree leaves are a "
            "neutral visual anchor. Genes without cached human variant trees reuse the "
            "existing ACTB subtree proxy and are marked in graph metadata."
        ),
        "datadivr_constraints": (
            "browser preview reads per-node radius sidecars; link thickness remains fixed "
            "in the current DataDiVR renderer"
        ),
        "node_radius_policy": (
            "leaves remain radius 1.0; internal tree nodes use log-scaled descendant-leaf "
            "support in [0.58, 1.50]; interlayer helper nodes stay deliberately small"
        ),
        "path_color_policy": (
            "links inherit the dominant descendant leaf group when that subtree reaches "
            f"{PATH_GROUP_DOMINANCE_THRESHOLD:.0%} dominance; ambiguous backbone links remain neutral"
        ),
    }
    for gene in pst.GENE_ORDER:
        graph = gene_graphs[gene]
        path_group_counts = Counter(
            assignment["path_group"] for assignment in path_groups_by_gene[gene].values()
        )
        audit["per_gene"][gene] = {
            "nodes": graph.number_of_nodes(),
            "edges": graph.number_of_edges(),
            "leaves": sum(1 for _, attrs in graph.nodes(data=True) if attrs.get("is_leaf")),
            "mammal_ortholog_leaves": pst.mammal_ortholog_count(graph),
            "human_subtree_leaves": pst.human_variant_leaf_count(graph),
            "tree_edge_crossings_xy": pst.count_tree_edge_crossings(graph, layouts[gene]),
            "is_tree": bool(nx.is_tree(pst.layout_graph(graph))),
            "path_group_link_counts": dict(sorted(path_group_counts.items())),
        }
    common_species = pst.common_inferred_homology_species()
    available_common_species = pst.available_common_inferred_homology_species()
    audit["common_inferred_species"] = len(common_species)
    audit["available_common_inferred_species"] = len(available_common_species)
    audit["max_displayed_ortholog_species"] = pst.MAX_DISPLAYED_ORTHOLOG_SPECIES
    total_leaves = {
        gene: sum(1 for _, attrs in gene_graphs[gene].nodes(data=True) if attrs.get("is_leaf"))
        for gene in pst.GENE_ORDER
    }
    audit["target_250_taxa_status"] = (
        "met_for_displayed_leaves_and_real_ortholog_species; "
        f"{len(common_species)} displayed ortholog species selected from "
        f"{len(available_common_species)} cached {'/'.join(pst.GENE_ORDER)} shared inferred species. "
        f"Displayed leaves per gene: {total_leaves}."
    )
    return scenes, audit


def encode_coordinate(value: float) -> tuple[int, int]:
    encoded = int(clamp(value, 0.0, 0.9999) * 65280)
    return encoded // 255, encoded % 255


def write_bmp_rgb(path: str, width: int, height: int, pixels: Iterable[tuple[int, int, int]]) -> None:
    rows = []
    pixel_list = list(pixels)
    row_stride = width * 3
    padding = (4 - (row_stride % 4)) % 4
    for row_index in range(height - 1, -1, -1):
        row = bytearray()
        for red, green, blue in pixel_list[row_index * width : (row_index + 1) * width]:
            row.extend([int(blue), int(green), int(red)])
        row.extend([0] * padding)
        rows.append(bytes(row))
    image_data = b"".join(rows)
    file_size = 14 + 40 + len(image_data)
    header = b"BM" + struct.pack("<IHHI", file_size, 0, 0, 54)
    dib = struct.pack("<IIIHHIIIIII", 40, width, height, 1, 24, 0, len(image_data), 2835, 2835, 0, 0)
    with open(path, "wb") as fh:
        fh.write(header)
        fh.write(dib)
        fh.write(image_data)


def png_chunk(kind: bytes, data: bytes) -> bytes:
    return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF)


def write_png_rgba(path: str, width: int, height: int, pixels: Iterable[tuple[int, int, int, int]]) -> None:
    pixel_list = list(pixels)
    rows = []
    for row_index in range(height):
        row = bytearray([0])
        for red, green, blue, alpha in pixel_list[row_index * width : (row_index + 1) * width]:
            row.extend([int(red), int(green), int(blue), int(alpha)])
        rows.append(bytes(row))
    raw = b"".join(rows)
    ihdr = struct.pack(">IIBBBBB", width, height, 8, 6, 0, 0, 0)
    with open(path, "wb") as fh:
        fh.write(b"\x89PNG\r\n\x1a\n")
        fh.write(png_chunk(b"IHDR", ihdr))
        fh.write(png_chunk(b"IDAT", zlib.compress(raw, 9)))
        fh.write(png_chunk(b"IEND", b""))


def layout_texture_pixels(positions: list[list[float]]) -> tuple[list[tuple[int, int, int]], list[tuple[int, int, int]], int]:
    height = 128 * (int(len(positions) / 16384) + 1)
    size = 128 * height
    high_pixels = [(0, 0, 0)] * size
    low_pixels = [(0, 0, 0)] * size
    for index, point in enumerate(positions):
        xh, xl = encode_coordinate(point[0])
        yh, yl = encode_coordinate(point[1])
        zh, zl = encode_coordinate(point[2])
        high_pixels[index] = (xh, yh, zh)
        low_pixels[index] = (xl, yl, zl)
    return high_pixels, low_pixels, height


def link_texture_pixels(edges: list[tuple[int, int]]) -> tuple[list[tuple[int, int, int]], int]:
    height = 64 * (int(len(edges) / 32768) + 1)
    pixels = [(0, 0, 0)] * (1024 * height)
    for index, (source, target) in enumerate(edges):
        pixels[index * 2] = (source % 128, int(source / 128) % 128, int(source / 16384))
        pixels[index * 2 + 1] = (target % 128, int(target / 128) % 128, int(target / 16384))
    return pixels, height


def rgba_texture_pixels(colors: list[tuple[int, int, int, int]], width: int, height: int) -> list[tuple[int, int, int, int]]:
    pixels = [(0, 0, 0, 0)] * (width * height)
    for index, color in enumerate(colors):
        pixels[index] = tuple(int(clamp(channel, 0, 255)) for channel in color)
    return pixels


def project_subdirs(project_dir: str) -> None:
    os.makedirs(project_dir, exist_ok=True)
    for dirname in ("layouts", "layoutsl", "layoutsRGB", "links", "linksRGB", "nodesizes", "legends", "analysis"):
        os.makedirs(os.path.join(project_dir, dirname), exist_ok=True)


def path_payload_from_scenes(scenes: list[nx.Graph], node_order: list[str]) -> dict:
    node_to_id = {node: index for index, node in enumerate(node_order)}
    records = []
    paths = []
    for index, record in enumerate(scenes[0].graph.get("path_records", [])):
        node_ids = [node_to_id[node] for node in record["path_nodes"] if node in node_to_id]
        if len(node_ids) < 2:
            continue
        path_record = {
            "id": index,
            "name": record.get("name", f"path {index}"),
            "kind": record.get("kind", ""),
            "gene_pair": record.get("gene_pair", ""),
            "group": record.get("group", ""),
            "flow_kind": record.get("flow_kind", ""),
            "context_level": record.get("context_level", 0),
            "label": record.get("label", ""),
            "score": record.get("score", None),
            "must_show": bool(record.get("must_show", False)),
            "nodes": node_ids,
            "source": node_ids[0],
            "target": node_ids[-1],
        }
        records.append(path_record)
        paths.append(node_ids)
    return {
        "projectname": PROJECT_NAME,
        "coordinate_system": "DataDiVR node ids; use layout textures or datadivr_coordinate_mappings.json for scene-specific coordinates",
        "animation_settings": PATH_ANIMATION_SETTINGS,
        "paths": paths,
        "path_records": records,
    }


def path_ids_payload(path_payload: dict) -> dict:
    return {"paths": path_payload["paths"]}


def path_connections_payload(path_payload: dict) -> dict:
    connections = []
    for record in path_payload["path_records"]:
        nodes = record["nodes"]
        connections.append(
            {
                "id": record["id"],
                "kind": record.get("kind", ""),
                "gene_pair": record.get("gene_pair", ""),
                "group": record.get("group", ""),
                "flow_kind": record.get("flow_kind", ""),
                "context_level": record.get("context_level", 0),
                "segments": [[source, target] for source, target in zip(nodes[:-1], nodes[1:])],
                "source": nodes[0],
                "target": nodes[-1],
            }
        )
    return {
        "projectname": PROJECT_NAME,
        "node_id_contract": "All path connections use numeric DataDiVR node ids.",
        "path_connections": connections,
        "path_records": path_payload["path_records"],
        "animation_settings": path_payload["animation_settings"],
    }


def coordinate_payload_from_scenes(
    scenes: list[nx.Graph],
    node_order: list[str],
    edge_order: list[tuple[str, str]],
) -> dict:
    node_to_id = {node: index for index, node in enumerate(node_order)}
    edge_to_id = {edge: index for index, edge in enumerate(edge_order)}
    base_nodes = []
    for node in node_order:
        attrs = scenes[0].nodes[node]
        annotation = attrs.get("annotation", {})
        base_nodes.append(
            {
                "id": node_to_id[node],
                "node_key": node,
                "name": attrs.get("name", node),
                "annotation": annotation,
                "gene": annotation.get("gene"),
                "source_node_id": annotation.get("node_id"),
                "node_radius": attrs.get("node_radius", 1.0),
            }
        )

    scene_payload = {}
    for scene in scenes:
        layout_name = scene.graph["layoutname"]
        scene_payload[layout_name] = {
            "nodes": [
                {
                    "id": node_to_id[node],
                    "node_key": node,
                    "position_unit_cube": scene.nodes[node]["pos"],
                    "position_preview": [
                        scene.nodes[node]["pos"][1] * -20.0,
                        scene.nodes[node]["pos"][2] * 20.0,
                        scene.nodes[node]["pos"][0] * 20.0,
                    ],
                    "rgba": scene.nodes[node]["nodecolor"],
                    "node_radius": scene.nodes[node].get("node_radius", 1.0),
                    "cluster": scene.nodes[node].get("cluster"),
                }
                for node in node_order
            ],
            "links": [
                {
                    "id": edge_to_id[(source, target)],
                    "source": node_to_id[source],
                    "target": node_to_id[target],
                    "source_key": source,
                    "target_key": target,
                    "rgba": scene.edges[source, target]["linkcolor"],
                }
                for source, target in edge_order
            ],
        }

    return {
        "projectname": PROJECT_NAME,
        "node_id_contract": "Numeric ids match nodes.json, pfile paths, paths.json, and all DataDiVR textures.",
        "position_unit_cube": "DataDiVR texture coordinates before preview transform.",
        "position_preview": "Preview-space coordinates after webGL_preview.js transform with scale=20.",
        "nodes": base_nodes,
        "scenes": scene_payload,
    }


def merged_json_from_scenes(
    scenes: list[nx.Graph],
    node_order: list[str],
    edge_order: list[tuple[str, str]],
    path_payload: dict,
    coordinate_payload: dict,
) -> dict:
    node_to_id = {node: index for index, node in enumerate(node_order)}
    return {
        "directed": False,
        "multigraph": False,
        "projectname": PROJECT_NAME,
        "info": INFO,
        "graphlayouts": [scene.graph["layoutname"] for scene in scenes],
        "annotationTypes": True,
        "paths": path_payload["paths"],
        "path_records": path_payload["path_records"],
        "path_animation_settings": path_payload["animation_settings"],
        "coordinate_mappings": coordinate_payload,
        "nodes": [
            {
                "id": node_to_id[node],
                "name": scenes[0].nodes[node].get("name", node),
                "annotation": scenes[0].nodes[node].get("annotation", {}),
                "node_radius": scenes[0].nodes[node].get("node_radius", 1.0),
            }
            for node in node_order
        ],
        "links": [
            {"id": index, "source": node_to_id[source], "target": node_to_id[target]}
            for index, (source, target) in enumerate(edge_order)
        ],
        "layouts": [
            {
                "layoutname": scene.graph["layoutname"],
                "nodes": [
                    {
                        "id": node_to_id[node],
                        "pos": scene.nodes[node]["pos"],
                        "nodecolor": scene.nodes[node]["nodecolor"],
                        "node_radius": scene.nodes[node].get("node_radius", 1.0),
                        "cluster": scene.nodes[node].get("cluster"),
                    }
                    for node in node_order
                ],
                "links": [
                    {
                        "source": node_to_id[source],
                        "target": node_to_id[target],
                        "linkcolor": scene.edges[source, target]["linkcolor"],
                    }
                    for source, target in edge_order
                ],
            }
            for scene in scenes
        ],
    }


def write_datadivr_project(scenes: list[nx.Graph]) -> dict:
    if os.path.exists(PROJECT_DIR):
        shutil.rmtree(PROJECT_DIR)
    project_subdirs(PROJECT_DIR)

    node_order = list(scenes[0].nodes())
    edge_order = list(scenes[0].edges())
    node_to_id = {node: index for index, node in enumerate(node_order)}
    edge_id_pairs = [(node_to_id[source], node_to_id[target]) for source, target in edge_order]
    path_payload = path_payload_from_scenes(scenes, node_order)
    paths_only_payload = path_ids_payload(path_payload)
    path_connections = path_connections_payload(path_payload)
    coordinate_payload = coordinate_payload_from_scenes(scenes, node_order, edge_order)

    nodes_json = {
        "nodes": [
            {
                "id": node_to_id[node],
                "n": scenes[0].nodes[node].get("name", node),
                "attrlist": scenes[0].nodes[node].get("annotation", {}),
            }
            for node in node_order
        ]
    }
    links_json = {
        "links": [
            {"id": index, "s": str(source), "e": str(target)}
            for index, (source, target) in enumerate(edge_id_pairs)
        ]
    }

    with open(os.path.join(PROJECT_DIR, "nodes.json"), "w") as fh:
        json.dump(nodes_json, fh, indent=2)
    with open(os.path.join(PROJECT_DIR, "links.json"), "w") as fh:
        json.dump(links_json, fh, indent=2)
    with open(os.path.join(PROJECT_DIR, "paths.json"), "w") as fh:
        json.dump(paths_only_payload, fh, indent=2)
    with open(os.path.join(PROJECT_DIR, "path_connections.json"), "w") as fh:
        json.dump(path_connections, fh, indent=2)
    with open(os.path.join(PROJECT_DIR, "coordinate_mappings.json"), "w") as fh:
        json.dump(coordinate_payload, fh, indent=2)

    for scene in scenes:
        name = scene.graph["layoutname"]
        positions = [scene.nodes[node]["pos"] for node in node_order]
        high_pixels, low_pixels, layout_height = layout_texture_pixels(positions)
        write_bmp_rgb(os.path.join(PROJECT_DIR, "layouts", f"{name}.bmp"), 128, layout_height, high_pixels)
        write_bmp_rgb(os.path.join(PROJECT_DIR, "layoutsl", f"{name}l.bmp"), 128, layout_height, low_pixels)

        node_colors = [scene.nodes[node]["nodecolor"] for node in node_order]
        node_color_pixels = rgba_texture_pixels(node_colors, 128, layout_height)
        write_png_rgba(os.path.join(PROJECT_DIR, "layoutsRGB", f"{name}.png"), 128, layout_height, node_color_pixels)

        node_sizes = [float(scene.nodes[node].get("node_radius", 1.0)) for node in node_order]
        with open(os.path.join(PROJECT_DIR, "nodesizes", f"{name}.json"), "w") as fh:
            json.dump(node_sizes, fh)

        link_color_height = 64 * (int(len(edge_order) / 32768) + 1)
        link_colors = [scene.edges[source, target]["linkcolor"] for source, target in edge_order]
        link_color_pixels = rgba_texture_pixels(link_colors, 512, link_color_height)
        write_png_rgba(os.path.join(PROJECT_DIR, "linksRGB", f"{name}.png"), 512, link_color_height, link_color_pixels)

    link_pixels, link_height = link_texture_pixels(edge_id_pairs)
    write_bmp_rgb(os.path.join(PROJECT_DIR, "links", f"{SCENE_NAMES[0]}.bmp"), 1024, link_height, link_pixels)

    pfile = {
        "name": PROJECT_NAME,
        "projectname": PROJECT_NAME,
        "layouts": SCENE_NAMES,
        "layoutsRGB": SCENE_NAMES,
        "links": [SCENE_NAMES[0]],
        "linksRGB": SCENE_NAMES,
        "nodeSizes": SCENE_NAMES,
        "paths": path_payload["paths"],
        "pathMetadataFile": "path_connections",
        "pathConnectionsFile": "path_connections",
        "pathAnimationSettings": PATH_ANIMATION_SETTINGS,
        "subtreeHighlightAnimation": SUBTREE_HIGHLIGHT_ANIMATION,
        "coordinateMappingsFile": "coordinate_mappings",
        "selections": named_group_selections(scenes, node_order),
        "linkSelections": named_link_group_selections(scenes, edge_order),
        "scenes": SCENE_NAMES,
        "labelcount": 0,
        "nodecount": len(node_order),
        "linkcount": len(edge_order),
        "labels": [len(node_order), 0],
        "annotationTypes": True,
        "info": INFO,
        "graphdesc": INFO,
        "legendfiles": [],
    }
    with open(os.path.join(PROJECT_DIR, "pfile.json"), "w") as fh:
        json.dump(pfile, fh, indent=2)
    with open(os.path.join(PROJECT_DIR, "pdata.json"), "w") as fh:
        json.dump(DEFAULT_PDATA, fh, indent=2)
    with open(os.path.join(PROJECT_DIR, "README.md"), "w") as fh:
        fh.write(
            "# Pangenome Housekeeping Stacked Trees\n\n"
            "This folder is a complete DataDiVR project. Copy the whole "
            "`Pangenome_Housekeeping_Stacked_Trees` directory into "
            "`DataDiVR_WebApp/static/projects/` and select it from the preview.\n\n"
            "Native DataDiVR files are in the project root plus the `layouts`, "
            "`layoutsl`, `layoutsRGB`, `links`, `linksRGB`, and `nodesizes` "
            "directories.\n\n"
            "Additional sidecars:\n\n"
            "- `paths.json` stores only explicit paths as numeric DataDiVR node IDs.\n"
            "- `path_connections.json` stores per-path segment pairs and metadata.\n"
            "- `coordinate_mappings.json` maps numeric node IDs to node keys, "
            "annotations, colors, and coordinates in each scene.\n"
            "- `pfile.json` includes `subtreeHighlightAnimation.presentation`, "
            "an 8-stage preview sequence: all layers, GAPDH, ENO1, RPLP0, "
            "the Ray-finned fish subtree in each layer, and finally the full "
            "Ray-finned fish subtree plus its inter-layer connections.\n"
            "- `analysis/` bundles the merged JSON export, NetworkX scene pickle, "
            "paths, coordinate mappings, and audit files for downstream agents.\n"
        )

    merged = merged_json_from_scenes(scenes, node_order, edge_order, path_payload, coordinate_payload)
    os.makedirs(os.path.dirname(OUTPUT_JSON), exist_ok=True)
    with open(OUTPUT_JSON, "w") as fh:
        json.dump(merged, fh, indent=2)
    with open(OUTPUT_PATHS_JSON, "w") as fh:
        json.dump(paths_only_payload, fh, indent=2)
    with open(OUTPUT_PATH_CONNECTIONS_JSON, "w") as fh:
        json.dump(path_connections, fh, indent=2)
    with open(OUTPUT_COORDINATES_JSON, "w") as fh:
        json.dump(coordinate_payload, fh, indent=2)
    with open(OUTPUT_SCENES_PICKLE, "wb") as fh:
        pickle.dump(scenes, fh)

    return {
        "node_count": len(node_order),
        "edge_count": len(edge_order),
        "scene_count": len(scenes),
        "project_dir": PROJECT_DIR,
        "merged_json": OUTPUT_JSON,
        "networkx_scenes": OUTPUT_SCENES_PICKLE,
        "paths_json": OUTPUT_PATHS_JSON,
        "path_connections_json": OUTPUT_PATH_CONNECTIONS_JSON,
        "coordinate_mappings": OUTPUT_COORDINATES_JSON,
        "path_count": len(path_payload["paths"]),
    }


def validate_project(summary: dict, scenes: list[nx.Graph]) -> dict:
    required = [
        os.path.join(PROJECT_DIR, "pfile.json"),
        os.path.join(PROJECT_DIR, "nodes.json"),
        os.path.join(PROJECT_DIR, "links.json"),
        os.path.join(PROJECT_DIR, "paths.json"),
        os.path.join(PROJECT_DIR, "path_connections.json"),
        os.path.join(PROJECT_DIR, "coordinate_mappings.json"),
    ]
    for scene_name in SCENE_NAMES:
        required.extend(
            [
                os.path.join(PROJECT_DIR, "layouts", f"{scene_name}.bmp"),
                os.path.join(PROJECT_DIR, "layoutsl", f"{scene_name}l.bmp"),
                os.path.join(PROJECT_DIR, "layoutsRGB", f"{scene_name}.png"),
                os.path.join(PROJECT_DIR, "linksRGB", f"{scene_name}.png"),
                os.path.join(PROJECT_DIR, "nodesizes", f"{scene_name}.json"),
            ]
        )
    required.append(os.path.join(PROJECT_DIR, "links", f"{SCENE_NAMES[0]}.bmp"))

    missing = [path for path in required if not os.path.exists(path)]
    position_values = []
    for scene in scenes:
        for _, attrs in scene.nodes(data=True):
            position_values.extend(attrs["pos"])
    finite_positions = all(math.isfinite(value) for value in position_values)
    in_unit_cube = all(0.0 <= value <= 1.0 for value in position_values)
    path_validation = validate_paths_file(
        os.path.join(PROJECT_DIR, "paths.json"),
        os.path.join(PROJECT_DIR, "path_connections.json"),
        os.path.join(PROJECT_DIR, "pfile.json"),
        scenes[0],
    )
    return {
        **summary,
        "missing_files": missing,
        "finite_positions": finite_positions,
        "positions_in_unit_cube": in_unit_cube,
        "min_position": min(position_values) if position_values else None,
        "max_position": max(position_values) if position_values else None,
        **path_validation,
    }


def validate_paths_file(paths_path: str, connections_path: str, pfile_path: str, scene: nx.Graph) -> dict:
    if not os.path.exists(paths_path) or not os.path.exists(connections_path) or not os.path.exists(pfile_path):
        return {
            "paths_valid": False,
            "path_validation_errors": ["paths.json, path_connections.json, or pfile.json missing"],
        }

    with open(paths_path) as fh:
        path_payload = json.load(fh)
    with open(connections_path) as fh:
        connection_payload = json.load(fh)
    with open(pfile_path) as fh:
        pfile = json.load(fh)

    node_order = list(scene.nodes())
    node_count = len(node_order)
    edge_pairs = {
        tuple(sorted((source, target)))
        for source, target in (
            (node_order.index(source), node_order.index(target))
            for source, target in scene.edges()
        )
    }
    errors = []
    paths = path_payload.get("paths", [])
    records = connection_payload.get("path_records", [])
    connections = connection_payload.get("path_connections", [])
    if pfile.get("paths") != paths:
        errors.append("pfile paths do not match paths.json paths")
    if pfile.get("pathMetadataFile") != "path_connections":
        errors.append("pfile pathMetadataFile does not point at path_connections")
    if len(paths) != len(records):
        errors.append("paths and path_records length mismatch")
    if len(paths) != len(connections):
        errors.append("paths and path_connections length mismatch")
    for index, path in enumerate(paths):
        if not isinstance(path, list) or len(path) < 2:
            errors.append(f"path {index} is not a list with at least two nodes")
            continue
        if not all(isinstance(node_id, int) for node_id in path):
            errors.append(f"path {index} contains non-integer node ids")
            continue
        missing_nodes = [node_id for node_id in path if node_id < 0 or node_id >= node_count]
        if missing_nodes:
            errors.append(f"path {index} references missing node ids: {missing_nodes[:5]}")
        for source, target in zip(path[:-1], path[1:]):
            if tuple(sorted((source, target))) not in edge_pairs:
                errors.append(f"path {index} segment {source}->{target} is not a project link")
        if index < len(records):
            record_nodes = records[index].get("nodes")
            if record_nodes != path:
                errors.append(f"path {index} record nodes do not match path")
            if records[index].get("source") != path[0] or records[index].get("target") != path[-1]:
                errors.append(f"path {index} source/target metadata does not match path")
        if index < len(connections):
            expected_segments = [[source, target] for source, target in zip(path[:-1], path[1:])]
            if connections[index].get("segments") != expected_segments:
                errors.append(f"path {index} connection segments do not match path")

    kind_counts = Counter(record.get("kind", "") for record in records)
    gene_pair_counts = Counter(record.get("gene_pair", "") for record in records)
    return {
        "paths_valid": not errors,
        "path_validation_errors": errors[:20],
        "path_validation_error_count": len(errors),
        "validated_path_count": len(paths),
        "validated_path_kind_counts": dict(kind_counts),
        "validated_path_gene_pair_counts": dict(gene_pair_counts),
    }


def write_project_analysis_bundle(audit: dict) -> None:
    analysis_dir = os.path.join(PROJECT_DIR, "analysis")
    os.makedirs(analysis_dir, exist_ok=True)
    manifest = {
        "projectname": PROJECT_NAME,
        "purpose": "Self-contained DataDiVR handoff files for downstream agents and notebooks.",
        "files": {
            f"{PROJECT_NAME}.json": "Merged human-readable NetworkX/DataDiVR export.",
            f"{PROJECT_NAME}_networkx_scenes.pkl": "Pickled NetworkX scene graphs with node, edge, layout, and flow metadata.",
            f"{PROJECT_NAME}_paths.json": "Explicit numeric DataDiVR node-id paths only.",
            f"{PROJECT_NAME}_path_connections.json": "Per-path node-id segment pairs plus path metadata.",
            "datadivr_coordinate_mappings.json": "Per-scene node ids, source node keys, coordinates, colors, and edge ids.",
            f"{PROJECT_NAME}_datadivr_audit.json": "Generation and validation audit.",
        },
    }
    with open(os.path.join(analysis_dir, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    with open(os.path.join(analysis_dir, "README.md"), "w") as fh:
        fh.write(
            "# DataDiVR analysis bundle\n\n"
            "This directory keeps the non-native analysis sidecars inside the "
            "DataDiVR project. DataDiVR renders the texture and JSON files in "
            "the project root; scripts and agents should use these sidecars for "
            "path semantics, coordinate lookup, and NetworkX-level inspection.\n\n"
            "- `../paths.json` and `Pangenome_Housekeeping_Stacked_Trees_paths.json` "
            "contain only numeric DataDiVR node-id paths.\n"
            "- `../path_connections.json` and `Pangenome_Housekeeping_Stacked_Trees_path_connections.json` "
            "contain segment pairs and path metadata.\n"
            "- `datadivr_coordinate_mappings.json` maps those ids back to node keys "
            "and per-scene coordinates.\n"
            "- `Pangenome_Housekeeping_Stacked_Trees_networkx_scenes.pkl` contains "
            "the full NetworkX scene objects.\n"
        )
    for path in (
        OUTPUT_JSON,
        OUTPUT_SCENES_PICKLE,
        OUTPUT_PATHS_JSON,
        OUTPUT_PATH_CONNECTIONS_JSON,
        OUTPUT_COORDINATES_JSON,
        OUTPUT_AUDIT,
    ):
        if os.path.exists(path):
            shutil.copy2(path, os.path.join(analysis_dir, os.path.basename(path)))
    with open(os.path.join(analysis_dir, "audit_summary.json"), "w") as fh:
        json.dump(audit, fh, indent=2)


def mirror_portable_datadivr_project() -> None:
    os.makedirs(os.path.dirname(PORTABLE_PROJECT_DIR), exist_ok=True)
    if os.path.exists(PORTABLE_PROJECT_DIR):
        shutil.rmtree(PORTABLE_PROJECT_DIR)
    shutil.copytree(PROJECT_DIR, PORTABLE_PROJECT_DIR)


def main() -> None:
    scenes, audit = build_combined_graphs()
    summary = write_datadivr_project(scenes)
    audit.update(validate_project(summary, scenes))
    audit["portable_project_dir"] = PORTABLE_PROJECT_DIR
    with open(OUTPUT_AUDIT, "w") as fh:
        json.dump(audit, fh, indent=2)
    write_project_analysis_bundle(audit)
    mirror_portable_datadivr_project()
    with open(OUTPUT_AUDIT, "w") as fh:
        json.dump(audit, fh, indent=2)
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
