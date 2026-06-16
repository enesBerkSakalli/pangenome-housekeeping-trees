"""
Render stacked 3D unrooted human-subtree gene trees.

This is deliberately tree-first:
- the topology is treated as an undirected unrooted tree
- coordinates are assigned with an equal-angle phylogenetic tree layout
- affine box fitting preserves the non-crossing tree embedding
- each gene occupies a stacked 3D cuboid; local z encodes branch-depth
"""

from __future__ import annotations

import json
import math
import os
import pickle
from collections import defaultdict

import networkx as nx
import numpy as np
import plotly.graph_objects as go
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform

ROOT = os.path.dirname(os.path.abspath(__file__))
GRAPH_BUNDLE_PATH = os.path.join(
    ROOT,
    "outputs_3d",
    "networks_population_expanded",
    "all_population_expanded_trees.pkl",
)

OUTPUT_MAIN = os.path.join(ROOT, "outputs_reference", "plotly_reference_population_stack.html")
OUTPUT_UNROOTED = os.path.join(ROOT, "outputs_3d", "plotly_unrooted_3d_variant_trees.html")
OUTPUT_STACKED = os.path.join(ROOT, "outputs_3d", "plotly_layered_stacked_trees.html")
OUTPUT_HYPERBOLIC = os.path.join(ROOT, "outputs_3d", "plotly_hyperbolic_stacked_trees.html")
LAYOUT_AUDIT_PATH = os.path.join(ROOT, "outputs_3d", "stacked_3d_unrooted_layout_audit.json")
HYPERBOLIC_AUDIT_PATH = os.path.join(ROOT, "outputs_3d", "hyperbolic_stacked_tree_layout_audit.json")
MAMMAL_CONTEXT_PATH = os.path.join(ROOT, "outputs_3d", "mammal_ortholog_context.json")
HOMOLOGY_ALIGNMENT_PATHS = {
    "GAPDH": os.path.join(ROOT, "outputs_3d", "gapdh_all_homologies.json"),
    "ACTB": os.path.join(ROOT, "outputs_3d", "actb_all_homologies.json"),
    "ENO1": os.path.join(ROOT, "outputs_3d", "eno1_all_homologies.json"),
    "RPLP0": os.path.join(ROOT, "outputs_3d", "rplp0_all_homologies.json"),
}
NCBI_HOMOLOGY_ALIGNMENT_PATHS = {
    "GAPDH": os.path.join(ROOT, "outputs_3d", "ncbi_gapdh_all_homologies.json"),
    "ACTB": os.path.join(ROOT, "outputs_3d", "ncbi_actb_all_homologies.json"),
    "ENO1": os.path.join(ROOT, "outputs_3d", "ncbi_eno1_all_homologies.json"),
    "RPLP0": os.path.join(ROOT, "outputs_3d", "ncbi_rplp0_all_homologies.json"),
}

GENE_ORDER = ["GAPDH", "ENO1", "RPLP0"]
LAYER_Z_SEPARATION = 980.0
GENE_Z = {"GAPDH": -LAYER_Z_SEPARATION, "ENO1": 0.0, "RPLP0": LAYER_Z_SEPARATION}
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

MAMMAL_CLADE_ORDER = [
    "Great apes",
    "Lesser apes",
    "Old World monkeys",
    "New World monkeys",
    "Strepsirrhines",
    "Tarsiiformes",
    "Scandentia",
    "Glires",
    "Carnivores",
    "Cetartiodactyla",
    "Perissodactyla",
    "Bats",
    "Eulipotyphla",
    "Afrotheria",
    "Xenarthra",
    "Monotremata",
    "Marsupialia",
    "Birds",
    "Reptiles",
    "Amphibians",
    "Ray-finned fish",
    "Cartilaginous fish",
    "Lobe-finned fish",
    "Jawless fish",
    "Fungi outgroup",
]
MAMMAL_CLADE_COLORS = {
    "Great apes": "#4ADE80",
    "Lesser apes": "#22C55E",
    "Old World monkeys": "#16A34A",
    "New World monkeys": "#06B6D4",
    "Strepsirrhines": "#2DD4BF",
    "Tarsiiformes": "#14B8A6",
    "Scandentia": "#84CC16",
    "Glires": "#38BDF8",
    "Carnivores": "#F97316",
    "Cetartiodactyla": "#FACC15",
    "Perissodactyla": "#FB923C",
    "Bats": "#A78BFA",
    "Eulipotyphla": "#60A5FA",
    "Afrotheria": "#F472B6",
    "Xenarthra": "#EF4444",
    "Monotremata": "#F9A8D4",
    "Marsupialia": "#C084FC",
    "Birds": "#B8F2E6",
    "Reptiles": "#A3E635",
    "Amphibians": "#10B981",
    "Ray-finned fish": "#22D3EE",
    "Cartilaginous fish": "#38BDF8",
    "Lobe-finned fish": "#0EA5E9",
    "Jawless fish": "#2563EB",
    "Fungi outgroup": "#D6D3D1",
    "Species ortholog context": "#E5E7EB",
}
SPECIAL_GROUP_COLORS = {
    "Human lineage": "#FFF2CC",
    "Human subtree": "#F4F8FF",
    "Gene-history flow": "#CFE6FF",
}
HUMAN_SUBTREE_COLOR = "#FFF2CC"
MAMMAL_LAYOUT_WEIGHT = 18

EXTRA_BIRD_SPECIES = {
    "meleagris_gallopavo": "Turkey",
    "parus_major": "Great tit",
    "serinus_canaria": "Canary",
}
EXTRA_REPTILE_SPECIES = {
    "notechis_scutatus": "Tiger snake",
}
EXTRA_FUNGAL_SPECIES = {
    "saccharomyces_cerevisiae": "Baker's yeast",
}
FOCAL_SPECIES_LABELS = {
    "pan_troglodytes",
    "pan_paniscus",
    "gorilla_gorilla",
    "pongo_abelii",
    "pongo_pygmaeus",
    "nomascus_leucogenys",
    "macaca_mulatta",
    "callithrix_jacchus",
    "mus_spicilegus",
    "rattus_norvegicus",
    "canis_lupus_familiaris",
    "panthera_leo",
    "bos_taurus",
    "tursiops_truncatus",
    "equus_asinus",
    "ornithorhynchus_anatinus",
    "notamacropus_eugenii",
    "anser_brachyrhynchus",
    "aquila_chrysaetos_chrysaetos",
    "crocodylus_porosus",
    "naja_naja",
    "xenopus_tropicalis",
    "danio_rerio",
    "salmo_salar",
    "latimeria_chalumnae",
    "petromyzon_marinus",
    "saccharomyces_cerevisiae",
}

AF_KEYS = [f"af_{pop.lower().replace(' ', '_').replace('-', '_')}" for pop in POP_ORDER]
MIX_RATIO_THRESHOLD = 0.35
MIX_MIN_SECONDARY_AF = 0.001
HOMOLOGY_MIN_COVERAGE = 0.25
HUMAN_SEQUENCE_LEAF = "homo_sapiens"
MAX_DISPLAYED_ORTHOLOG_SPECIES = 700
FAST_LINKAGE_TREE_THRESHOLD = 360

LAYER_PLANE_RADIUS = 430.0
TREE_TARGET_RADIUS = 390.0
BOX_HALF_SIZE = 430.0
BOX_FILL_FRACTION = 0.92
LOCAL_Z_SPAN = 360.0
HYPERBOLIC_RADIUS = 430.0
HYPERBOLIC_FILL_FRACTION = 0.982
FLOW_CURVE_STEPS = 42
FLOW_OUTWARD_OFFSET = 68.0
HUMAN_MARKERS_PER_POPULATION = 42
HUMAN_SUBTREE_LEAVES_PER_POPULATION = 19
HISTORY_FLOW_TOP_PER_PAIR = 36
HISTORY_FLOW_SCORE_IDENTITY_WEIGHT = 7.0
NODE_MARKER_SIZE = 6.2
INTERNAL_NODE_MARKER_MIN = 3.8
INTERNAL_NODE_MARKER_MAX = 8.9
HYPERBOLIC_LEAF_Z_STAGGER = 34.0


def rgba_from_hex(hex_color: str, alpha: float) -> str:
    value = hex_color.lstrip("#")
    red = int(value[0:2], 16)
    green = int(value[2:4], 16)
    blue = int(value[4:6], 16)
    return f"rgba({red}, {green}, {blue}, {alpha})"


def load_graphs() -> dict[str, nx.DiGraph]:
    if not os.path.exists(GRAPH_BUNDLE_PATH):
        import population_expanded_trees

        population_expanded_trees.main()
    with open(GRAPH_BUNDLE_PATH, "rb") as fh:
        return pickle.load(fh)


def load_mammal_context() -> dict:
    if not os.path.exists(MAMMAL_CONTEXT_PATH):
        import mammal_ortholog_context

        mammal_ortholog_context.main()
    with open(MAMMAL_CONTEXT_PATH) as fh:
        return json.load(fh)


def common_available_mammal_species(context: dict) -> set[str]:
    species_sets = []
    for gene in GENE_ORDER:
        species_sets.append(
            {
                str(record.get("species"))
                for record in context.get("genes", {}).get(gene, [])
                if record.get("available") and record.get("species")
            }
        )
    if not species_sets:
        return set()
    return set.intersection(*species_sets)


def load_homology_payload(gene: str) -> list[dict]:
    path = NCBI_HOMOLOGY_ALIGNMENT_PATHS.get(gene)
    if not path or not os.path.exists(path):
        path = HOMOLOGY_ALIGNMENT_PATHS.get(gene)
    if not path or not os.path.exists(path):
        return []
    with open(path) as fh:
        payload = json.load(fh)
    data = payload.get("data") or []
    if not data:
        return []
    return list(data[0].get("homologies") or [])


def ungapped_sequence(sequence: str) -> str:
    return "".join(char for char in sequence if char != "-")


def homology_coverage(homology: dict) -> float:
    source = (homology.get("source") or {}).get("align_seq") or ""
    target = (homology.get("target") or {}).get("align_seq") or ""
    source_positions = sum(1 for char in source if char != "-")
    if source_positions <= 0:
        return 0.0
    shared_positions = sum(
        1
        for source_char, target_char in zip(source, target)
        if source_char != "-" and target_char != "-"
    )
    return shared_positions / source_positions


def best_homologies_for_gene(gene: str) -> dict[str, dict]:
    best: dict[str, tuple[tuple[float, ...], dict]] = {}
    for homology in load_homology_payload(gene):
        target = homology.get("target") or {}
        species = str(target.get("species") or "")
        if not species or not target.get("align_seq"):
            continue
        coverage = homology_coverage(homology)
        rank = (
            0.0 if homology.get("type") == "ortholog_one2one" else 1.0,
            -coverage,
            -float(target.get("perc_id") or 0.0),
            -float(target.get("perc_pos") or 0.0),
            -float(len(target.get("align_seq") or "")),
        )
        if species not in best or rank < best[species][0]:
            best[species] = (rank, homology)
    return {species: homology for species, (_rank, homology) in best.items()}


def available_common_inferred_homology_species() -> set[str]:
    species_sets = []
    for gene in GENE_ORDER:
        best = best_homologies_for_gene(gene)
        species_sets.append(
            {
                species
                for species, homology in best.items()
                if homology_coverage(homology) >= HOMOLOGY_MIN_COVERAGE
            }
        )
    if not species_sets:
        return set()
    return set.intersection(*species_sets)


def common_inferred_homology_species() -> set[str]:
    best_by_gene = {gene: best_homologies_for_gene(gene) for gene in GENE_ORDER}
    common = {
        species
        for species in set.intersection(
            *[
                {
                    species
                    for species, homology in best.items()
                    if homology_coverage(homology) >= HOMOLOGY_MIN_COVERAGE
                }
                for best in best_by_gene.values()
            ]
        )
    } if best_by_gene else set()
    if len(common) <= MAX_DISPLAYED_ORTHOLOG_SPECIES:
        return common

    def clade_for_species(species: str) -> str:
        target = (best_by_gene[GENE_ORDER[0]].get(species) or {}).get("target") or {}
        return str(target.get("clade") or "Species ortholog context")

    def evidence_score(species: str) -> tuple[float, float, str]:
        identities = []
        coverages = []
        for gene in GENE_ORDER:
            homology = best_by_gene[gene][species]
            target = homology.get("target") or {}
            identities.append(float(target.get("perc_id") or 0.0))
            coverages.append(homology_coverage(homology))
        return (min(identities), min(coverages), species)

    buckets: dict[str, list[str]] = defaultdict(list)
    for species in common:
        buckets[clade_for_species(species)].append(species)
    for bucket in buckets.values():
        bucket.sort(key=evidence_score, reverse=True)

    selected = []
    selected_set = set()
    focal_common = sorted(species for species in FOCAL_SPECIES_LABELS if species in common)
    for species in focal_common:
        selected.append(species)
        selected_set.add(species)
    clade_order = [*MAMMAL_CLADE_ORDER, "Species ortholog context"]
    while len(selected) < MAX_DISPLAYED_ORTHOLOG_SPECIES:
        added = False
        for clade in clade_order:
            bucket = buckets.get(clade, [])
            while bucket and bucket[0] in selected_set:
                bucket.pop(0)
            if not bucket:
                continue
            species = bucket.pop(0)
            selected.append(species)
            selected_set.add(species)
            added = True
            if len(selected) >= MAX_DISPLAYED_ORTHOLOG_SPECIES:
                break
        if not added:
            break
    return set(selected)


def context_species_metadata(context: dict) -> dict[str, dict]:
    metadata = {}
    for records in context.get("genes", {}).values():
        for record in records:
            species = record.get("species")
            if species and species not in metadata:
                metadata[str(species)] = record
    return metadata


def species_label_from_id(species: str) -> str:
    return " ".join(part.capitalize() for part in species.split("_"))


def species_display_metadata(species: str, homology: dict, context_metadata: dict[str, dict]) -> dict:
    target = homology.get("target") or {}
    context_record = context_metadata.get(species, {})
    if target.get("label") or target.get("clade"):
        label = target.get("label") or context_record.get("label") or species_label_from_id(species)
        clade = target.get("clade") or context_record.get("clade") or "Species ortholog context"
    elif species in EXTRA_BIRD_SPECIES:
        label = EXTRA_BIRD_SPECIES[species]
        clade = "Birds"
    elif species in EXTRA_REPTILE_SPECIES:
        label = EXTRA_REPTILE_SPECIES[species]
        clade = "Reptiles"
    elif species in EXTRA_FUNGAL_SPECIES:
        label = EXTRA_FUNGAL_SPECIES[species]
        clade = "Fungi outgroup"
    else:
        label = context_record.get("label") or species_label_from_id(species)
        clade = context_record.get("clade") or "Ray-finned fish"
    return {
        "species": species,
        "label": label,
        "clade": clade,
        "available": True,
        "orthology_type": homology.get("type", ""),
        "taxonomy_level": homology.get("taxonomy_level", ""),
        "ensembl_gene_id": target.get("id", ""),
        "protein_id": target.get("protein_id", ""),
        "taxon_id": target.get("taxon_id", ""),
        "percent_identity": float(target.get("perc_id") or 0.0),
        "percent_positive": float(target.get("perc_pos") or 0.0),
        "alignment_coverage": homology_coverage(homology),
    }


def human_reference_sequence(homologies: list[dict]) -> str:
    source_sequences = [
        ungapped_sequence((homology.get("source") or {}).get("align_seq") or "")
        for homology in homologies
        if (homology.get("source") or {}).get("align_seq")
    ]
    if not source_sequences:
        return ""
    return max(source_sequences, key=lambda sequence: (len(sequence), sequence))


def target_sequence_on_human_axis(homology: dict, human_length: int) -> str:
    source = (homology.get("source") or {}).get("align_seq") or ""
    target = (homology.get("target") or {}).get("align_seq") or ""
    mapped = ["-"] * human_length
    source_position = -1
    for source_char, target_char in zip(source, target):
        if source_char == "-":
            continue
        source_position += 1
        if source_position >= human_length:
            continue
        mapped[source_position] = target_char if target_char and target_char != "-" else "-"
    return "".join(mapped)


def sequence_distance(left: str, right: str) -> float:
    comparable = 0
    differences = 0
    for left_char, right_char in zip(left, right):
        if left_char == "-" or right_char == "-":
            continue
        comparable += 1
        if left_char != right_char:
            differences += 1
    if comparable == 0:
        return 1.0
    return differences / comparable


def build_neighbor_joining_graph(
    gene: str,
    leaf_attrs: dict[str, dict],
    sequences: dict[str, str],
) -> tuple[nx.Graph, str]:
    """Infer a fully resolved unrooted distance tree with neighbor joining."""
    nodes = sorted(sequences)
    distances: dict[tuple[str, str], float] = {}

    def key(left: str, right: str) -> tuple[str, str]:
        return (left, right) if left < right else (right, left)

    def get_distance(left: str, right: str) -> float:
        if left == right:
            return 0.0
        return distances[key(left, right)]

    def set_distance(left: str, right: str, value: float) -> None:
        if left != right:
            distances[key(left, right)] = max(float(value), 1e-9)

    for index, left in enumerate(nodes):
        for right in nodes[index + 1 :]:
            set_distance(left, right, sequence_distance(sequences[left], sequences[right]))

    tree = nx.Graph()
    for node in nodes:
        tree.add_node(node, **leaf_attrs[node])

    active = list(nodes)
    internal_index = 0
    while len(active) > 3:
        n_active = len(active)
        totals = {
            node: sum(get_distance(node, other) for other in active if other != node)
            for node in active
        }
        best_pair = None
        best_score = math.inf
        for left_index, left in enumerate(active):
            for right in active[left_index + 1 :]:
                score = (n_active - 2) * get_distance(left, right) - totals[left] - totals[right]
                tie_break = (score, get_distance(left, right), left, right)
                if best_pair is None or tie_break < (best_score, get_distance(best_pair[0], best_pair[1]), best_pair[0], best_pair[1]):
                    best_score = score
                    best_pair = (left, right)
        left, right = best_pair
        left_right_distance = get_distance(left, right)
        delta = (totals[left] - totals[right]) / max(n_active - 2, 1)
        left_length = max(1e-9, 0.5 * (left_right_distance + delta))
        right_length = max(1e-9, 0.5 * (left_right_distance - delta))
        internal = f"{gene}_ortholog_nj_{internal_index}"
        internal_index += 1
        tree.add_node(
            internal,
            label=f"{gene} inferred split {internal_index}",
            is_leaf=False,
            is_mammal_context=True,
            is_species_backbone=True,
            mammal_clade="Species ortholog context",
        )
        tree.add_edge(internal, left, weight=left_length)
        tree.add_edge(internal, right, weight=right_length)
        for other in active:
            if other in (left, right):
                continue
            set_distance(
                internal,
                other,
                0.5 * (get_distance(left, other) + get_distance(right, other) - left_right_distance),
            )
        active = [node for node in active if node not in (left, right)]
        active.append(internal)
        active.sort()

    root = f"{gene}_ortholog_nj_root"
    tree.add_node(
        root,
        label=f"{gene} unrooted inferred center",
        is_leaf=False,
        is_mammal_context=True,
        is_species_backbone=True,
        mammal_clade="Species ortholog context",
    )
    if len(active) == 2:
        left, right = active
        length = max(1e-9, get_distance(left, right) / 2.0)
        tree.add_edge(root, left, weight=length)
        tree.add_edge(root, right, weight=length)
    else:
        first, second, third = active
        first_length = max(
            1e-9,
            0.5 * (get_distance(first, second) + get_distance(first, third) - get_distance(second, third)),
        )
        second_length = max(
            1e-9,
            0.5 * (get_distance(first, second) + get_distance(second, third) - get_distance(first, third)),
        )
        third_length = max(
            1e-9,
            0.5 * (get_distance(first, third) + get_distance(second, third) - get_distance(first, second)),
        )
        tree.add_edge(root, first, weight=first_length)
        tree.add_edge(root, second, weight=second_length)
        tree.add_edge(root, third, weight=third_length)
    return tree, root


def build_linkage_distance_graph(
    gene: str,
    leaf_attrs: dict[str, dict],
    sequences: dict[str, str],
) -> tuple[nx.Graph, str]:
    """Infer a scalable fully resolved distance tree for large ortholog sets."""
    nodes = sorted(sequences)
    matrix = np.zeros((len(nodes), len(nodes)), dtype=float)
    for left_index, left in enumerate(nodes):
        for right_index in range(left_index + 1, len(nodes)):
            right = nodes[right_index]
            distance = sequence_distance(sequences[left], sequences[right])
            matrix[left_index, right_index] = distance
            matrix[right_index, left_index] = distance
    condensed = squareform(matrix, checks=False)
    linked = linkage(condensed, method="average", optimal_ordering=True)
    scipy_root = to_tree(linked, rd=False)

    tree = nx.Graph()
    for node in nodes:
        tree.add_node(node, **leaf_attrs[node])

    internal_index = {"value": 0}

    def add_cluster(cluster) -> str:
        if cluster.is_leaf():
            return nodes[cluster.id]
        internal_index["value"] += 1
        internal = f"{gene}_ortholog_linkage_{internal_index['value']}"
        tree.add_node(
            internal,
            label=f"{gene} inferred distance split {internal_index['value']}",
            is_leaf=False,
            is_mammal_context=True,
            is_species_backbone=True,
            mammal_clade="Species ortholog context",
        )
        left = add_cluster(cluster.left)
        right = add_cluster(cluster.right)
        left_length = max(1e-9, float(cluster.dist - getattr(cluster.left, "dist", 0.0)) / 2.0)
        right_length = max(1e-9, float(cluster.dist - getattr(cluster.right, "dist", 0.0)) / 2.0)
        tree.add_edge(internal, left, weight=left_length)
        tree.add_edge(internal, right, weight=right_length)
        return internal

    root = add_cluster(scipy_root)
    tree.nodes[root].update(label=f"{gene} inferred distance center")
    return tree, root


def orient_unrooted_tree(tree: nx.Graph, root: str) -> nx.DiGraph:
    oriented = nx.DiGraph()
    stack = [(root, None)]
    while stack:
        node, parent = stack.pop()
        oriented.add_node(node, **tree.nodes[node])
        if parent is not None:
            oriented.add_edge(parent, node, weight=float(tree.edges[parent, node].get("weight", 1e-9)))
        children = [neighbor for neighbor in tree.neighbors(node) if neighbor != parent]
        for child in sorted(children, reverse=True):
            stack.append((child, node))
    return oriented


def selected_human_variant_records(source_graph: nx.DiGraph, population: str, limit: int) -> list[tuple[str, dict]]:
    records = []
    for node, attrs in source_graph.nodes(data=True):
        if not attrs.get("is_leaf") or attrs.get("is_mammal_ortholog"):
            continue
        if attrs.get("assigned_population") != population:
            continue
        records.append(
            (
                str(node),
                dict(attrs),
            )
        )
    records.sort(
        key=lambda item: (
            -float(item[1].get("assigned_af", 0.0)),
            str(item[1].get("variant_id", item[1].get("label", item[0]))),
        )
    )
    return records[:limit]


def attach_human_population_summaries(graph: nx.DiGraph, source_graph: nx.DiGraph, gene: str, human_anchor: str) -> None:
    graph.nodes[human_anchor].update(
        label="Human variation subtree",
        is_leaf=False,
        is_population_root=False,
        is_human_variant_context=True,
        is_human_sequence_anchor=True,
        taxon_count=len(POP_ORDER),
        assigned_population="",
        population_label="",
        layout_weight=sum(population_summary_stats(source_graph, population)["count"] for population in POP_ORDER),
    )
    split_index = {"value": 0}
    variant_index = {"value": 0}

    def add_variant_leaf(parent: str, population: str, source_node: str, source_attrs: dict) -> str:
        variant_index["value"] += 1
        variant_node = f"{gene}_human_variant_{variant_index['value']:03d}"
        assigned_af = float(source_attrs.get("assigned_af", 0.0))
        variant_id = str(source_attrs.get("variant_id", source_attrs.get("label", source_node)))
        graph.add_node(
            variant_node,
            label=f"Human variant leaf {variant_index['value']:03d}",
            is_leaf=True,
            taxon_count=1,
            is_population_root=False,
            is_human_variant_leaf=True,
            assigned_population=population,
            population_label="",
            source_variant_node=source_node,
            variant_id=variant_id,
            assigned_af=assigned_af,
            layout_weight=6.0,
        )
        branch_length = 0.055 + (0.105 * (1.0 - min(assigned_af * 40.0, 1.0)))
        graph.add_edge(parent, variant_node, weight=branch_length)
        return variant_node

    def attach_variant_records(parent: str, population: str, records: list[tuple[str, dict]], depth: int = 0) -> None:
        if not records:
            return
        if len(records) == 1:
            source_node, source_attrs = records[0]
            add_variant_leaf(parent, population, source_node, source_attrs)
            return
        midpoint = len(records) // 2
        split_index["value"] += 1
        split_node = f"{gene}_human_variant_split_{split_index['value']}"
        graph.add_node(
            split_node,
            label=f"Human variant split {split_index['value']}",
            is_leaf=False,
            taxon_count=len(records),
            is_human_variant_context=True,
            mammal_clade="Primates",
            layout_weight=sum(float(attrs.get("layout_weight", 6.0)) for _, attrs in records),
        )
        graph.add_edge(parent, split_node, weight=0.070 + (0.010 * depth))
        attach_variant_records(split_node, population, records[:midpoint], depth + 1)
        attach_variant_records(split_node, population, records[midpoint:], depth + 1)

    def attach_binary(parent: str, populations: list[str], depth: int = 0) -> None:
        if not populations:
            return
        if len(populations) == 1:
            records = selected_human_variant_records(
                source_graph,
                populations[0],
                HUMAN_SUBTREE_LEAVES_PER_POPULATION,
            )
            attach_variant_records(parent, populations[0], records, depth)
            return
        midpoint = len(populations) // 2
        for side, group in (("left", populations[:midpoint]), ("right", populations[midpoint:])):
            if len(group) == 1:
                records = selected_human_variant_records(
                    source_graph,
                    group[0],
                    HUMAN_SUBTREE_LEAVES_PER_POPULATION,
                )
                attach_variant_records(parent, group[0], records, depth)
                continue
            split_index["value"] += 1
            split_node = f"{gene}_human_population_split_{split_index['value']}_{side}"
            graph.add_node(
                split_node,
                label=f"Human variation split {split_index['value']}",
                is_leaf=False,
                taxon_count=len(group),
                is_human_variant_context=True,
                mammal_clade="Primates",
                layout_weight=sum(population_summary_stats(source_graph, population)["count"] for population in group),
            )
            graph.add_edge(parent, split_node, weight=0.12 + (0.02 * depth))
            attach_binary(split_node, group, depth + 1)

    attach_binary(human_anchor, POP_ORDER)


def compose_inferred_taxon_aligned_tree(graph: nx.DiGraph, gene: str, context: dict) -> nx.DiGraph:
    """Build a same-taxa tree inferred from aligned Ensembl ortholog proteins."""
    best = best_homologies_for_gene(gene)
    common_species = sorted(common_inferred_homology_species())
    if not common_species:
        return compose_taxon_aligned_teaching_tree(graph, gene, context)

    context_metadata = context_species_metadata(context)
    homologies = [best[species] for species in common_species if species in best]
    human_sequence = human_reference_sequence(homologies)
    if not human_sequence:
        return compose_taxon_aligned_teaching_tree(graph, gene, context)

    sequences: dict[str, str] = {}
    leaf_attrs: dict[str, dict] = {}
    human_anchor = f"{gene}_human_sequence_anchor"
    sequences[human_anchor] = human_sequence
    leaf_attrs[human_anchor] = {
        "label": "Human sequence anchor",
        "is_leaf": True,
        "taxon_count": 1,
        "is_human_sequence_anchor": True,
        "is_mammal_context": True,
        "mammal_clade": "Great apes",
        "species": HUMAN_SEQUENCE_LEAF,
        "species_label": "Human",
        "layout_weight": MAMMAL_LAYOUT_WEIGHT,
    }

    records = []
    for species in common_species:
        homology = best.get(species)
        if not homology or homology_coverage(homology) < HOMOLOGY_MIN_COVERAGE:
            continue
        record = species_display_metadata(species, homology, context_metadata)
        species_node = f"{gene}_mammal_{species}"
        records.append(record)
        sequences[species_node] = target_sequence_on_human_axis(homology, len(human_sequence))
        leaf_attrs[species_node] = {
            "label": record.get("label", species),
            "is_leaf": True,
            "taxon_count": 1,
            "is_mammal_ortholog": True,
            "is_mammal_context": True,
            "mammal_clade": record.get("clade", "Species ortholog context"),
            "species": species,
            "species_label": record.get("label", species),
            "orthology_type": record.get("orthology_type", ""),
            "taxonomy_level": record.get("taxonomy_level", ""),
            "ensembl_gene_id": record.get("ensembl_gene_id", ""),
            "protein_id": record.get("protein_id", ""),
            "taxon_id": record.get("taxon_id", ""),
            "percent_identity": float(record.get("percent_identity", 0.0)),
            "percent_positive": float(record.get("percent_positive", 0.0)),
            "alignment_coverage": float(record.get("alignment_coverage", 0.0)),
            "layout_weight": MAMMAL_LAYOUT_WEIGHT,
        }

    if len(sequences) > FAST_LINKAGE_TREE_THRESHOLD:
        inferred_tree, inferred_root = build_linkage_distance_graph(gene, leaf_attrs, sequences)
        tree_method = (
            "average-linkage hierarchical clustering on p-distance over shared aligned amino-acid positions"
        )
    else:
        inferred_tree, inferred_root = build_neighbor_joining_graph(gene, leaf_attrs, sequences)
        tree_method = (
            "neighbor joining on p-distance over shared aligned amino-acid positions from homology alignments"
        )
    inferred = orient_unrooted_tree(inferred_tree, inferred_root)
    attach_human_population_summaries(inferred, graph, gene, human_anchor)

    inferred.graph.update(graph.graph)
    inferred.graph["root"] = inferred_root
    inferred.graph["plot_mode"] = "taxon_aligned_inferred_ortholog_tree"
    inferred.graph["species_backbone"] = tree_method
    inferred.graph["tree_inference_method"] = tree_method
    inferred.graph["homology_min_coverage"] = HOMOLOGY_MIN_COVERAGE
    inferred.graph["displayed_taxa_per_gene"] = len(records) + (len(POP_ORDER) * HUMAN_SUBTREE_LEAVES_PER_POPULATION)
    inferred.graph["mammal_orthologs"] = len(records)
    inferred.graph["human_variant_taxa"] = len(POP_ORDER) * HUMAN_SUBTREE_LEAVES_PER_POPULATION
    inferred.graph["supporting_human_variant_leaves"] = human_variant_leaf_count(graph)
    inferred.graph["supporting_mixed_candidates"] = len(mixed_variant_candidates(graph))
    inferred.graph["mammal_taxon_filter"] = (
        f"balanced display subset of same species with usable aligned homology sequences in {', '.join(GENE_ORDER)}"
    )
    inferred.graph["available_common_ortholog_species"] = len(available_common_inferred_homology_species())
    inferred.graph["max_displayed_ortholog_species"] = MAX_DISPLAYED_ORTHOLOG_SPECIES
    inferred.graph["human_taxon_filter"] = "same neutral human subtree leaves expanded from the inferred human sequence leaf"
    inferred.graph["mammal_context_source"] = context.get("source", "")
    inferred.graph["mammal_context_source_url"] = context.get("source_url", "")
    inferred.graph["alignment_source_files"] = {
        item: os.path.relpath(
            NCBI_HOMOLOGY_ALIGNMENT_PATHS[item]
            if os.path.exists(NCBI_HOMOLOGY_ALIGNMENT_PATHS[item])
            else path,
            ROOT,
        )
        for item, path in HOMOLOGY_ALIGNMENT_PATHS.items()
    }
    return inferred


def compose_mammal_context_tree(graph: nx.DiGraph, gene: str, context: dict) -> nx.DiGraph:
    """Attach shared species ortholog leaves to the human variant tree."""
    common_species = common_available_mammal_species(context)
    records = [
        record
        for record in context.get("genes", {}).get(gene, [])
        if record.get("available") and record.get("species") in common_species
    ]
    if not records:
        return graph

    composite = graph.copy()
    original_root = root_node(graph)
    context_root = f"{gene}_mammal_context_root"
    human_anchor = f"{gene}_human_variant_context"

    composite.add_node(
        context_root,
        label=f"{gene} species ortholog context",
        is_leaf=False,
        taxon_count=int(graph.graph.get("selected_taxa", 0)) + len(records),
        is_mammal_context=True,
        mammal_clade="Species ortholog context",
        layout_weight=MAMMAL_LAYOUT_WEIGHT * len(records),
    )
    composite.add_node(
        human_anchor,
        label="Human variation subtree",
        is_leaf=False,
        taxon_count=int(graph.graph.get("selected_taxa", 0)),
        is_human_variant_context=True,
        mammal_clade="Primates",
        layout_weight=graph.graph.get("selected_taxa", 1),
    )
    composite.add_edge(context_root, human_anchor, weight=0.075)
    composite.add_edge(human_anchor, original_root, weight=0.045)

    for clade in MAMMAL_CLADE_ORDER:
        clade_records = [record for record in records if record.get("clade") == clade]
        if not clade_records:
            continue
        clade_node = f"{gene}_mammal_clade_{clade}"
        composite.add_node(
            clade_node,
            label=clade,
            is_leaf=False,
            taxon_count=len(clade_records),
            is_mammal_context=True,
            mammal_clade=clade,
            layout_weight=MAMMAL_LAYOUT_WEIGHT * len(clade_records),
        )
        composite.add_edge(context_root, clade_node, weight=0.16 + (0.025 * MAMMAL_CLADE_ORDER.index(clade)))

        for record in clade_records:
            species = str(record.get("species"))
            species_node = f"{gene}_mammal_{species}"
            percent_identity = float(record.get("percent_identity", 0.0))
            branch_length = max(0.025, (100.0 - percent_identity) / 165.0)
            composite.add_node(
                species_node,
                label=record.get("label", species),
                is_leaf=True,
                taxon_count=1,
                is_mammal_ortholog=True,
                is_mammal_context=True,
                mammal_clade=clade,
                species=species,
                species_label=record.get("label", species),
                orthology_type=record.get("orthology_type", ""),
                taxonomy_level=record.get("taxonomy_level", ""),
                ensembl_gene_id=record.get("ensembl_gene_id", ""),
                protein_id=record.get("protein_id", ""),
                percent_identity=percent_identity,
                percent_positive=float(record.get("percent_positive", 0.0)),
                layout_weight=MAMMAL_LAYOUT_WEIGHT,
            )
            composite.add_edge(clade_node, species_node, weight=branch_length)

    composite.graph.update(graph.graph)
    composite.graph["root"] = context_root
    composite.graph["mammal_orthologs"] = len(records)
    composite.graph["mammal_taxon_filter"] = "common available species orthologs across all plotted genes"
    composite.graph["mammal_context_source"] = context.get("source", "")
    composite.graph["mammal_context_source_url"] = context.get("source_url", "")
    return composite


def population_summary_stats(graph: nx.DiGraph, population: str) -> dict:
    leaves = [
        attrs
        for _, attrs in graph.nodes(data=True)
        if attrs.get("is_leaf")
        and not attrs.get("is_mammal_ortholog")
        and attrs.get("assigned_population") == population
    ]
    af_values = [float(attrs.get("assigned_af", 0.0)) for attrs in leaves]
    mixed_count = sum(1 for item in mixed_variant_candidates(graph) if item["dominant"] == population)
    return {
        "count": len(leaves),
        "mean_af": float(np.mean(af_values)) if af_values else 0.0,
        "max_af": float(max(af_values)) if af_values else 0.0,
        "mixed_count": mixed_count,
    }


def compose_taxon_aligned_teaching_tree(graph: nx.DiGraph, gene: str, context: dict) -> nx.DiGraph:
    """Build a same-taxa teaching tree on a bifurcating species backbone."""
    records = [
        record
        for record in context.get("genes", {}).get(gene, [])
        if record.get("available") and record.get("species") in common_available_mammal_species(context)
    ]
    teaching = nx.DiGraph()
    root = f"{gene}_taxon_aligned_root"

    teaching.add_node(
        root,
        label=f"{gene} fully resolved same-taxa comparison",
        is_leaf=False,
        taxon_count=len(records) + len(POP_ORDER),
        is_species_backbone=True,
    )

    def safe_id(label: str) -> str:
        return (
            label.lower()
            .replace(" ", "_")
            .replace("-", "_")
            .replace("/", "_")
            .replace("(", "")
            .replace(")", "")
        )

    def add_internal(label: str, clade: str = "Species ortholog context") -> str:
        node = f"{gene}_backbone_{safe_id(label)}"
        if node in teaching:
            return node
        teaching.add_node(
            node,
            label=label,
            is_leaf=False,
            is_mammal_context=True,
            is_species_backbone=True,
            mammal_clade=clade,
        )
        return node

    def connect(parent: str, child: str, weight: float) -> None:
        if not teaching.has_edge(parent, child):
            teaching.add_edge(parent, child, weight=weight)

    def add_species_leaf(parent: str, record: dict, edge_weight: float | None = None) -> str:
        clade = str(record.get("clade", "Species ortholog context"))
        species = str(record.get("species"))
        species_node = f"{gene}_mammal_{species}"
        percent_identity = float(record.get("percent_identity", 0.0))
        branch_length = edge_weight if edge_weight is not None else max(0.025, (100.0 - percent_identity) / 165.0)
        teaching.add_node(
            species_node,
            label=record.get("label", species),
            is_leaf=True,
            taxon_count=1,
            is_mammal_ortholog=True,
            is_mammal_context=True,
            mammal_clade=clade,
            species=species,
            species_label=record.get("label", species),
            orthology_type=record.get("orthology_type", ""),
            taxonomy_level=record.get("taxonomy_level", ""),
            ensembl_gene_id=record.get("ensembl_gene_id", ""),
            protein_id=record.get("protein_id", ""),
            percent_identity=percent_identity,
            percent_positive=float(record.get("percent_positive", 0.0)),
            layout_weight=MAMMAL_LAYOUT_WEIGHT,
        )
        connect(parent, species_node, branch_length)
        return species_node

    def add_population_leaf(parent: str, population: str) -> str:
        stats = population_summary_stats(graph, population)
        population_node = f"{gene}_human_population_{population}"
        count = max(1, stats["count"])
        teaching.add_node(
            population_node,
            label=f"Human variation leaf {POP_ORDER.index(population) + 1}",
            is_leaf=True,
            taxon_count=1,
            is_population_root=True,
            is_human_population_summary=True,
            assigned_population=population,
            population_label=population,
            supporting_variant_count=stats["count"],
            mean_assigned_af=stats["mean_af"],
            max_assigned_af=stats["max_af"],
            supporting_mixed_candidates=stats["mixed_count"],
            layout_weight=count,
        )
        branch_length = 0.10 + (0.18 * (1.0 - min(stats["mean_af"] * 200.0, 1.0)))
        connect(parent, population_node, branch_length)
        return population_node

    def split_records(items: list[dict], label: str, clade: str) -> tuple[list[dict], list[dict]]:
        ordered = sorted(
            items,
            key=lambda item: (
                -float(item.get("percent_identity", 0.0)),
                str(item.get("label", item.get("species", ""))),
            ),
        )
        midpoint = max(1, len(ordered) // 2)
        return ordered[:midpoint], ordered[midpoint:]

    def attach_resolved_records(parent: str, clade: str, items: list[dict], depth: int = 0) -> None:
        if not items:
            return
        if len(items) == 1:
            add_species_leaf(parent, items[0])
            return
        node = add_internal(f"{clade} split {depth + 1}", clade)
        connect(parent, node, 0.12 + (0.025 * depth))
        left, right = split_records(items, node, clade)
        attach_resolved_records(node, clade, left, depth + 1)
        attach_resolved_records(node, clade, right, depth + 1)

    def attach_resolved_populations(parent: str, populations: list[str], depth: int = 0) -> None:
        if len(populations) == 1:
            add_population_leaf(parent, populations[0])
            return
        node = add_internal(f"Human variation split {depth + 1}", "Primates")
        connect(parent, node, 0.12 + (0.02 * depth))
        midpoint = len(populations) // 2
        attach_resolved_populations(node, populations[:midpoint], depth + 1)
        attach_resolved_populations(node, populations[midpoint:], depth + 1)

    records_by_clade: dict[str, list[dict]] = {
        clade: [record for record in records if record.get("clade") == clade]
        for clade in MAMMAL_CLADE_ORDER
    }

    def attach_clade(parent: str, clade: str, edge_weight: float = 0.16) -> None:
        items = records_by_clade.get(clade, [])
        if not items:
            return
        node = add_internal(clade, clade)
        connect(parent, node, edge_weight)
        attach_resolved_records(node, clade, items)

    # Bifurcating species backbone. Terminal clades are then internally resolved
    # into deterministic binary subtrees, so no clade is rendered as a star.
    jawless = add_internal("Jawless fish", "Jawless fish")
    gnathostomes = add_internal("Jawed vertebrates")
    connect(root, jawless, 0.24)
    connect(root, gnathostomes, 0.18)
    attach_clade(jawless, "Jawless fish")

    ray_finned = add_internal("Ray-finned fish", "Ray-finned fish")
    sarcopterygians = add_internal("Sarcopterygians")
    connect(gnathostomes, ray_finned, 0.18)
    connect(gnathostomes, sarcopterygians, 0.18)
    attach_resolved_records(ray_finned, "Ray-finned fish", records_by_clade.get("Ray-finned fish", []))

    lobe_finned = add_internal("Lobe-finned fish", "Lobe-finned fish")
    tetrapods = add_internal("Tetrapods")
    connect(sarcopterygians, lobe_finned, 0.18)
    connect(sarcopterygians, tetrapods, 0.18)
    attach_clade(lobe_finned, "Lobe-finned fish")

    amphibians = add_internal("Amphibians", "Amphibians")
    amniotes = add_internal("Amniotes")
    connect(tetrapods, amphibians, 0.18)
    connect(tetrapods, amniotes, 0.18)
    attach_resolved_records(amphibians, "Amphibians", records_by_clade.get("Amphibians", []))

    sauropsids = add_internal("Sauropsids")
    mammals = add_internal("Mammals")
    connect(amniotes, sauropsids, 0.16)
    connect(amniotes, mammals, 0.16)

    reptiles = add_internal("Reptiles", "Reptiles")
    birds = add_internal("Birds", "Birds")
    connect(sauropsids, reptiles, 0.16)
    connect(sauropsids, birds, 0.16)
    attach_resolved_records(reptiles, "Reptiles", records_by_clade.get("Reptiles", []))
    attach_resolved_records(birds, "Birds", records_by_clade.get("Birds", []))

    monotremes = add_internal("Monotremata", "Monotremata")
    therians = add_internal("Therian mammals")
    connect(mammals, monotremes, 0.15)
    connect(mammals, therians, 0.15)
    attach_resolved_records(monotremes, "Monotremata", records_by_clade.get("Monotremata", []))

    marsupials = add_internal("Marsupialia", "Marsupialia")
    placentals = add_internal("Placental mammals")
    connect(therians, marsupials, 0.15)
    connect(therians, placentals, 0.15)
    attach_resolved_records(marsupials, "Marsupialia", records_by_clade.get("Marsupialia", []))

    atlantogenata = add_internal("Atlantogenata")
    boreoeutheria = add_internal("Boreoeutheria")
    connect(placentals, atlantogenata, 0.14)
    connect(placentals, boreoeutheria, 0.14)

    afrotheria = add_internal("Afrotheria", "Afrotheria")
    xenarthra = add_internal("Xenarthra", "Xenarthra")
    connect(atlantogenata, afrotheria, 0.14)
    connect(atlantogenata, xenarthra, 0.14)
    attach_resolved_records(afrotheria, "Afrotheria", records_by_clade.get("Afrotheria", []))
    attach_resolved_records(xenarthra, "Xenarthra", records_by_clade.get("Xenarthra", []))

    euarchontoglires = add_internal("Euarchontoglires")
    laurasiatheria = add_internal("Laurasiatheria")
    connect(boreoeutheria, euarchontoglires, 0.14)
    connect(boreoeutheria, laurasiatheria, 0.14)

    primatomorpha = add_internal("Primatomorpha")
    glires = add_internal("Glires", "Glires")
    connect(euarchontoglires, primatomorpha, 0.13)
    connect(euarchontoglires, glires, 0.13)
    attach_resolved_records(glires, "Glires", records_by_clade.get("Glires", []))

    scandentia = add_internal("Scandentia", "Scandentia")
    primates = add_internal("Primates")
    connect(primatomorpha, scandentia, 0.13)
    connect(primatomorpha, primates, 0.13)
    attach_clade(scandentia, "Scandentia")

    strepsirrhines = add_internal("Strepsirrhines", "Strepsirrhines")
    haplorrhines = add_internal("Haplorrhines")
    connect(primates, strepsirrhines, 0.13)
    connect(primates, haplorrhines, 0.13)
    attach_resolved_records(strepsirrhines, "Strepsirrhines", records_by_clade.get("Strepsirrhines", []))

    tarsiers = add_internal("Tarsiiformes", "Tarsiiformes")
    simians = add_internal("Simians")
    connect(haplorrhines, tarsiers, 0.12)
    connect(haplorrhines, simians, 0.12)
    attach_clade(tarsiers, "Tarsiiformes")

    new_world = add_internal("New World monkeys", "New World monkeys")
    catarrhini = add_internal("Catarrhini")
    connect(simians, new_world, 0.12)
    connect(simians, catarrhini, 0.12)
    attach_resolved_records(new_world, "New World monkeys", records_by_clade.get("New World monkeys", []))

    old_world = add_internal("Old World monkeys", "Old World monkeys")
    hominoids = add_internal("Hominoids")
    connect(catarrhini, old_world, 0.12)
    connect(catarrhini, hominoids, 0.12)
    attach_resolved_records(old_world, "Old World monkeys", records_by_clade.get("Old World monkeys", []))

    lesser_apes = add_internal("Lesser apes", "Lesser apes")
    great_ape_human = add_internal("Great apes + human subtree", "Great apes")
    connect(hominoids, lesser_apes, 0.11)
    connect(hominoids, great_ape_human, 0.11)
    attach_clade(lesser_apes, "Lesser apes")

    great_apes = add_internal("Great apes", "Great apes")
    human_summaries = add_internal("Human variation subtree", "Primates")
    connect(great_ape_human, great_apes, 0.11)
    connect(great_ape_human, human_summaries, 0.11)
    attach_resolved_records(great_apes, "Great apes", records_by_clade.get("Great apes", []))
    attach_resolved_populations(human_summaries, POP_ORDER)

    eulipotyphla = add_internal("Eulipotyphla", "Eulipotyphla")
    scrotifera = add_internal("Scrotifera")
    connect(laurasiatheria, eulipotyphla, 0.13)
    connect(laurasiatheria, scrotifera, 0.13)
    attach_resolved_records(eulipotyphla, "Eulipotyphla", records_by_clade.get("Eulipotyphla", []))

    bats = add_internal("Bats", "Bats")
    ferungulata = add_internal("Ferungulata")
    connect(scrotifera, bats, 0.13)
    connect(scrotifera, ferungulata, 0.13)
    attach_resolved_records(bats, "Bats", records_by_clade.get("Bats", []))

    carnivores = add_internal("Carnivores", "Carnivores")
    ungulates = add_internal("Ungulates")
    connect(ferungulata, carnivores, 0.12)
    connect(ferungulata, ungulates, 0.12)
    attach_resolved_records(carnivores, "Carnivores", records_by_clade.get("Carnivores", []))

    perissodactyla = add_internal("Perissodactyla", "Perissodactyla")
    cetartiodactyla = add_internal("Cetartiodactyla", "Cetartiodactyla")
    connect(ungulates, perissodactyla, 0.12)
    connect(ungulates, cetartiodactyla, 0.12)
    attach_resolved_records(perissodactyla, "Perissodactyla", records_by_clade.get("Perissodactyla", []))
    attach_resolved_records(cetartiodactyla, "Cetartiodactyla", records_by_clade.get("Cetartiodactyla", []))

    teaching.graph.update(graph.graph)
    teaching.graph["root"] = root
    teaching.graph["plot_mode"] = "taxon_aligned_teaching"
    teaching.graph["species_backbone"] = "bifurcating vertebrate taxonomy backbone with deterministic binary terminal clade resolution"
    teaching.graph["displayed_taxa_per_gene"] = len(records) + len(POP_ORDER)
    teaching.graph["mammal_orthologs"] = len(records)
    teaching.graph["human_population_taxa"] = len(POP_ORDER)
    teaching.graph["supporting_human_variant_leaves"] = human_variant_leaf_count(graph)
    teaching.graph["supporting_mixed_candidates"] = len(mixed_variant_candidates(graph))
    teaching.graph["mammal_taxon_filter"] = "common available species orthologs across all plotted genes"
    teaching.graph["human_taxon_filter"] = "same neutral human subtree leaves across all plotted genes"
    teaching.graph["mammal_context_source"] = context.get("source", "")
    teaching.graph["mammal_context_source_url"] = context.get("source_url", "")
    return teaching


def population_angle(population: str) -> float:
    start_angle = -math.pi / 2.0
    return start_angle + (POP_ORDER.index(population) * math.tau / len(POP_ORDER))


def square_boundary_point(angle: float, half_size: float = BOX_HALF_SIZE) -> np.ndarray:
    x_direction = math.cos(angle)
    y_direction = math.sin(angle)
    scale = half_size / max(abs(x_direction), abs(y_direction), 1e-9)
    return np.array([x_direction * scale, y_direction * scale, 0.0], dtype=float)


def population_anchors(geometry: str = "box") -> dict[str, np.ndarray]:
    anchors = {}
    for population in POP_ORDER:
        angle = population_angle(population)
        if geometry == "hyperbolic":
            radius = HYPERBOLIC_RADIUS * 0.78
            anchors[population] = np.array(
                [radius * math.cos(angle), radius * math.sin(angle), 0.0],
                dtype=float,
            )
        else:
            anchors[population] = square_boundary_point(angle, BOX_HALF_SIZE * BOX_FILL_FRACTION)
    return anchors


def root_node(graph: nx.DiGraph) -> str:
    return str(graph.graph.get("root") or next(node for node in graph if graph.in_degree(node) == 0))


def population_roots(graph: nx.DiGraph) -> dict[str, str]:
    roots = {}
    for node, attrs in graph.nodes(data=True):
        if attrs.get("is_population_root"):
            roots[str(attrs["population_label"])] = node
    return roots


def leaf_afs(attrs: dict) -> dict[str, float]:
    return {population: float(attrs.get(key, 0.0)) for population, key in zip(POP_ORDER, AF_KEYS)}


def mixed_variant_candidates(graph: nx.DiGraph) -> list[dict]:
    candidates = []
    for node, attrs in graph.nodes(data=True):
        if not attrs.get("is_leaf"):
            continue
        values = leaf_afs(attrs)
        dominant = str(attrs.get("assigned_population"))
        dominant_af = values.get(dominant, 0.0)
        if dominant_af <= 0:
            continue
        secondary = max((pop for pop in POP_ORDER if pop != dominant), key=lambda pop: values[pop])
        secondary_af = values[secondary]
        ratio = secondary_af / dominant_af
        if ratio < MIX_RATIO_THRESHOLD or secondary_af < MIX_MIN_SECONDARY_AF:
            continue
        candidates.append(
            {
                "node": node,
                "dominant": dominant,
                "secondary": secondary,
                "dominant_af": dominant_af,
                "secondary_af": secondary_af,
                "ratio": ratio,
                "variant_id": attrs.get("variant_id", attrs.get("label", node)),
            }
        )
    return sorted(candidates, key=lambda item: (-item["ratio"], item["dominant"], item["secondary"]))


def displayed_human_leaf_nodes(graph: nx.DiGraph) -> set[str]:
    selected = {item["node"] for item in mixed_variant_candidates(graph)}
    for population in POP_ORDER:
        leaves = []
        for node, attrs in graph.nodes(data=True):
            if not attrs.get("is_leaf") or attrs.get("is_mammal_ortholog"):
                continue
            if attrs.get("assigned_population") != population:
                continue
            leaves.append(
                (
                    float(attrs.get("assigned_af", 0.0)),
                    str(attrs.get("variant_id", attrs.get("label", node))),
                    node,
                )
            )
        leaves.sort(key=lambda item: (-item[0], item[1]))
        selected.update(node for _, _, node in leaves[:HUMAN_MARKERS_PER_POPULATION])
    return selected


def layout_graph(graph: nx.DiGraph) -> nx.Graph:
    undirected = nx.Graph()
    undirected.add_nodes_from(graph.nodes(data=True))
    lengths = [float(data.get("weight", 0.0)) for _, _, data in graph.edges(data=True)]
    positive = [value for value in lengths if value > 0]
    floor = max(min(positive) if positive else 1.0, 1e-5)
    for parent, child, data in graph.edges(data=True):
        length = max(float(data.get("weight", 0.0)), floor)
        undirected.add_edge(
            parent,
            child,
            layout_length=length,
        )
    return undirected


def leaf_count_in_component(tree: nx.Graph, node: str, parent: str | None, cache: dict[tuple[str, str | None], int]) -> int:
    key = (node, parent)
    if key in cache:
        return cache[key]
    attrs = tree.nodes[node]
    if attrs.get("is_leaf"):
        cache[key] = max(1, int(attrs.get("layout_weight", 1)))
        return cache[key]
    total = 0
    for neighbor in tree.neighbors(node):
        if neighbor == parent:
            continue
        total += leaf_count_in_component(tree, neighbor, node, cache)
    cache[key] = max(1, total)
    return cache[key]


def component_population_score(tree: nx.Graph, node: str, parent: str | None) -> float:
    scores = []
    stack = [(node, parent)]
    seen = set()
    while stack:
        current, previous = stack.pop()
        if current in seen:
            continue
        seen.add(current)
        attrs = tree.nodes[current]
        population = attrs.get("assigned_population") or attrs.get("population_label")
        if population in POP_ORDER:
            scores.append(float(POP_ORDER.index(population)))
        for neighbor in tree.neighbors(current):
            if neighbor != previous:
                stack.append((neighbor, current))
    return sum(scores) / len(scores) if scores else len(POP_ORDER)


def choose_layout_center(tree: nx.Graph) -> str:
    centers = nx.center(tree)
    return max(centers, key=lambda node: (tree.degree(node), str(node)))


def transformed_edge_length(tree: nx.Graph, left: str, right: str, max_length: float) -> float:
    raw = float(tree.edges[left, right].get("layout_length", 0.0))
    scaled = math.sqrt(raw / max_length) if max_length > 0 else 0.0
    is_terminal = bool(tree.nodes[left].get("is_leaf") or tree.nodes[right].get("is_leaf"))
    if is_terminal:
        return 24.0 + (70.0 * scaled)
    return 12.0 + (48.0 * scaled)


def assign_equal_angle_layout(
    tree: nx.Graph,
    node: str,
    parent: str | None,
    theta_min: float,
    theta_max: float,
    positions: dict[str, np.ndarray],
    count_cache: dict[tuple[str, str | None], int],
    max_length: float,
) -> None:
    children = [neighbor for neighbor in tree.neighbors(node) if neighbor != parent]
    if not children:
        return

    children.sort(
        key=lambda child: (
            component_population_score(tree, child, node),
            -leaf_count_in_component(tree, child, node, count_cache),
            str(child),
        )
    )
    total = sum(leaf_count_in_component(tree, child, node, count_cache) for child in children)
    cursor = theta_min
    span = theta_max - theta_min
    for child in children:
        child_count = leaf_count_in_component(tree, child, node, count_cache)
        child_span = span * (child_count / total) if total else span / len(children)
        child_min = cursor
        child_max = cursor + child_span
        child_angle = (child_min + child_max) / 2.0
        length = transformed_edge_length(tree, node, child, max_length)
        positions[child] = positions[node] + np.array(
            [length * math.cos(child_angle), length * math.sin(child_angle)],
            dtype=float,
        )
        assign_equal_angle_layout(tree, child, node, child_min, child_max, positions, count_cache, max_length)
        cursor = child_max


def align_layout_to_population_anchors(
    graph: nx.DiGraph,
    positions_2d: dict[str, np.ndarray],
    geometry: str = "box",
) -> dict[str, np.ndarray]:
    roots = population_roots(graph)
    anchors = population_anchors(geometry)
    available = [population for population in POP_ORDER if population in roots]
    source = np.array([positions_2d[roots[population]] for population in available], dtype=float)
    target = np.array([anchors[population][:2] for population in available], dtype=float)
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    source_centered = source - source_center
    target_centered = target - target_center

    best_transform = None
    best_error = math.inf
    for reflection in (np.eye(2), np.array([[-1.0, 0.0], [0.0, 1.0]])):
        reflected = source_centered @ reflection
        covariance = reflected.T @ target_centered
        left, _singular, right_t = np.linalg.svd(covariance)
        rotation = left @ right_t
        scale = np.linalg.norm(target_centered) / max(np.linalg.norm(reflected), 1e-9)
        transformed = reflected @ rotation * scale
        error = float(np.mean(np.linalg.norm(transformed - target_centered, axis=1)))
        if error < best_error:
            best_error = error
            best_transform = (reflection, rotation, scale)

    reflection, rotation, scale = best_transform
    aligned = {}
    for node, point in positions_2d.items():
        centered = point - source_center
        aligned[node] = ((centered @ reflection) @ rotation) * scale + target_center
    return aligned


def normalize_layout_radius(positions_2d: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    matrix = np.array(list(positions_2d.values()), dtype=float)
    center = matrix.mean(axis=0)
    centered = {node: point - center for node, point in positions_2d.items()}
    max_radius = max(float(np.linalg.norm(point)) for point in centered.values()) or 1.0
    scale = TREE_TARGET_RADIUS / max_radius
    return {node: point * scale for node, point in centered.items()}


def fit_layout_to_box(positions_2d: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Use the full cuboid footprint without changing edge-crossing topology."""
    matrix = np.array(list(positions_2d.values()), dtype=float)
    low = matrix.min(axis=0)
    high = matrix.max(axis=0)
    center = (low + high) / 2.0
    span = np.maximum(high - low, 1e-9)
    target_span = 2.0 * BOX_HALF_SIZE * BOX_FILL_FRACTION
    scale = np.array([target_span / span[0], target_span / span[1]], dtype=float)
    return {node: (point - center) * scale for node, point in positions_2d.items()}


def project_layout_to_hyperbolic_disk(positions_2d: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Fit the equal-angle tree into a Poincare-style disk without bending topology."""
    matrix = np.array(list(positions_2d.values()), dtype=float)
    center = matrix.mean(axis=0)
    centered = {node: point - center for node, point in positions_2d.items()}
    max_radius = max(float(np.linalg.norm(point)) for point in centered.values()) or 1.0
    scale = (HYPERBOLIC_RADIUS * HYPERBOLIC_FILL_FRACTION) / max_radius
    return {node: point * scale for node, point in centered.items()}


def is_inferred_ortholog_tree(graph: nx.DiGraph) -> bool:
    return graph.graph.get("plot_mode") == "taxon_aligned_inferred_ortholog_tree"


def rotate_layout_to_target_vector(
    positions_2d: dict[str, np.ndarray],
    current_vector: np.ndarray,
    target_angle: float,
) -> dict[str, np.ndarray]:
    current_angle = math.atan2(float(current_vector[1]), float(current_vector[0]))
    angle = target_angle - current_angle
    rotation = np.array(
        [
            [math.cos(angle), -math.sin(angle)],
            [math.sin(angle), math.cos(angle)],
        ],
        dtype=float,
    )
    return {node: point @ rotation.T for node, point in positions_2d.items()}


def orient_inferred_ortholog_layout(
    graph: nx.DiGraph,
    positions_2d: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Orient inferred trees so human/great apes sit opposite the outgroups."""
    human_nodes = [
        node
        for node, attrs in graph.nodes(data=True)
        if attrs.get("is_human_sequence_anchor")
    ]
    if not human_nodes:
        return positions_2d
    human_point = positions_2d[human_nodes[0]]
    outgroup_points = [
        positions_2d[node]
        for node, attrs in graph.nodes(data=True)
        if attrs.get("is_mammal_ortholog")
        and attrs.get("mammal_clade") in {"Fungi outgroup", "Jawless fish", "Lobe-finned fish"}
        and node in positions_2d
    ]
    if not outgroup_points:
        return positions_2d
    outgroup_centroid = np.mean(np.array(outgroup_points, dtype=float), axis=0)
    vector = human_point - outgroup_centroid
    if float(np.linalg.norm(vector)) < 1e-9:
        return positions_2d
    return rotate_layout_to_target_vector(positions_2d, vector, math.radians(35.0))


def branch_depth_offsets(tree: nx.Graph, center: str, max_length: float) -> dict[str, float]:
    depth_graph = nx.Graph()
    depth_graph.add_nodes_from(tree.nodes)
    for left, right in tree.edges:
        depth_graph.add_edge(left, right, depth_length=transformed_edge_length(tree, left, right, max_length))

    depths = nx.single_source_dijkstra_path_length(depth_graph, center, weight="depth_length")
    max_depth = max(depths.values()) if depths else 1.0
    if max_depth <= 0:
        return {node: 0.0 for node in tree.nodes}
    return {node: ((depths.get(node, 0.0) / max_depth) - 0.5) * LOCAL_Z_SPAN for node in tree.nodes}


def hyperbolic_terminal_leaf_z_offsets(
    graph: nx.DiGraph,
    positions_2d: dict[str, np.ndarray],
) -> dict[str, float]:
    leaves = []
    displayed_human_nodes = displayed_human_leaf_nodes(graph)
    for node, attrs in graph.nodes(data=True):
        if node not in positions_2d or not attrs.get("is_leaf"):
            continue
        if not attrs.get("is_mammal_ortholog") and node not in displayed_human_nodes:
            continue
        point = positions_2d[node]
        leaves.append((math.atan2(float(point[1]), float(point[0])), float(np.linalg.norm(point)), str(node)))
    leaves.sort()
    pattern = [-2.0, 0.0, 2.0, -1.0, 1.0]
    return {
        node: pattern[index % len(pattern)] * HYPERBOLIC_LEAF_Z_STAGGER
        for index, (_angle, _radius, node) in enumerate(leaves)
    }


def descendant_leaf_supports(graph: nx.DiGraph) -> dict[str, int]:
    support: dict[str, int] = {}
    for node in reversed(list(nx.topological_sort(graph))):
        node = str(node)
        children = [str(child) for child in graph.successors(node)]
        attrs = graph.nodes[node]
        if attrs.get("is_leaf") or not children:
            support[node] = 1
            continue
        support[node] = max(1, sum(support.get(child, 1) for child in children))
    return support


def internal_radius_multiplier(support: int, max_support: int) -> float:
    if max_support <= 1:
        return 0.72
    normalized = math.log1p(max(1, support)) / math.log1p(max_support)
    return 0.58 + (0.92 * (normalized**0.72))


def compute_unrooted_3d_layout(graph: nx.DiGraph, gene: str, geometry: str = "box") -> dict[str, np.ndarray]:
    """Equal-angle unrooted tree layout fitted into a stacked 3D volume."""
    tree = layout_graph(graph)
    center = choose_layout_center(tree)
    count_cache: dict[tuple[str, str | None], int] = {}
    max_length = max((float(data.get("layout_length", 0.0)) for _, _, data in tree.edges(data=True)), default=1.0)
    positions_2d = {center: np.array([0.0, 0.0], dtype=float)}
    assign_equal_angle_layout(
        tree=tree,
        node=center,
        parent=None,
        theta_min=-math.pi,
        theta_max=math.pi,
        positions=positions_2d,
        count_cache=count_cache,
        max_length=max_length,
    )
    positions_2d = normalize_layout_radius(positions_2d)
    if is_inferred_ortholog_tree(graph):
        positions_2d = orient_inferred_ortholog_layout(graph, positions_2d)
        positions_2d = normalize_layout_radius(positions_2d)
    else:
        positions_2d = align_layout_to_population_anchors(graph, positions_2d, geometry)
        positions_2d = normalize_layout_radius(positions_2d)
    if geometry == "hyperbolic":
        positions_2d = project_layout_to_hyperbolic_disk(positions_2d)
    else:
        positions_2d = fit_layout_to_box(positions_2d)
    z_offsets = branch_depth_offsets(tree, center, max_length)
    leaf_z_offsets = hyperbolic_terminal_leaf_z_offsets(graph, positions_2d) if geometry == "hyperbolic" else {}

    return {
        node: np.array([point[0], point[1], GENE_Z[gene] + z_offsets[node] + leaf_z_offsets.get(node, 0.0)], dtype=float)
        for node, point in positions_2d.items()
    }


def plane_trace(gene: str) -> go.Surface:
    z = GENE_Z[gene]
    return go.Surface(
        x=[[-BOX_HALF_SIZE, BOX_HALF_SIZE], [-BOX_HALF_SIZE, BOX_HALF_SIZE]],
        y=[[-BOX_HALF_SIZE, -BOX_HALF_SIZE], [BOX_HALF_SIZE, BOX_HALF_SIZE]],
        z=[[z, z], [z, z]],
        opacity=0.08,
        showscale=False,
        hoverinfo="skip",
        surfacecolor=[[0, 0], [0, 0]],
        colorscale=[[0, "rgba(160,190,235,0.06)"], [1, "rgba(160,190,235,0.06)"]],
        showlegend=False,
        name=f"{gene} reference plane",
    )


def box_frame_trace(gene: str) -> go.Scatter3d:
    z_low = GENE_Z[gene] - (LOCAL_Z_SPAN / 2.0)
    z_high = GENE_Z[gene] + (LOCAL_Z_SPAN / 2.0)
    h = BOX_HALF_SIZE
    corners = {
        "a": np.array([-h, -h, z_low]),
        "b": np.array([h, -h, z_low]),
        "c": np.array([h, h, z_low]),
        "d": np.array([-h, h, z_low]),
        "e": np.array([-h, -h, z_high]),
        "f": np.array([h, -h, z_high]),
        "g": np.array([h, h, z_high]),
        "i": np.array([-h, h, z_high]),
    }
    segments = [
        ("a", "b"),
        ("b", "c"),
        ("c", "d"),
        ("d", "a"),
        ("e", "f"),
        ("f", "g"),
        ("g", "i"),
        ("i", "e"),
        ("a", "e"),
        ("b", "f"),
        ("c", "g"),
        ("d", "i"),
    ]
    x_values, y_values, z_values = [], [], []
    for start_key, end_key in segments:
        start = corners[start_key]
        end = corners[end_key]
        x_values.extend([start[0], end[0], None])
        y_values.extend([start[1], end[1], None])
        z_values.extend([start[2], end[2], None])
    return go.Scatter3d(
        x=x_values,
        y=y_values,
        z=z_values,
        mode="lines",
        line=dict(color="rgba(185,210,242,0.24)", width=1.15),
        hoverinfo="skip",
        showlegend=False,
        name=f"{gene} cuboid frame",
    )


def hyperbolic_disk_trace(gene: str) -> go.Surface:
    z = GENE_Z[gene]
    radial_steps = np.linspace(0.0, HYPERBOLIC_RADIUS, 28)
    theta_steps = np.linspace(0.0, math.tau, 96)
    x_values = []
    y_values = []
    z_values = []
    surface = []
    for radius in radial_steps:
        x_values.append([radius * math.cos(theta) for theta in theta_steps])
        y_values.append([radius * math.sin(theta) for theta in theta_steps])
        z_values.append([z for _ in theta_steps])
        surface.append([radius / HYPERBOLIC_RADIUS for _ in theta_steps])
    return go.Surface(
        x=x_values,
        y=y_values,
        z=z_values,
        surfacecolor=surface,
        colorscale=[
            [0.0, "rgba(21,42,70,0.18)"],
            [0.72, "rgba(42,76,112,0.10)"],
            [1.0, "rgba(132,184,235,0.20)"],
        ],
        opacity=0.18,
        showscale=False,
        hoverinfo="skip",
        showlegend=False,
        name=f"{gene} hyperbolic disk",
    )


def hyperbolic_grid_traces(gene: str) -> list[go.Scatter3d]:
    z = GENE_Z[gene] - 2.0
    traces = []
    theta_values = np.linspace(0.0, math.tau, 144)
    for radius_fraction in (0.25, 0.50, 0.72, 0.88, 0.982):
        radius = HYPERBOLIC_RADIUS * radius_fraction
        traces.append(
            go.Scatter3d(
                x=[radius * math.cos(theta) for theta in theta_values],
                y=[radius * math.sin(theta) for theta in theta_values],
                z=[z for _ in theta_values],
                mode="lines",
                line=dict(color="rgba(184,213,245,0.16)", width=0.8 if radius_fraction < 0.98 else 1.6),
                hoverinfo="skip",
                showlegend=False,
                name=f"{gene} hyperbolic radius {radius_fraction:.2f}",
            )
        )

    for angle in np.linspace(0.0, math.tau, 16, endpoint=False):
        traces.append(
            go.Scatter3d(
                x=[0.0, HYPERBOLIC_RADIUS * math.cos(angle)],
                y=[0.0, HYPERBOLIC_RADIUS * math.sin(angle)],
                z=[z, z],
                mode="lines",
                line=dict(color="rgba(184,213,245,0.10)", width=0.7),
                hoverinfo="skip",
                showlegend=False,
                name=f"{gene} hyperbolic spoke",
            )
        )
    return traces


def gene_label_trace(gene: str) -> go.Scatter3d:
    return go.Scatter3d(
        x=[-BOX_HALF_SIZE * 1.06],
        y=[-BOX_HALF_SIZE * 1.06],
        z=[GENE_Z[gene] + (LOCAL_Z_SPAN / 2.0) + 42.0],
        mode="text",
        text=[f"<b>{gene}</b>"],
        textfont=dict(size=24, color="#F4F8FF", family="Avenir Next, IBM Plex Sans, sans-serif"),
        hoverinfo="skip",
        showlegend=False,
        name=f"{gene} label",
    )


def edge_traces(graph: nx.DiGraph, layout: dict[str, np.ndarray], gene: str) -> list[go.Scatter3d]:
    buckets: dict[str, dict[str, list[float | None]]] = defaultdict(lambda: {"x": [], "y": [], "z": []})
    for parent, child in graph.edges():
        attrs = graph.nodes[child]
        if (
            attrs.get("is_human_sequence_anchor")
            or attrs.get("is_human_variant_context")
            or attrs.get("is_human_population_summary")
        ):
            population = "Human subtree"
        elif attrs.get("is_mammal_context"):
            population = str(attrs.get("mammal_clade") or "Species ortholog context")
        else:
            population = str(attrs.get("population_label") or attrs.get("assigned_population") or "topology")
        start = layout[parent]
        end = layout[child]
        buckets[population]["x"].extend([start[0], end[0], None])
        buckets[population]["y"].extend([start[1], end[1], None])
        buckets[population]["z"].extend([start[2], end[2], None])

    traces = []
    for population, coords in buckets.items():
        if population in SPECIAL_GROUP_COLORS:
            color = rgba_from_hex(SPECIAL_GROUP_COLORS[population], 0.58)
        elif population in MAMMAL_CLADE_COLORS:
            color = rgba_from_hex(MAMMAL_CLADE_COLORS[population], 0.70)
        else:
            color = "rgba(226,236,250,0.34)"
        traces.append(
            go.Scatter3d(
                x=coords["x"],
                y=coords["y"],
                z=coords["z"],
                mode="lines",
                line=dict(color=color, width=0.82),
                hoverinfo="skip",
                showlegend=False,
                name=f"{gene} {population} tree edges",
            )
        )
    return traces


def node_traces(graph: nx.DiGraph, layout: dict[str, np.ndarray], gene: str) -> list[go.Scatter3d]:
    traces = []
    displayed_human_nodes = displayed_human_leaf_nodes(graph)
    supports = descendant_leaf_supports(graph)
    max_support = max(supports.values(), default=1)

    internal_x, internal_y, internal_z, internal_size, internal_color, internal_hover = [], [], [], [], [], []
    for node, attrs in graph.nodes(data=True):
        if attrs.get("is_leaf") or node not in layout:
            continue
        point = layout[node]
        support = supports.get(str(node), 1)
        radius = internal_radius_multiplier(support, max_support)
        clade = str(attrs.get("mammal_clade") or "")
        if attrs.get("is_human_sequence_anchor") or attrs.get("is_human_variant_context"):
            color = rgba_from_hex(HUMAN_SUBTREE_COLOR, 0.42)
        elif clade in MAMMAL_CLADE_COLORS:
            color = rgba_from_hex(MAMMAL_CLADE_COLORS[clade], 0.40)
        else:
            color = "rgba(202,220,250,0.30)"
        internal_x.append(point[0])
        internal_y.append(point[1])
        internal_z.append(point[2])
        internal_size.append(INTERNAL_NODE_MARKER_MIN + ((INTERNAL_NODE_MARKER_MAX - INTERNAL_NODE_MARKER_MIN) * ((radius - 0.58) / 0.92)))
        internal_color.append(color)
        internal_hover.append(
            f"<b>{gene}</b><br>{attrs.get('label', 'internal split')}<br>"
            f"internal split supporting {support} descendant leaves<extra></extra>"
        )
    if internal_x:
        traces.append(
            go.Scatter3d(
                x=internal_x,
                y=internal_y,
                z=internal_z,
                mode="markers",
                marker=dict(
                    size=internal_size,
                    symbol="circle",
                    color=internal_color,
                    opacity=0.72,
                    line=dict(color="rgba(255,255,255,0.30)", width=0.45),
                ),
                customdata=internal_hover,
                hovertemplate="%{customdata}",
                showlegend=False,
                name=f"{gene} topology support nodes",
            )
        )

    for clade in MAMMAL_CLADE_ORDER:
        xs, ys, zs, labels, text = [], [], [], [], []
        for node, attrs in graph.nodes(data=True):
            if not attrs.get("is_mammal_ortholog") or attrs.get("mammal_clade") != clade:
                continue
            point = layout[node]
            xs.append(point[0])
            ys.append(point[1])
            zs.append(point[2])
            species_label = attrs.get("species_label", attrs.get("label", node))
            species = str(attrs.get("species", ""))
            labels.append(species_label if species in FOCAL_SPECIES_LABELS else "")
            text.append(
                f"<b>{gene}</b><br>{species_label}<br>"
                f"{attrs.get('orthology_type', '')}<br>"
                f"{float(attrs.get('percent_identity', 0.0)):.1f}% protein identity<br>"
                f"{float(attrs.get('alignment_coverage', 0.0)):.0%} aligned coverage<br>"
                f"{attrs.get('taxonomy_level', '')}<br>"
                f"{attrs.get('ensembl_gene_id', '')}<extra></extra>"
            )
        if not xs:
            continue
        traces.append(
            go.Scatter3d(
                x=xs,
                y=ys,
                z=zs,
                mode="markers",
                marker=dict(
                    size=NODE_MARKER_SIZE,
                    symbol="circle",
                    color=MAMMAL_CLADE_COLORS[clade],
                    opacity=0.94,
                    line=dict(color="rgba(255,255,255,0.72)", width=0.82),
                ),
                customdata=text,
                hovertemplate="%{customdata}",
                showlegend=(gene == GENE_ORDER[0]),
                name=f"Species ortholog · {clade}",
            )
        )

    human_x, human_y, human_z, human_hover = [], [], [], []
    for node, attrs in graph.nodes(data=True):
        if not attrs.get("is_leaf") or not attrs.get("is_human_variant_leaf"):
            continue
        point = layout[node]
        human_x.append(point[0])
        human_y.append(point[1])
        human_z.append(point[2])
        human_hover.append(
            f"<b>{gene}</b><br>Human variation subtree leaf<br>"
            "real gene-specific variant leaf attached at the inferred human sequence position<br>"
            f"{attrs.get('variant_id', '')}<br>"
            f"assigned allele frequency {float(attrs.get('assigned_af', 0.0)):.4g}<extra></extra>"
        )
    if human_x:
        traces.append(
            go.Scatter3d(
                x=human_x,
                y=human_y,
                z=human_z,
                mode="markers",
                marker=dict(
                    size=NODE_MARKER_SIZE,
                    symbol="circle",
                    color=HUMAN_SUBTREE_COLOR,
                    opacity=0.96,
                    line=dict(color="rgba(255,255,255,0.78)", width=0.9),
                ),
                customdata=human_hover,
                hovertemplate="%{customdata}",
                showlegend=(gene == GENE_ORDER[0]),
                name="Human variation subtree",
            )
        )

    mixed = [item for item in mixed_variant_candidates(graph) if item["node"] in layout]
    if mixed:
        traces.append(
            go.Scatter3d(
                x=[layout[item["node"]][0] for item in mixed],
                y=[layout[item["node"]][1] for item in mixed],
                z=[layout[item["node"]][2] + 8.0 for item in mixed],
                mode="markers",
                marker=dict(
                    size=NODE_MARKER_SIZE,
                    symbol="circle",
                    color=HUMAN_SUBTREE_COLOR,
                    opacity=0.98,
                    line=dict(color="rgba(255,255,255,0.85)", width=0.8),
                ),
                text=[
                    f"<b>{gene}</b><br>{item['variant_id']}<br>"
                    f"mixed-frequency evidence ratio {item['ratio']:.2f}<extra></extra>"
                    for item in mixed
                ],
                hovertemplate="%{text}",
                showlegend=(gene == GENE_ORDER[0]),
                name="Human mixed-frequency evidence",
            )
        )
    return traces


def centroid_for_nodes(
    graph: nx.DiGraph,
    layout: dict[str, np.ndarray],
    predicate,
) -> tuple[np.ndarray, int] | None:
    points = []
    count = 0
    for node, attrs in graph.nodes(data=True):
        if node not in layout or not predicate(attrs):
            continue
        points.append(layout[node])
        count += int(attrs.get("supporting_variant_count", attrs.get("taxon_count", 1)))
    if not points:
        return None
    return np.mean(np.array(points, dtype=float), axis=0), count


def flow_curve(start: np.ndarray, end: np.ndarray, index: int, total: int) -> np.ndarray:
    midpoint = (start[:2] + end[:2]) / 2.0
    norm = float(np.linalg.norm(midpoint))
    if norm < 1e-9:
        angle = math.tau * (index / max(total, 1))
        unit = np.array([math.cos(angle), math.sin(angle)], dtype=float)
    else:
        unit = midpoint / norm
    offset_xy = unit * FLOW_OUTWARD_OFFSET
    delta_z = float(end[2] - start[2])
    control_a = start + np.array([offset_xy[0], offset_xy[1], delta_z * 0.25], dtype=float)
    control_b = end + np.array([offset_xy[0], offset_xy[1], -delta_z * 0.25], dtype=float)
    t_values = np.linspace(0.0, 1.0, FLOW_CURVE_STEPS)
    points = []
    for t in t_values:
        inv = 1.0 - t
        point = (
            (inv**3) * start
            + (3.0 * inv * inv * t) * control_a
            + (3.0 * inv * t * t) * control_b
            + (t**3) * end
        )
        points.append(point)
    return np.array(points, dtype=float)


def append_flow_segment(
    coords: dict[str, list],
    curve: np.ndarray,
    hover_text: str,
) -> None:
    coords["x"].extend(curve[:, 0].tolist())
    coords["y"].extend(curve[:, 1].tolist())
    coords["z"].extend(curve[:, 2].tolist())
    coords["text"].extend([hover_text for _ in range(curve.shape[0])])
    coords["x"].append(None)
    coords["y"].append(None)
    coords["z"].append(None)
    coords["text"].append(None)


def taxon_flow_key(attrs: dict) -> str | None:
    if attrs.get("is_mammal_ortholog") and attrs.get("species"):
        return f"species::{attrs['species']}"
    return None


def taxon_flow_label(attrs: dict) -> str:
    if attrs.get("is_mammal_ortholog"):
        return str(attrs.get("species_label") or attrs.get("label") or attrs.get("species"))
    return str(attrs.get("label") or "Human variation")


def taxon_flow_color(attrs: dict) -> str:
    return MAMMAL_CLADE_COLORS.get(
        str(attrs.get("mammal_clade")),
        SPECIAL_GROUP_COLORS["Gene-history flow"],
    )


def taxon_flow_index(graph: nx.DiGraph) -> dict[str, str]:
    index = {}
    for node, attrs in graph.nodes(data=True):
        key = taxon_flow_key(attrs)
        if key:
            index[key] = node
    return index


def human_anchor_node(graph: nx.DiGraph) -> str | None:
    for node, attrs in graph.nodes(data=True):
        if attrs.get("is_human_sequence_anchor"):
            return node
    return None


def weighted_distance_between(graph: nx.DiGraph, source: str | None, target: str | None) -> float:
    if source is None or target is None:
        return 0.0
    tree = layout_graph(graph)
    try:
        return float(nx.shortest_path_length(tree, source, target, weight="layout_length"))
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return 0.0


def history_flow_records(
    graphs: dict[str, nx.DiGraph],
    layouts: dict[str, dict[str, np.ndarray]],
) -> list[dict]:
    records = []
    for source_gene, target_gene in zip(GENE_ORDER[:-1], GENE_ORDER[1:]):
        source_index = taxon_flow_index(graphs[source_gene])
        target_index = taxon_flow_index(graphs[target_gene])
        source_human = human_anchor_node(graphs[source_gene])
        target_human = human_anchor_node(graphs[target_gene])
        candidates = []
        for key in sorted(set(source_index) & set(target_index)):
            source_node = source_index[key]
            target_node = target_index[key]
            source_attrs = graphs[source_gene].nodes[source_node]
            target_attrs = graphs[target_gene].nodes[target_node]
            start = layouts[source_gene][source_node]
            end = layouts[target_gene][target_node]
            screen_shift = float(np.linalg.norm(end[:2] - start[:2]))
            source_human_distance = weighted_distance_between(graphs[source_gene], source_human, source_node)
            target_human_distance = weighted_distance_between(graphs[target_gene], target_human, target_node)
            human_distance_delta = abs(target_human_distance - source_human_distance)
            identity_delta = 0.0
            if source_attrs.get("is_mammal_ortholog"):
                identity_delta = abs(
                    float(target_attrs.get("percent_identity", 0.0))
                    - float(source_attrs.get("percent_identity", 0.0))
                )
            support_delta = 0
            if source_attrs.get("is_human_population_summary"):
                support_delta = abs(
                    int(target_attrs.get("supporting_variant_count", 0))
                    - int(source_attrs.get("supporting_variant_count", 0))
                )
            score = (
                screen_shift
                + (human_distance_delta * 260.0)
                + (identity_delta * HISTORY_FLOW_SCORE_IDENTITY_WEIGHT)
                + min(float(support_delta), 400.0) * 0.12
            )
            must_show = bool(
                source_attrs.get("species") in FOCAL_SPECIES_LABELS
            )
            candidates.append(
                {
                    "key": key,
                    "source_gene": source_gene,
                    "target_gene": target_gene,
                    "source_node": source_node,
                    "target_node": target_node,
                    "source_attrs": source_attrs,
                    "target_attrs": target_attrs,
                    "start": start,
                    "end": end,
                    "screen_shift": screen_shift,
                    "human_distance_delta": human_distance_delta,
                    "identity_delta": identity_delta,
                    "support_delta": support_delta,
                    "score": score,
                    "must_show": must_show,
                }
            )
        must = [item for item in candidates if item["must_show"]]
        optional = sorted(
            [item for item in candidates if not item["must_show"]],
            key=lambda item: (-item["score"], item["key"]),
        )
        selected = must + optional[: max(0, HISTORY_FLOW_TOP_PER_PAIR - len(must))]
        records.extend(sorted(selected, key=lambda item: (item["source_gene"], -item["score"], item["key"])))
    return records


def gene_history_flow_traces(
    graphs: dict[str, nx.DiGraph],
    layouts: dict[str, dict[str, np.ndarray]],
) -> list[go.Scatter3d]:
    traces = []
    records = history_flow_records(graphs, layouts)
    for index, record in enumerate(records):
        source_attrs = record["source_attrs"]
        target_attrs = record["target_attrs"]
        color = taxon_flow_color(source_attrs)
        curve = flow_curve(record["start"], record["end"], index, max(1, len(records)))
        hover = (
            f"<b>Same species across gene trees</b><br>{taxon_flow_label(source_attrs)}<br>"
            f"{record['source_gene']} -> {record['target_gene']}<br>"
            f"human-branch distance change: {record['human_distance_delta']:.4f}<br>"
            f"protein identity change: {record['identity_delta']:.2f} percentage points<br>"
            f"screen movement: {record['screen_shift']:.1f}<br>"
            "The curve links the same taxon between independently inferred gene trees.<extra></extra>"
        )
        width = min(1.45, 0.48 + (record["score"] / 360.0))
        alpha = 0.18 if not record["must_show"] else 0.30
        traces.append(
            go.Scatter3d(
                x=curve[:, 0],
                y=curve[:, 1],
                z=curve[:, 2],
                mode="lines",
                line=dict(color=rgba_from_hex(color, alpha), width=width),
                text=[hover for _ in range(curve.shape[0])],
                hovertemplate="%{text}",
                showlegend=False,
                name=f"Gene-history flow · {taxon_flow_label(source_attrs)}",
            )
        )
    return traces


def inter_gene_flow_traces(
    graphs: dict[str, nx.DiGraph],
    layouts: dict[str, dict[str, np.ndarray]],
) -> list[go.Scatter3d]:
    traces = []
    adjacent_pairs = list(zip(GENE_ORDER[:-1], GENE_ORDER[1:]))

    represented_clades = [
        clade
        for clade in MAMMAL_CLADE_ORDER
        if any(mammal_clade_counts(graphs[gene]).get(clade, 0) for gene in GENE_ORDER)
    ]
    for clade_index, clade in enumerate(represented_clades):
        coords = {"x": [], "y": [], "z": [], "text": []}
        counts = []
        for source_gene, target_gene in adjacent_pairs:
            source = centroid_for_nodes(
                graphs[source_gene],
                layouts[source_gene],
                lambda attrs, group=clade: (
                    attrs.get("is_mammal_ortholog") and attrs.get("mammal_clade") == group
                ),
            )
            target = centroid_for_nodes(
                graphs[target_gene],
                layouts[target_gene],
                lambda attrs, group=clade: (
                    attrs.get("is_mammal_ortholog") and attrs.get("mammal_clade") == group
                ),
            )
            if not source or not target:
                continue
            source_centroid, source_count = source
            target_centroid, target_count = target
            counts.extend([source_count, target_count])
            curve = flow_curve(source_centroid, target_centroid, clade_index, len(represented_clades))
            shift = float(np.linalg.norm(target_centroid[:2] - source_centroid[:2]))
            hover = (
                f"<b>Species clade flow</b><br>{clade}<br>"
                f"{source_gene} ortholog leaves: {source_count}<br>"
                f"{target_gene} ortholog leaves: {target_count}<br>"
                f"screen centroid shift: {shift:.1f}<br>"
                "Aggregated ortholog-context guide; not physical gene transfer.<extra></extra>"
            )
            append_flow_segment(coords, curve, hover)
        if not coords["x"]:
            continue
        mean_count = float(np.mean(counts)) if counts else 0.0
        line_width = min(1.25, 0.55 + (mean_count / 50.0))
        traces.append(
            go.Scatter3d(
                x=coords["x"],
                y=coords["y"],
                z=coords["z"],
                mode="lines",
                line=dict(color=rgba_from_hex(MAMMAL_CLADE_COLORS[clade], 0.18), width=line_width),
                text=coords["text"],
                hovertemplate="%{text}",
                showlegend=False,
                name=f"Species flow · {clade}",
            )
        )
    return traces


def population_leaf_counts(graph: nx.DiGraph) -> dict[str, int]:
    return {
        population: int(
            sum(
                int(attrs.get("supporting_variant_count", 1))
                for _, attrs in graph.nodes(data=True)
                if attrs.get("is_leaf")
                and not attrs.get("is_mammal_ortholog")
                and attrs.get("assigned_population") == population
            )
        )
        for population in POP_ORDER
    }


def human_variant_leaf_count(graph: nx.DiGraph) -> int:
    return sum(
        1
        for _, attrs in graph.nodes(data=True)
        if attrs.get("is_leaf") and not attrs.get("is_mammal_ortholog")
    )


def supporting_human_variant_count(graph: nx.DiGraph) -> int:
    if graph.graph.get("supporting_human_variant_leaves"):
        return int(graph.graph["supporting_human_variant_leaves"])
    return human_variant_leaf_count(graph)


def mixed_candidate_count(graph: nx.DiGraph) -> int:
    if graph.graph.get("supporting_mixed_candidates") is not None:
        return int(graph.graph["supporting_mixed_candidates"])
    return len(mixed_variant_candidates(graph))


def mammal_ortholog_count(graph: nx.DiGraph) -> int:
    return sum(1 for _, attrs in graph.nodes(data=True) if attrs.get("is_mammal_ortholog"))


def mammal_clade_counts(graph: nx.DiGraph) -> dict[str, int]:
    return {
        clade: sum(
            1
            for _, attrs in graph.nodes(data=True)
            if attrs.get("is_mammal_ortholog") and attrs.get("mammal_clade") == clade
        )
        for clade in MAMMAL_CLADE_ORDER
    }


def taxon_color_audit(graphs: dict[str, nx.DiGraph]) -> dict:
    represented_clades = [
        clade
        for clade in MAMMAL_CLADE_ORDER
        if any(mammal_clade_counts(graphs[gene]).get(clade, 0) for gene in GENE_ORDER)
    ]
    represented_colors = {
        clade: MAMMAL_CLADE_COLORS.get(clade, MAMMAL_CLADE_COLORS["Species ortholog context"])
        for clade in represented_clades
    }
    reverse: dict[str, list[str]] = defaultdict(list)
    for clade, color in represented_colors.items():
        reverse[color.lower()].append(clade)
    duplicate_colors = {
        color: clades
        for color, clades in reverse.items()
        if len(clades) > 1
    }
    return {
        "ortholog_color_unit": "taxonomic clade",
        "human_subtree_color": HUMAN_SUBTREE_COLOR,
        "represented_clade_count": len(represented_clades),
        "represented_clade_colors": represented_colors,
        "duplicate_represented_clade_colors": duplicate_colors,
        "human_color_reused_by_ortholog_clade": HUMAN_SUBTREE_COLOR.lower()
        in {color.lower() for color in represented_colors.values()},
    }


def leaf_rendering_audit(graphs: dict[str, nx.DiGraph]) -> dict:
    return {
        "all_visible_node_symbols": "circle",
        "all_visible_node_marker_size": NODE_MARKER_SIZE,
        "ortholog_size_policy": "constant radius for every rendered species/outgroup node",
        "human_subtree_size_policy": "constant radius matching species/outgroup nodes",
        "mixed_evidence_overlay_size_policy": "constant radius matching all other rendered node markers",
        "hyperbolic_terminal_z_stagger": HYPERBOLIC_LEAF_Z_STAGGER,
    }


def story_metric_payload(graphs: dict[str, nx.DiGraph]) -> dict:
    genes = {}
    totals = {
        "selected_human_taxa": 0,
        "displayed_aligned_taxa": 0,
        "displayed_human_subtree_leaves": 0,
        "source_variants": 0,
        "mixed_candidates": 0,
        "mammal_orthologs": 0,
        "ape_orthologs": 0,
    }
    represented_clades = set()

    for gene in GENE_ORDER:
        graph = graphs[gene]
        selected = supporting_human_variant_count(graph)
        human_taxa = human_variant_leaf_count(graph)
        mammals = mammal_ortholog_count(graph)
        clade_counts = mammal_clade_counts(graph)
        apes = clade_counts.get("Great apes", 0) + clade_counts.get("Lesser apes", 0)
        source = graph.graph.get("total_variants_in_source") or "?"
        mixed_count = mixed_candidate_count(graph)
        displayed_taxa = mammals + human_taxa

        for clade, count in clade_counts.items():
            if count:
                represented_clades.add(clade)

        totals["selected_human_taxa"] += selected
        totals["displayed_aligned_taxa"] += displayed_taxa
        totals["displayed_human_subtree_leaves"] += human_taxa
        totals["mammal_orthologs"] += mammals
        totals["ape_orthologs"] += apes
        totals["mixed_candidates"] += mixed_count
        if isinstance(source, int):
            totals["source_variants"] += source
            source_label = f"{source:,}"
        else:
            source_label = str(source)

        genes[gene] = {
            "selected_human_taxa": selected,
            "displayed_aligned_taxa": displayed_taxa,
            "displayed_human_subtree_leaves": human_taxa,
            "source_variants": source,
            "source_label": source_label,
            "mammal_orthologs": mammals,
            "ape_orthologs": apes,
            "mammal_clade_counts": clade_counts,
            "mixed_candidates": mixed_count,
        }

    return {
        "genes": genes,
        "totals": totals,
        "represented_mammal_clades": [
            clade for clade in MAMMAL_CLADE_ORDER if clade in represented_clades
        ],
    }


def inter_gene_flow_summary(graphs: dict[str, nx.DiGraph]) -> dict:
    adjacent_pairs = list(zip(GENE_ORDER[:-1], GENE_ORDER[1:]))
    clade_counts = {gene: mammal_clade_counts(graphs[gene]) for gene in GENE_ORDER}
    mammal_segments = 0
    represented_clades = [
        clade
        for clade in MAMMAL_CLADE_ORDER
        if any(clade_counts[gene].get(clade, 0) for gene in GENE_ORDER)
    ]
    for source_gene, target_gene in adjacent_pairs:
        mammal_segments += sum(
            1
            for clade in represented_clades
            if clade_counts[source_gene].get(clade, 0)
            and clade_counts[target_gene].get(clade, 0)
        )
    return {
        "mode": "enabled",
        "style": "curved aggregate centroid flow guides",
        "adjacent_gene_pairs": [f"{a}->{b}" for a, b in adjacent_pairs],
        "mammal_clade_groups": len(represented_clades),
        "mammal_clade_flow_segments": mammal_segments,
        "total_flow_segments": mammal_segments,
        "aggregation": (
            "species flows connect clade centroids for the shared ortholog species panel"
        ),
        "interpretation_caveat": (
            "flows are correspondence and screen-shift guides across gene trees, not physical gene transfer "
            "or proof of recombination"
        ),
    }


def gene_history_flow_summary(
    graphs: dict[str, nx.DiGraph],
    layouts: dict[str, dict[str, np.ndarray]],
) -> dict:
    records = history_flow_records(graphs, layouts)
    species_records = [record for record in records if record["source_attrs"].get("is_mammal_ortholog")]
    top_records = sorted(records, key=lambda item: (-item["score"], item["key"]))[:10]
    return {
        "mode": "enabled",
        "style": "curved same-taxon gene-history correspondence links",
        "adjacent_gene_pairs": [f"{a}->{b}" for a, b in zip(GENE_ORDER[:-1], GENE_ORDER[1:])],
        "visible_flow_segments": len(records),
        "species_or_outgroup_segments": len(species_records),
        "selection_policy": (
            "always show focal species/outgroups; fill remaining slots with the largest "
            "same-taxon shifts by branch-distance-to-human, protein identity change, and screen movement"
        ),
        "real_data_basis": (
            "species links use independently inferred gene-tree positions plus protein identity and branch-distance "
            "changes from Ensembl homology alignments"
        ),
        "interpretation_caveat": (
            "curves show how the same taxon lands in different gene trees; they are not physical gene-transfer paths"
        ),
        "top_shift_examples": [
            {
                "taxon": taxon_flow_label(record["source_attrs"]),
                "pair": f"{record['source_gene']}->{record['target_gene']}",
                "score": round(float(record["score"]), 3),
                "human_branch_distance_delta": round(float(record["human_distance_delta"]), 5),
                "protein_identity_delta": round(float(record["identity_delta"]), 3),
                "support_delta": int(record["support_delta"]),
            }
            for record in top_records
        ],
    }


def conclusion_sentence(metrics: dict) -> str:
    return (
        "Each housekeeping gene has its own inferred ortholog tree on the same species/outgroups; the neutral "
        "human subtree keeps the human sequence context expanded with one consistent visual style."
    )


def story_annotations(graphs: dict[str, nx.DiGraph], geometry: str = "box") -> list[dict]:
    return []


def scene_ranges(layouts: dict[str, dict[str, np.ndarray]], geometry: str = "box") -> tuple[list[float], list[float], list[float]]:
    points = []
    for layout in layouts.values():
        points.extend(layout.values())
    if geometry != "hyperbolic":
        for z in GENE_Z.values():
            z_low = z - (LOCAL_Z_SPAN / 2.0)
            z_high = z + (LOCAL_Z_SPAN / 2.0)
            points.extend(
                [
                    np.array([-BOX_HALF_SIZE, -BOX_HALF_SIZE, z_low], dtype=float),
                    np.array([BOX_HALF_SIZE, BOX_HALF_SIZE, z_high], dtype=float),
                ]
            )
    matrix = np.array(points, dtype=float)

    def padded(values: np.ndarray, minimum: float) -> list[float]:
        low = float(np.min(values))
        high = float(np.max(values))
        pad = max(minimum, (high - low) * 0.10)
        return [low - pad, high + pad]

    return padded(matrix[:, 0], 48.0), padded(matrix[:, 1], 48.0), padded(matrix[:, 2], 80.0)


def scene_aspect_ratio(geometry: str = "box") -> dict[str, float]:
    if geometry == "hyperbolic":
        return dict(x=2.05, y=1.58, z=2.18)
    return dict(x=1.18, y=1.18, z=2.05)


def build_figure(graphs: dict[str, nx.DiGraph], geometry: str = "box") -> tuple[go.Figure, dict[str, dict[str, np.ndarray]]]:
    layouts = {gene: compute_unrooted_3d_layout(graphs[gene], gene, geometry) for gene in GENE_ORDER}
    figure = go.Figure()

    for gene in GENE_ORDER:
        if geometry == "hyperbolic":
            pass
        else:
            figure.add_trace(plane_trace(gene))
        figure.add_trace(gene_label_trace(gene))

    for gene in GENE_ORDER:
        for trace in edge_traces(graphs[gene], layouts[gene], gene):
            figure.add_trace(trace)

    for trace in gene_history_flow_traces(graphs, layouts):
        figure.add_trace(trace)

    for gene in GENE_ORDER:
        for trace in node_traces(graphs[gene], layouts[gene], gene):
            figure.add_trace(trace)

    x_range, y_range, z_range = scene_ranges(layouts, geometry)
    aligned_species_count = mammal_ortholog_count(graphs[GENE_ORDER[0]])
    aligned_human_count = human_variant_leaf_count(graphs[GENE_ORDER[0]])
    if geometry == "hyperbolic":
        title_text = (
            "<b>Inferred stacked 3D housekeeping-gene trees</b><br>"
            f"<sup>{aligned_species_count} shared ortholog species + "
            f"{aligned_human_count} neutral human-subtree leaves in every gene layer</sup>"
        )
    else:
        title_text = (
            "<b>Inferred stacked 3D housekeeping-gene trees</b><br>"
            f"<sup>{aligned_species_count} shared ortholog species + "
            f"{aligned_human_count} neutral human-subtree leaves in every gene layer</sup>"
        )
    figure.update_layout(
        title=dict(
            text=title_text,
            x=0.5,
            xanchor="center",
            font=dict(size=24, color="#F4F8FF", family="Avenir Next, IBM Plex Sans, sans-serif"),
        ),
        template="plotly_dark",
        paper_bgcolor="#06101E",
        plot_bgcolor="#06101E",
        font=dict(family="Avenir Next, IBM Plex Sans, sans-serif"),
        height=940,
        margin=dict(l=0, r=0, t=104, b=0),
        legend=dict(
            title="Taxon clade",
            bgcolor="rgba(6,16,30,0.78)",
            bordercolor="rgba(177,203,238,0.22)",
            borderwidth=1,
            x=0.01,
            y=0.98,
        ),
        annotations=story_annotations(graphs, geometry),
        scene=dict(
            bgcolor="#06101E",
            xaxis=dict(title="", range=x_range, visible=False),
            yaxis=dict(title="", range=y_range, visible=False),
            zaxis=dict(title="", range=z_range, visible=False),
            camera=dict(eye=dict(x=1.92, y=1.62, z=1.18), up=dict(x=0, y=0, z=1)),
            aspectmode="manual",
            aspectratio=scene_aspect_ratio(geometry),
        ),
    )
    return figure, layouts


def orientation(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    return float((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def proper_segment_intersection(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> bool:
    eps = 1e-9
    o1 = orientation(a, b, c)
    o2 = orientation(a, b, d)
    o3 = orientation(c, d, a)
    o4 = orientation(c, d, b)
    return (o1 * o2 < -eps) and (o3 * o4 < -eps)


def count_tree_edge_crossings(graph: nx.DiGraph, layout: dict[str, np.ndarray]) -> int:
    edges = list(graph.edges())
    crossings = 0
    for idx, (a_node, b_node) in enumerate(edges):
        a = layout[a_node][:2]
        b = layout[b_node][:2]
        for c_node, d_node in edges[idx + 1 :]:
            if len({a_node, b_node, c_node, d_node}) < 4:
                continue
            c = layout[c_node][:2]
            d = layout[d_node][:2]
            if proper_segment_intersection(a, b, c, d):
                crossings += 1
    return crossings


def tree_structure_metrics(graph: nx.DiGraph) -> dict:
    undirected = layout_graph(graph)
    internal_nodes = [
        node
        for node, attrs in graph.nodes(data=True)
        if not attrs.get("is_leaf")
    ]
    return {
        "is_undirected_tree": bool(nx.is_tree(undirected)),
        "is_directed_arborescence": bool(nx.is_arborescence(graph)),
        "self_loop_edges": len(list(nx.selfloop_edges(graph))),
        "nodes_minus_one_equals_edges": bool(graph.number_of_edges() == graph.number_of_nodes() - 1),
        "max_internal_out_degree": max((graph.out_degree(node) for node in internal_nodes), default=0),
        "max_internal_undirected_degree": max((undirected.degree(node) for node in internal_nodes), default=0),
        "unrooted_binary_check": bool(
            nx.is_tree(undirected)
            and not list(nx.selfloop_edges(graph))
            and max((undirected.degree(node) for node in internal_nodes), default=0) <= 3
        ),
    }


def rendered_leaf_spacing_metrics(
    graph: nx.DiGraph,
    layout: dict[str, np.ndarray],
    x_range: list[float],
    y_range: list[float],
    z_range: list[float],
    aspectratio: dict[str, float],
) -> dict:
    points_xy = []
    points_xyz = []
    displayed_human_nodes = displayed_human_leaf_nodes(graph)
    for node, attrs in graph.nodes(data=True):
        if not attrs.get("is_leaf"):
            continue
        if attrs.get("is_mammal_ortholog") or node in displayed_human_nodes:
            point = layout[node]
            x_span = max(float(x_range[1] - x_range[0]), 1e-9)
            y_span = max(float(y_range[1] - y_range[0]), 1e-9)
            z_span = max(float(z_range[1] - z_range[0]), 1e-9)
            xy = np.array(
                [
                    ((point[0] - x_range[0]) / x_span) * aspectratio["x"],
                    ((point[1] - y_range[0]) / y_span) * aspectratio["y"],
                ],
                dtype=float,
            )
            xyz = np.array(
                [
                    xy[0],
                    xy[1],
                    ((point[2] - z_range[0]) / z_span) * aspectratio["z"],
                ],
                dtype=float,
            )
            points_xy.append(xy)
            points_xyz.append(xyz)
    if len(points_xy) < 2:
        return {
            "rendered_leaf_count": len(points_xy),
            "nearest_neighbor_xy_min": None,
            "nearest_neighbor_xy_p05": None,
            "nearest_neighbor_xy_median": None,
            "nearest_neighbor_3d_min": None,
            "nearest_neighbor_3d_p05": None,
            "nearest_neighbor_3d_median": None,
        }

    def nearest_summary(points: list[np.ndarray]) -> dict[str, float]:
        matrix = np.array(points, dtype=float)
        nearest = []
        for index, point in enumerate(matrix):
            deltas = matrix - point
            distances = np.sqrt(np.sum(deltas * deltas, axis=1))
            distances[index] = np.inf
            nearest.append(float(np.min(distances)))
        nearest_array = np.array(nearest, dtype=float)
        return {
            "min": round(float(np.min(nearest_array)), 5),
            "p05": round(float(np.quantile(nearest_array, 0.05)), 5),
            "median": round(float(np.median(nearest_array)), 5),
        }

    xy_summary = nearest_summary(points_xy)
    xyz_summary = nearest_summary(points_xyz)
    return {
        "rendered_leaf_count": len(points_xy),
        "nearest_neighbor_xy_min": xy_summary["min"],
        "nearest_neighbor_xy_p05": xy_summary["p05"],
        "nearest_neighbor_xy_median": xy_summary["median"],
        "nearest_neighbor_3d_min": xyz_summary["min"],
        "nearest_neighbor_3d_p05": xyz_summary["p05"],
        "nearest_neighbor_3d_median": xyz_summary["median"],
        "normalized_by_scene_aspect": True,
    }


def layer_separation_metrics(layouts: dict[str, dict[str, np.ndarray]]) -> dict:
    ranges = {}
    for gene in GENE_ORDER:
        matrix = np.array(list(layouts[gene].values()), dtype=float)
        ranges[gene] = [float(np.min(matrix[:, 2])), float(np.max(matrix[:, 2]))]
    gaps = {}
    for lower_gene, upper_gene in zip(GENE_ORDER[:-1], GENE_ORDER[1:]):
        gaps[f"{lower_gene}->{upper_gene}"] = float(ranges[upper_gene][0] - ranges[lower_gene][1])
    return {
        "configured_layer_z_separation": LAYER_Z_SEPARATION,
        "gene_z_centers": GENE_Z,
        "layer_z_ranges": ranges,
        "clear_z_gaps_between_adjacent_layers": gaps,
    }


def export_layout_audit(
    graphs: dict[str, nx.DiGraph],
    layouts: dict[str, dict[str, np.ndarray]],
    output_path: str = LAYOUT_AUDIT_PATH,
    geometry: str = "box",
) -> None:
    x_range, y_range, z_range = scene_ranges(layouts, geometry)
    aspectratio = scene_aspect_ratio(geometry)
    if geometry == "hyperbolic":
        layout_algorithm = (
            "equal-angle unrooted phylogenetic tree layout with disjoint leaf-count angular wedges; "
            "topology-preserving uniform fit inside a bounded Poincare-style disk; "
            "local z encodes normalized cumulative branch distance from the selected unrooted center"
        )
        geometry_payload = {
            "mode": "hyperbolic",
            "hyperbolic_radius": HYPERBOLIC_RADIUS,
            "hyperbolic_fill_fraction": HYPERBOLIC_FILL_FRACTION,
            "local_z_span": LOCAL_Z_SPAN,
            "visible_substrate": "none; Poincare disk surface and circular grid traces are hidden",
        }
    else:
        layout_algorithm = (
            "equal-angle unrooted phylogenetic tree layout with disjoint leaf-count angular wedges; "
            "affine cuboid fitting preserves the planar no-crossing embedding; "
            "local z encodes normalized cumulative branch distance from the selected unrooted center"
        )
        geometry_payload = {
            "mode": "box",
            "box_half_size": BOX_HALF_SIZE,
            "local_z_span": LOCAL_Z_SPAN,
            "box_fill_fraction": BOX_FILL_FRACTION,
        }
    metrics = story_metric_payload(graphs)
    flow_summary = gene_history_flow_summary(graphs, layouts)
    tree_methods_by_gene = {
        gene: graphs[gene].graph.get("tree_inference_method", "")
        for gene in GENE_ORDER
    }
    tree_methods = sorted({method for method in tree_methods_by_gene.values() if method})
    audit = {
        "layout_algorithm": layout_algorithm,
        "graph_interpretation": (
            "undirected unrooted NetworkX gene trees inferred from aligned ortholog protein sequences; "
            "the inferred human sequence leaf is expanded into the same neutral human subtree; "
            "same-species history-flow links connect adjacent gene layers as correspondence guides only"
        ),
        "tree_inference": {
            "method": tree_methods[0] if len(tree_methods) == 1 else "; ".join(tree_methods),
            "methods_by_gene": tree_methods_by_gene,
            "distance": "p-distance over shared aligned amino-acid positions",
            "sequence_source": "NCBI Datasets ortholog protein packages aligned with MAFFT when outputs_3d/ncbi_*_all_homologies.json exists; otherwise Ensembl REST homology align_seq payloads",
            "minimum_target_coverage": HOMOLOGY_MIN_COVERAGE,
            "human_expansion": "the inferred human protein leaf is expanded into neutral human-subtree leaves",
        },
        "geometry": geometry_payload,
        "inter_tree_connections": flow_summary,
        "color_encoding": taxon_color_audit(graphs),
        "leaf_rendering": leaf_rendering_audit(graphs),
        "internal_clade_markers": "hidden",
        "story": {
            "headline": "Same taxa, different gene-version histories",
            "conclusion": conclusion_sentence(metrics),
            "reading_path": [
                "each gene layer displays the same aligned taxa: shared ortholog species plus neutral human-subtree leaves",
                "ortholog species topology is independently inferred per gene from aligned protein homology sequences",
                "the inferred human protein leaf is expanded into an unlabeled local human subtree",
                "supporting human-variant counts are retained in hover/audit metadata using a neutral display",
                "curved flow guides connect the same species or outgroup between adjacent gene trees",
            ],
            "taxon_meaning": (
                "displayed taxa are aligned across all genes: ortholog leaves are the same Ensembl species; "
                "human leaves are neutral human-subtree leaves attached at the inferred human sequence position"
            ),
            "human_display_policy": (
                "human variants are summarized into neutral subtree leaves; ortholog layers without "
                "cached gene-specific human variant trees reuse the ACTB subtree as an explicitly marked proxy"
            ),
            "mammal_taxon_filter": "only species with available ortholog records for all plotted genes are rendered",
            "human_taxon_filter": "neutral human-subtree leaves are expanded from the inferred human sequence leaf for every displayed layer",
            "flow_meaning": (
                "curved lines are same-species/outgroup correspondence guides connected between adjacent independently inferred gene trees"
            ),
            "color_meaning": "ortholog leaves use unique clade colors; the human subtree uses one neutral color",
            "leaf_marker_meaning": "all rendered leaves use the same circular marker shape",
            "interpretation_caveat": (
                "same-taxon flow guides show gene-history shifts, not physical gene-transfer paths"
            ),
            "represented_mammal_clades": metrics["represented_mammal_clades"],
            "totals": metrics["totals"],
            "mammal_context_source": "Ensembl REST homology/symbol endpoint",
        },
        "scene_ranges": {"x": x_range, "y": y_range, "z": z_range},
        "scene_aspectratio": aspectratio,
        "layer_separation": layer_separation_metrics(layouts),
        "genes": {},
    }
    for gene in GENE_ORDER:
        graph = graphs[gene]
        layout = layouts[gene]
        matrix = np.array(list(layout.values()), dtype=float)
        xy_span = matrix[:, :2].max(axis=0) - matrix[:, :2].min(axis=0)
        xy_extent = 2.0 * (HYPERBOLIC_RADIUS if geometry == "hyperbolic" else BOX_HALF_SIZE)
        z_range_gene = [float(np.min(matrix[:, 2])), float(np.max(matrix[:, 2]))]
        audit["genes"][gene] = {
            "nodes": graph.number_of_nodes(),
            "edges": graph.number_of_edges(),
            "leaves": sum(1 for _, attrs in graph.nodes(data=True) if attrs.get("is_leaf")),
            "human_variant_leaves": human_variant_leaf_count(graph),
            "displayed_human_subtree_leaves": human_variant_leaf_count(graph),
            "supporting_human_variant_leaves": supporting_human_variant_count(graph),
            "visible_human_markers": len(displayed_human_leaf_nodes(graph)),
            "mammal_ortholog_leaves": mammal_ortholog_count(graph),
            "mammal_clade_counts": mammal_clade_counts(graph),
            "source_variants": graph.graph.get("total_variants_in_source", ""),
            "source_selection": graph.graph.get("source_selection", ""),
            "tree_inference_method": graph.graph.get("tree_inference_method", ""),
            "species_backbone": graph.graph.get("species_backbone", ""),
            "homology_min_coverage": graph.graph.get("homology_min_coverage", ""),
            "mammal_context_source": graph.graph.get("mammal_context_source", ""),
            "mammal_context_source_url": graph.graph.get("mammal_context_source_url", ""),
            "mixed_candidates": len(mixed_variant_candidates(graph)),
            "coordinates_finite": bool(np.isfinite(matrix).all()),
            "coordinates_inside_scene": bool(
                (matrix[:, 0] >= x_range[0]).all()
                and (matrix[:, 0] <= x_range[1]).all()
                and (matrix[:, 1] >= y_range[0]).all()
                and (matrix[:, 1] <= y_range[1]).all()
                and (matrix[:, 2] >= z_range[0]).all()
                and (matrix[:, 2] <= z_range[1]).all()
            ),
            "tree_edge_crossings_xy": count_tree_edge_crossings(graph, layout),
            "tree_structure": tree_structure_metrics(graph),
            "rendered_leaf_spacing": rendered_leaf_spacing_metrics(graph, layout, x_range, y_range, z_range, aspectratio),
            "layout_center": choose_layout_center(layout_graph(graph)),
            "xy_span": [float(xy_span[0]), float(xy_span[1])],
            "xy_occupancy": [
                float(xy_span[0] / xy_extent),
                float(xy_span[1] / xy_extent),
            ],
            "local_z_range": z_range_gene,
            "local_z_span_used": float(z_range_gene[1] - z_range_gene[0]),
        }
    with open(output_path, "w") as fh:
        json.dump(audit, fh, indent=2)


def main() -> None:
    graphs = load_graphs()
    mammal_context = load_mammal_context()
    human_subtree_source_gene = "ACTB"
    gene_source_graphs = {}
    for gene in GENE_ORDER:
        if gene in graphs:
            gene_source_graphs[gene] = graphs[gene]
            continue
        proxy = graphs[human_subtree_source_gene].copy()
        proxy.graph.update(graphs[human_subtree_source_gene].graph)
        proxy.graph["human_subtree_proxy_source_gene"] = human_subtree_source_gene
        proxy.graph["human_subtree_proxy_note"] = (
            f"{gene} ortholog layer uses the existing {human_subtree_source_gene} neutral human-subtree "
            "variant display because no gene-specific human variant tree is cached for this gene"
        )
        gene_source_graphs[gene] = proxy
    selected = {
        gene: compose_inferred_taxon_aligned_tree(gene_source_graphs[gene], gene, mammal_context)
        for gene in GENE_ORDER
    }
    figure, layouts = build_figure(selected, geometry="box")
    export_layout_audit(selected, layouts, LAYOUT_AUDIT_PATH, geometry="box")

    for output_path in (OUTPUT_MAIN, OUTPUT_UNROOTED, OUTPUT_STACKED):
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        figure.write_html(output_path, include_plotlyjs=True, full_html=True)
        print(f"Wrote {output_path}")
    print(f"Wrote {LAYOUT_AUDIT_PATH}")

    hyperbolic_figure, hyperbolic_layouts = build_figure(selected, geometry="hyperbolic")
    export_layout_audit(selected, hyperbolic_layouts, HYPERBOLIC_AUDIT_PATH, geometry="hyperbolic")
    hyperbolic_figure.write_html(OUTPUT_HYPERBOLIC, include_plotlyjs=True, full_html=True)
    print(f"Wrote {OUTPUT_HYPERBOLIC}")
    print(f"Wrote {HYPERBOLIC_AUDIT_PATH}")


if __name__ == "__main__":
    main()
