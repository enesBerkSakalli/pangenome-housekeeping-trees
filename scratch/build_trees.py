"""
Process gnomAD population allele frequency data and build UPGMA Newick trees per gene.
Outputs data.js for the Pangenome Tree Visualizer.
"""
import json
import math
import os
import glob

SCRATCH_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(os.path.dirname(SCRATCH_DIR))

POPULATIONS = {
    "afr": "African",
    "amr": "American",
    "asj": "Ashkenazi Jewish",
    "eas": "East Asian",
    "fin": "Finnish",
    "mid": "Middle Eastern",
    "nfe": "Non-Finnish European",
    "sas": "South Asian",
    "remaining": "Remaining",
}

POP_ORDER = ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas"]  # exclude 'remaining' for clarity

GENES = {
    "LCT":   {"file": "LCT_variants.json",   "description": "Lactase (dairy tolerance adaptation)"},
    "EPAS1": {"file": "EPAS1_variants.json",  "description": "HIF2α (high-altitude adaptation)"},
    "HLA-A": {"file": "HLA-A_variants.json",  "description": "HLA Class I (immune recognition)"},
    "GAPDH": {"file": "GAPDH_variants.json",  "description": "GAPDH (housekeeping control)"},
    "MAPT":  {"file": "MAPT_variants.json",   "description": "Microtubule-associated protein tau"},
}


def load_variants(filepath):
    with open(filepath) as f:
        data = json.load(f)
    # Try gene key first, then region
    try:
        variants = data["data"]["gene"]["variants"]
    except (KeyError, TypeError):
        try:
            variants = data["data"]["region"]["variants"]
        except (KeyError, TypeError):
            variants = []
    return variants


def compute_pop_af_vectors(variants):
    """For each population, compute a list of allele frequencies across all variants."""
    pop_ac = {p: 0 for p in POP_ORDER}
    pop_an = {p: 0 for p in POP_ORDER}

    for variant in variants:
        for source in ["exome", "genome"]:
            src = variant.get(source)
            if not src:
                continue
            pops = src.get("populations", [])
            for pop in pops:
                pid = pop.get("id")
                if pid in POP_ORDER:
                    pop_ac[pid] += pop.get("ac", 0)
                    pop_an[pid] += pop.get("an", 0)

    # Compute AF per population
    pop_af = {}
    for p in POP_ORDER:
        if pop_an[p] > 0:
            pop_af[p] = pop_ac[p] / pop_an[p]
        else:
            pop_af[p] = 0.0
    return pop_af


def fst_distance(af1, af2):
    """
    Simple pairwise distance based on allele frequency difference.
    Uses a vector of per-variant AFs to compute average absolute difference.
    """
    return abs(af1 - af2)


def compute_distance_matrix(pop_afs):
    """Compute pairwise distance matrix between populations."""
    n = len(POP_ORDER)
    matrix = [[0.0] * n for _ in range(n)]
    for i, p1 in enumerate(POP_ORDER):
        for j, p2 in enumerate(POP_ORDER):
            if i == j:
                matrix[i][j] = 0.0
            elif i < j:
                d = fst_distance(pop_afs[p1], pop_afs[p2])
                matrix[i][j] = d
                matrix[j][i] = d
    return matrix


def upgma(matrix, labels):
    """
    UPGMA clustering to build a Newick tree.
    Returns a Newick string.
    """
    import copy
    n = len(labels)
    dists = copy.deepcopy(matrix)
    # clusters: list of (newick_str, size)
    clusters = [(labels[i], 1) for i in range(n)]
    heights = [0.0] * n

    while len(clusters) > 1:
        # Find minimum distance
        best_d = float("inf")
        best_i, best_j = 0, 1
        for i in range(len(clusters)):
            for j in range(i+1, len(clusters)):
                if dists[i][j] < best_d:
                    best_d = dists[i][j]
                    best_i, best_j = i, j

        # Heights of new node
        new_height = best_d / 2.0
        bl_i = new_height - heights[best_i]
        bl_j = new_height - heights[best_j]

        # Format branch lengths to 4 decimal places
        ni, si = clusters[best_i]
        nj, sj = clusters[best_j]
        new_node = f"({ni}:{bl_i:.4f},{nj}:{bl_j:.4f})"
        new_size = si + sj

        # Compute new distances
        new_dists_row = []
        for k in range(len(clusters)):
            if k == best_i or k == best_j:
                continue
            d_new = (dists[best_i][k] * si + dists[best_j][k] * sj) / new_size
            new_dists_row.append(d_new)

        # Rebuild clusters and matrix
        new_clusters = []
        new_heights_list = []
        new_matrix_indices = []
        for k in range(len(clusters)):
            if k == best_i or k == best_j:
                continue
            new_clusters.append(clusters[k])
            new_heights_list.append(heights[k])
            new_matrix_indices.append(k)

        new_clusters.append((new_node, new_size))
        new_heights_list.append(new_height)

        # Rebuild distance matrix
        new_size_m = len(new_clusters)
        new_matrix = [[0.0] * new_size_m for _ in range(new_size_m)]
        for a in range(new_size_m - 1):
            for b in range(a+1, new_size_m - 1):
                oa = new_matrix_indices[a]
                ob = new_matrix_indices[b]
                new_matrix[a][b] = dists[oa][ob]
                new_matrix[b][a] = dists[oa][ob]
            # Distance to new merged cluster
            new_matrix[a][new_size_m - 1] = new_dists_row[a]
            new_matrix[new_size_m - 1][a] = new_dists_row[a]

        clusters = new_clusters
        heights = new_heights_list
        dists = new_matrix

    newick = clusters[0][0] + ";"
    return newick


def build_all_trees():
    results = {}
    pop_labels = [POPULATIONS[p] for p in POP_ORDER]

    for gene, meta in GENES.items():
        filepath = os.path.join(SCRATCH_DIR, meta["file"])
        if not os.path.exists(filepath):
            print(f"  WARNING: {filepath} not found, skipping {gene}")
            continue

        print(f"Processing {gene}...")
        variants = load_variants(filepath)
        print(f"  Loaded {len(variants)} variants")

        pop_afs = compute_pop_af_vectors(variants)
        print(f"  Population AFs: { {p: f'{v:.5f}' for p,v in pop_afs.items()} }")

        matrix = compute_distance_matrix(pop_afs)
        newick = upgma(matrix, pop_labels)

        # Also build distance matrix for display (rounded)
        dist_display = []
        for i in range(len(POP_ORDER)):
            row = []
            for j in range(len(POP_ORDER)):
                row.append(round(matrix[i][j], 6))
            dist_display.append(row)

        results[gene] = {
            "newick": newick,
            "description": meta["description"],
            "populations": pop_labels,
            "distanceMatrix": dist_display,
            "popAFs": {POPULATIONS[p]: round(pop_afs[p], 6) for p in POP_ORDER},
        }
        print(f"  Newick: {newick[:80]}...")

    return results


if __name__ == "__main__":
    trees = build_all_trees()

    # Write data.js
    out_path = os.path.join(OUT_DIR, "data.js")
    js_data = json.dumps(trees, indent=2)
    with open(out_path, "w") as f:
        f.write(f"// Auto-generated from gnomAD population allele frequency data\n")
        f.write(f"// gnomAD v4 | GRCh38 | Data used under https://gnomad.broadinstitute.org/policies\n\n")
        f.write(f"const PANGENOME_DATA = {js_data};\n")

    print(f"\nWrote {out_path}")
    print("Genes processed:", list(trees.keys()))
