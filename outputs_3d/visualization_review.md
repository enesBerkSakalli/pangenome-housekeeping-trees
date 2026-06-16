# Visualization Review

## Problems in the Previous View

- The main view overloaded one 3D scene with tree topology, population grouping,
  inter-gene Sankey-like flows, within-gene reticulation chords, and literature
  context. That made the plot visually dense and hard to interpret.
- The 3D stack was useful as an explorer, but it was not a good primary
  explanatory figure because perspective and occlusion made the trees and links
  compete with each other.
- The phrase "recombination" was too strong for the current data. The local
  input is population-stratified allele frequencies, not phased haplotypes or
  inferred ancestral recombination graphs. The corrected language is
  "mixed ancestry candidate" or "reticulation candidate."
- The old link layer mixed two different meanings: same-population links between
  gene layers and within-gene secondary-population signals. The redesigned
  dashboard separates those concepts.

## Revised Reading Path

1. Start with gene-level cards: number of variant taxa, population subtrees,
   mixed candidates, and mixed population-pair links.
2. Read the stacked bar chart to compare population-expanded subtree size across
   GAPDH, ACTB, and RPLP0.
3. Read the heatmap matrix to see which dominant-to-secondary population signals
   occur in each gene.
4. Use the layered stacked tree as the primary tree view. Population blocks are
   fixed across genes, intra-population subtrees are planar orthogonal trees,
   and same-population connectors run vertically between gene layers.
5. Use the 3D stacked explorer only as a secondary diagnostic view.

## Current Data Interpretation

- Each gene tree has 256 variant leaves and 8 explicit population subtrees.
- A leaf is assigned to the population where its allele frequency is highest.
- A mixed candidate is flagged when the secondary population allele frequency is
  at least 0.35 of the dominant population allele frequency and at least 0.001.
- These signals suggest allele-frequency sharing or local-ancestry-like
  structure. They do not prove recombination without phased haplotypes, local
  ancestry calls, LD blocks, or ARG inference.

## Coordinate Debug Result

- All node coordinates are finite.
- All reticulation curves are inside the Plotly scene bounds.
- Every reticulation curve endpoint maps back to its intended population subtree
  root with zero x/y endpoint error in the coordinate audit.
- The main layered tree layout now has a separate planar interval audit. It
  checks that every subtree occupies a contiguous leaf interval and that child
  intervals do not overlap. The current audit reports no violations.

See `mixed_ancestry_link_coordinate_audit.json` for the machine-readable audit.
See `planar_tree_layout_audit.json` for the crossing-free layered-tree audit.
