# Pangenome Housekeeping-Gene Trees

This repository builds exploratory pangenome/ortholog visualizations for stacked
housekeeping-gene trees. The current figure compares GAPDH, ENO1, and RPLP0 as
taxon-aligned 3D tree layers, with shared ortholog species and a neutral human
subtree in every layer.

## What is included

- `plotly_stacked_trees.py` builds the Plotly 3D HTML views and layout audits.
- `datadivr_housekeeping_project.py` exports the same tree construction into the
  DataDiVR project format.
- `variant_3d_analysis.py`, `population_expanded_trees.py`,
  `reference_population_trees.py`, and related scripts build intermediate
  NetworkX tree data.
- `outputs_3d/plotly_*.html` are the current reviewable Plotly artifacts.
- `outputs_3d/*audit.json` records machine-readable construction checks.
- `outputs_3d/datadivr_coordinate_mappings.json` exposes explicit per-scene
  DataDiVR node coordinates in the normalized 0..1 unit cube.
- `datadivr_project/Pangenome_Housekeeping_Stacked_Trees/` is a ready-to-copy
  DataDiVR project folder. This is the final handoff artifact.
- `external/DataDiVR_WebApp` is an external DataDiVR checkout used for the local
  WebGL/VR preview.

Large local caches, virtual environments, raw variant bundles, and regenerated
NetworkX intermediates are intentionally ignored by git. They can be rebuilt
from the scripts when the corresponding cached source data are available.

## Quickstart

Use Python 3.9 through 3.12 for the generator. In this local workspace,
`/usr/bin/python3` creates the known-good Python 3.9 virtual environment; the
Homebrew `python3` currently points at a newer interpreter and should not be
used for this project unless the pinned scientific packages are refreshed.

```bash
/usr/bin/python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python plotly_stacked_trees.py
python datadivr_housekeeping_project.py
```

Keep this root `.venv` separate from `external/DataDiVR_WebApp/venv`. The root
environment uses the pangenome generator pins such as `numpy==2.0.2`; the
DataDiVR app environment uses its own older runtime pins for the preview server.

The Plotly HTML outputs are written to `outputs_3d/`. The DataDiVR export is
written to `outputs_3d/Pangenome_Housekeeping_Stacked_Trees.json`,
`external/DataDiVR_WebApp/static/projects/Pangenome_Housekeeping_Stacked_Trees/`,
and the portable top-level copy under `datadivr_project/`.

## Review the Plotly figure

Open either of these files in a browser:

- `outputs_3d/plotly_hyperbolic_stacked_trees.html`
- `outputs_3d/plotly_unrooted_3d_variant_trees.html`

The audits to check after regeneration are:

- `outputs_3d/stacked_3d_unrooted_layout_audit.json`
- `outputs_3d/hyperbolic_stacked_tree_layout_audit.json`
- `outputs_3d/Pangenome_Housekeeping_Stacked_Trees_datadivr_audit.json`

Healthy current outputs have finite coordinates, no XY tree-edge crossings, and
valid DataDiVR unit-cube positions.

## Review in DataDiVR

Initialize the DataDiVR checkout, install its dependencies as described in
`external/DataDiVR_WebApp/README.md`, then run:

```bash
python datadivr_housekeeping_project.py
scripts/start_datadivr.sh
```

Open `http://127.0.0.1:3000/preview` and select
`Pangenome_Housekeeping_Stacked_Trees`.

The project starts a data-driven presentation in the WebGL preview:

1. all stacked tree layers visible,
2. GAPDH, ENO1, and RPLP0 layers shown one after another with stage labels,
3. the Ray-finned fish subtree highlighted inside each layer,
4. the full Ray-finned fish subtree and its inter-layer connections restored
   together at the end.

For a ready-made project copy, use:

```text
datadivr_project/Pangenome_Housekeeping_Stacked_Trees/
```

The project includes non-native sidecars for agents and downstream tools:

- `paths.json` and `pfile.json["paths"]` store only explicit paths as numeric
  DataDiVR node IDs, for example `[[376, 5109, 2239], ...]`.
- `path_connections.json` stores the segment pairs and path metadata for those
  numeric paths.
- `ancestor_context_flow` paths connect higher ancestor levels around direct
  MRCA corridors; they are intentionally not labeled as direct MRCAs.
- `coordinate_mappings.json` maps numeric DataDiVR node IDs to source node keys,
  annotations, colors, and per-scene coordinates.
- `analysis/Pangenome_Housekeeping_Stacked_Trees_networkx_scenes.pkl` stores the
  full NetworkX scene objects.
- `analysis/manifest.json` lists the bundled analysis files.

The local workspace includes small preview customizations for node radii and
cleaner legends. They are captured in
`patches/datadivr-preview-node-radius.patch` so they can be reapplied to a fresh
DataDiVR checkout.

## Repository Notes

The project currently has no explicit open-source license selected. Add a
`LICENSE` file before telling others they may reuse or redistribute the code.
