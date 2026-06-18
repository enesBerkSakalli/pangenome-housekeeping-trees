# Agent Guide

This project builds pangenome/ortholog visualizations for stacked housekeeping
gene trees and exports them to Plotly and DataDiVR.

## High-Value Entry Points

- `plotly_stacked_trees.py`: Builds the inferred GAPDH, ENO1, and RPLP0 stacked
  3D trees, Plotly HTML views, and layout audits.
- `datadivr_housekeeping_project.py`: Converts the inferred NetworkX trees into
  a DataDiVR project. This is the main exporter for `pfile.json`, texture files,
  `nodes.json`, `links.json`, `paths.json`, NetworkX scene pickles, and the
  merged JSON export.
- `external/DataDiVR_WebApp/static/js/webGL_preview.js`: Local DataDiVR preview
  customization. It reads node-size sidecars and explicit path metadata, then
  renders path curves and animated pulses.

## Generated Data

- `datadivr_project/Pangenome_Housekeeping_Stacked_Trees/`: Portable, ready to
  copy DataDiVR project folder. Treat this as the final product.
- `external/DataDiVR_WebApp/static/projects/Pangenome_Housekeeping_Stacked_Trees/`:
  Runtime copy used by the local DataDiVR checkout.
- `outputs_3d/Pangenome_Housekeeping_Stacked_Trees.json`: Merged DataDiVR-style
  graph/layout JSON.
- `outputs_3d/Pangenome_Housekeeping_Stacked_Trees_networkx_scenes.pkl`: Four
  NetworkX scene graphs.
- `outputs_3d/Pangenome_Housekeeping_Stacked_Trees_paths.json`: Explicit
  path metadata.
- `outputs_3d/datadivr_coordinate_mappings.json`: Per-scene normalized
  DataDiVR coordinate mappings.
- `datadivr_project/Pangenome_Housekeeping_Stacked_Trees/analysis/`: The same
  non-native sidecars bundled inside the DataDiVR project, including the
  NetworkX scene pickle.

## Path Model

DataDiVR still receives ordinary pairwise links in `links.json` and link
textures. The additional path layer is independent:

- `pfile.json["paths"]` is a list of paths as numeric DataDiVR node IDs:
  `[[1, 2, 3, 4], [4, 6, 7], ...]`.
- `paths.json["path_records"]` adds metadata for each path, including kind,
  gene pair, clade/group, label, and the same numeric `nodes` list.
- Current path kinds are `clade_flow` and `selected_taxon_flow`.
- Do not infer paths by connected components of interlayer edges. Several
  corridors share MRCA nodes, so component inference can merge distinct paths.
  Use `paths.json` as the source of truth.

## Animation Model

`pfile.json["pathAnimationSettings"]` controls the WebGL preview animation:

- `enabled`
- `drawPathCurves`
- `maxVisiblePaths`
- `curveSegments`
- `pulseRadius`
- `pulseSpeed`
- `pulseStagger`
- `curveOpacity`
- `pulseOpacity`
- `focusSceneSpeedBoost`

The preview uses Three.js curves over the listed node IDs and moves pulse
markers along the curve with a normalized 0..1 progress value. Keep animation
settings data-driven in `pfile.json`; avoid hard-coding project-specific values
inside `webGL_preview.js` unless they are stable defaults.

## Validation Commands

Use the project virtualenv:

```bash
.venv/bin/python datadivr_housekeeping_project.py
.venv/bin/python -m py_compile datadivr_housekeeping_project.py
```

Important audit checks:

- `missing_files` is empty.
- `finite_positions` is `true`.
- `positions_in_unit_cube` is `true`.
- `path_count` matches the number of explicit path records.

## Publishing

The public repo is:

```text
https://github.com/enesBerkSakalli/pangenome-housekeeping-trees
```

This workspace may be published from a temporary checkout of the public repo if
the top-level project directory is not itself a Git working tree. Keep the
published scope focused on scripts, docs, `datadivr_project/`, `outputs_3d/`
handoff JSON/pickle files, and `patches/datadivr-preview-node-radius.patch`.

The DataDiVR checkout is external. Local preview changes must be captured in
`patches/datadivr-preview-node-radius.patch` so a fresh DataDiVR checkout can
replay them.
