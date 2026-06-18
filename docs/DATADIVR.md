# DataDiVR Local Review

## Build the project export

From the repository root:

```bash
source .venv/bin/activate
python datadivr_housekeeping_project.py
```

Expected audit highlights:

- `missing_files` is empty.
- `finite_positions` is `true`.
- `positions_in_unit_cube` is `true`.
- The current project has 8 pre-rendered DataDiVR scenes:
  `00_all_tree_layers`, `01_GAPDH_tree_layer`, `02_ENO1_tree_layer`,
  `03_RPLP0_tree_layer`, `04_GAPDH_ray_finned_fish_subtree`,
  `05_ENO1_ray_finned_fish_subtree`, `06_RPLP0_ray_finned_fish_subtree`, and
  `07_ray_finned_fish_subtree_connections`.

## Start the server

```bash
scripts/start_datadivr.sh
```

Then open:

```text
http://127.0.0.1:3000/preview
```

Select `Pangenome_Housekeeping_Stacked_Trees` in the preview UI.

The generated `pfile.json` includes `subtreeHighlightAnimation.presentation`.
Each step points to a pre-rendered DataDiVR scene, so the node and link colors
are already compiled into scene-specific texture files:

1. `All tree layers`
2. `GAPDH tree layer`
3. `ENO1 tree layer`
4. `RPLP0 tree layer`
5. `Ray-finned fish subtree in GAPDH`
6. `Ray-finned fish subtree in ENO1`
7. `Ray-finned fish subtree in RPLP0`
8. `Ray-finned fish subtree and inter-layer connections`

## Fresh DataDiVR checkout notes

The upstream DataDiVR preview does not know about this project sidecar:

```text
static/projects/Pangenome_Housekeeping_Stacked_Trees/nodesizes/*.json
```

Apply `patches/datadivr-preview-node-radius.patch` to preserve the local preview
behavior where internal tree support affects node radius, the legend is
deduplicated by active layout, explicit path sidecars are loaded, and the
pre-rendered subtree presentation is played from `pfile.json`.
