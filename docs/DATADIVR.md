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
- The current project has 4 scenes:
  `00_all_genes_stacked`, `01_GAPDH_focus`, `02_ENO1_focus`, and
  `03_RPLP0_focus`.

## Start the server

```bash
scripts/start_datadivr.sh
```

Then open:

```text
http://127.0.0.1:3000/preview
```

Select `Pangenome_Housekeeping_Stacked_Trees` in the preview UI.

## Fresh DataDiVR checkout notes

The upstream DataDiVR preview does not know about this project sidecar:

```text
static/projects/Pangenome_Housekeeping_Stacked_Trees/nodesizes/*.json
```

Apply `patches/datadivr-preview-node-radius.patch` to preserve the local preview
behavior where internal tree support affects node radius and the legend is
deduplicated by active layout.
