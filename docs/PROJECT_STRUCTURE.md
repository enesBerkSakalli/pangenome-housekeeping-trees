# Project Structure

## Primary scripts

- `plotly_stacked_trees.py`: Builds the stacked 3D Plotly tree figures. It also
  writes layout audits that verify coordinate finiteness, scene bounds, tree
  structure, edge crossings, and layer separation.
- `datadivr_housekeeping_project.py`: Converts the inferred tree layouts into a
  DataDiVR project with one all-genes scene and one focus scene per gene.
- `variant_3d_analysis.py`: Builds initial 3D variant tree representations.
- `population_expanded_trees.py`: Expands human population summaries into
  consistent tree layers.
- `ncbi_ortholog_sequences.py`: Fetches and prepares NCBI ortholog sequence
  material when source data need to be refreshed.
- `mammal_ortholog_context.py`: Builds mammal/taxon context metadata used by
  the tree layers.

## Generated directories

- `scratch/`: Raw local variant JSON input bundles.
- `networks/`: Early NetworkX tree exports.
- `outputs/`: Matplotlib/layered 2D tree outputs.
- `outputs_reference/`: Reference population tree outputs.
- `outputs_3d/`: Main 3D tree outputs, audits, and intermediate caches.
- `external/DataDiVR_WebApp/static/projects/`: Generated DataDiVR project files.

Most generated directories are ignored by git because they are large and can be
recreated. Selected HTML/audit artifacts are kept for review.

## External dependency

`external/DataDiVR_WebApp` points at
`https://github.com/menchelab/DataDiVR_WebApp`. Treat it as an external project,
not as code owned by this repository. Local preview changes are documented as a
patch under `patches/`.
