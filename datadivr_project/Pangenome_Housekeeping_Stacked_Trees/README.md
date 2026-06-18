# Pangenome Housekeeping Stacked Trees

This folder is a complete DataDiVR project. Copy the whole `Pangenome_Housekeeping_Stacked_Trees` directory into `DataDiVR_WebApp/static/projects/` and select it from the preview.

Native DataDiVR files are in the project root plus the `layouts`, `layoutsl`, `layoutsRGB`, `links`, `linksRGB`, and `nodesizes` directories.

Additional sidecars:

- `paths.json` stores only explicit paths as numeric DataDiVR node IDs.
- `path_connections.json` stores per-path segment pairs and metadata.
- `coordinate_mappings.json` maps numeric node IDs to node keys, annotations, colors, and coordinates in each scene.
- `pfile.json` includes `subtreeHighlightAnimation`, which starts with the Ray-finned fish subtree highlighted across all layers and then per active layer while non-focused subtrees are dimmed.
- `analysis/` bundles the merged JSON export, NetworkX scene pickle, paths, coordinate mappings, and audit files for downstream agents.
