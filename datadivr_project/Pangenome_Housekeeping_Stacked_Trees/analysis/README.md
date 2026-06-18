# DataDiVR analysis bundle

This directory keeps the non-native analysis sidecars inside the DataDiVR project. DataDiVR renders the texture and JSON files in the project root; scripts and agents should use these sidecars for path semantics, coordinate lookup, and NetworkX-level inspection.

- `../paths.json` and `Pangenome_Housekeeping_Stacked_Trees_paths.json` contain only numeric DataDiVR node-id paths.
- `../path_connections.json` and `Pangenome_Housekeeping_Stacked_Trees_path_connections.json` contain segment pairs and path metadata.
- `datadivr_coordinate_mappings.json` maps those ids back to node keys and per-scene coordinates.
- `Pangenome_Housekeeping_Stacked_Trees_networkx_scenes.pkl` contains the full NetworkX scene objects.
