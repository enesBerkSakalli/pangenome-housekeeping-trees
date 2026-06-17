# Pangenome_Housekeeping_Stacked_Trees

Fertiges DataDiVR-Projekt fuer die gestapelten housekeeping-gene trees.

## Inhalt

- `00_all_genes_stacked`: alle drei Gen-Layer gemeinsam.
- `01_GAPDH_focus`: Fokus auf GAPDH.
- `02_ENO1_focus`: Fokus auf ENO1.
- `03_RPLP0_focus`: Fokus auf RPLP0.

Das Projekt enthaelt 5.345 Nodes und 5.460 Edges. Die Tree-Layer nutzen GAPDH,
ENO1 und RPLP0 mit denselben 700 Ortholog-Spezies plus neutralem Human-Subtree.

## Nutzung

Lege diesen Ordner unter `DataDiVR_WebApp/static/projects/` ab, starte den
DataDiVR Flask-Server und oeffne:

```text
http://127.0.0.1:3000/preview
```

Dann im Preview-UI `Pangenome_Housekeeping_Stacked_Trees` auswaehlen.

## Validierung

Der Export wurde mit `datadivr_housekeeping_project.py` erzeugt. Das zugehoerige
Audit liegt hier:

```text
outputs_3d/Pangenome_Housekeeping_Stacked_Trees_datadivr_audit.json
```

Der aktuelle Audit-Status: keine fehlenden Dateien, finite Positionen,
Positionen im Unit-Cube.
