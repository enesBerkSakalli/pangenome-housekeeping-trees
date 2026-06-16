"""
Fetch NCBI Ortholog protein sets for the three housekeeping genes.

The output mirrors the Ensembl homology JSON shape already consumed by
plotly_stacked_trees.py, but the source is NCBI Datasets + local MAFFT
protein alignments.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import urllib.request
import zipfile
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "outputs_3d"
NCBI_DIR = OUT_DIR / "ncbi_orthologs"
DATASETS = ROOT / "external" / "DataDiVR_WebApp" / "venv" / "bin" / "datasets"

GENE_IDS = {"GAPDH": "2597", "ACTB": "60", "ENO1": "2023", "RPLP0": "6175"}
GENE_ORDER = ["GAPDH", "ENO1", "RPLP0"]
HUMAN_TAX_ID = "9606"

CLADE_BY_LINEAGE = [
    ("Great apes", {"9604"}),
    ("Lesser apes", {"9577"}),
    ("Old World monkeys", {"9527"}),
    ("New World monkeys", {"9479"}),
    ("Strepsirrhines", {"376911", "9445"}),
    ("Tarsiiformes", {"9476"}),
    ("Scandentia", {"9392"}),
    ("Glires", {"314147"}),
    ("Carnivores", {"33554"}),
    ("Cetartiodactyla", {"91561"}),
    ("Perissodactyla", {"9787"}),
    ("Bats", {"9397"}),
    ("Eulipotyphla", {"9362"}),
    ("Afrotheria", {"311790"}),
    ("Xenarthra", {"9348"}),
    ("Monotremata", {"9255"}),
    ("Marsupialia", {"9263"}),
    ("Birds", {"8782"}),
    ("Reptiles", {"8457"}),
    ("Amphibians", {"8292"}),
    ("Lobe-finned fish", {"118072"}),
    ("Ray-finned fish", {"7898"}),
    ("Cartilaginous fish", {"7777"}),
    ("Jawless fish", {"1476529", "7745", "7744"}),
]


def run(command: list[str], cwd: Path | None = None) -> None:
    subprocess.run(command, cwd=cwd, check=True)


def slugify_species(name: str, tax_id: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
    return slug or f"taxon_{tax_id}"


def read_jsonl(path: Path) -> list[dict]:
    records = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    return records


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current = ""
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                sequences[current] = []
            elif current:
                sequences[current].append(line)
    return {key: "".join(parts) for key, parts in sequences.items()}


def write_fasta(records: list[dict], path: Path) -> None:
    with path.open("w") as fh:
        for record in records:
            fh.write(f">{record['node_id']}\n")
            sequence = record["sequence"]
            for index in range(0, len(sequence), 80):
                fh.write(sequence[index : index + 80] + "\n")


def datasets_download(gene: str, refresh: bool) -> Path:
    NCBI_DIR.mkdir(parents=True, exist_ok=True)
    zip_path = NCBI_DIR / f"{gene.lower()}_orthologs.zip"
    extract_dir = NCBI_DIR / gene.lower()
    if refresh:
        zip_path.unlink(missing_ok=True)
        if extract_dir.exists():
            shutil.rmtree(extract_dir)
    if not zip_path.exists():
        run(
            [
                str(DATASETS),
                "download",
                "gene",
                "symbol",
                gene,
                "--taxon",
                "human",
                "--ortholog",
                "vertebrates",
                "--include",
                "protein,product-report",
                "--filename",
                str(zip_path),
                "--no-progressbar",
            ]
        )
    if not extract_dir.exists():
        with zipfile.ZipFile(zip_path) as archive:
            archive.extractall(extract_dir)
    return extract_dir / "ncbi_dataset" / "data"


def best_protein_accession(record: dict, fasta_sequences: dict[str, str], human_length: int | None) -> str | None:
    candidates = []
    for transcript in record.get("transcripts") or []:
        protein = transcript.get("protein") or {}
        accession = protein.get("accessionVersion")
        if not accession or accession not in fasta_sequences:
            continue
        length = int(protein.get("length") or len(fasta_sequences[accession]))
        length_penalty = abs(length - human_length) if human_length else 0
        candidates.append((length_penalty, -length, accession))
    if not candidates:
        return None
    candidates.sort()
    return candidates[0][2]


def taxon_metadata(tax_ids: list[str]) -> dict[str, dict]:
    result = {}
    for index in range(0, len(tax_ids), 120):
        chunk = tax_ids[index : index + 120]
        url = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/" + ",".join(chunk)
        with urllib.request.urlopen(url, timeout=45) as response:
            payload = json.load(response)
        for node in payload.get("taxonomy_nodes") or []:
            taxonomy = node.get("taxonomy") or {}
            tax_id = str(taxonomy.get("tax_id") or "")
            if tax_id:
                result[tax_id] = taxonomy
    return result


def classify_clade(tax_id: str, taxonomy: dict) -> str:
    lineage = {str(item) for item in taxonomy.get("lineage") or []}
    lineage.add(str(tax_id))
    for clade, ids in CLADE_BY_LINEAGE:
        if lineage & ids:
            return clade
    return "Species ortholog context"


def choose_records(gene: str, data_dir: Path) -> list[dict]:
    product_records = read_jsonl(data_dir / "product_report.jsonl")
    fasta_sequences = read_fasta(data_dir / "protein.faa")
    human_record = next((item for item in product_records if str(item.get("taxId")) == HUMAN_TAX_ID), None)
    if not human_record:
        raise RuntimeError(f"{gene}: no human record in NCBI product report")
    human_accession = best_protein_accession(human_record, fasta_sequences, None)
    if not human_accession:
        raise RuntimeError(f"{gene}: no human protein sequence in NCBI FASTA")
    human_length = len(fasta_sequences[human_accession])

    taxids = sorted({str(record.get("taxId")) for record in product_records if record.get("taxId")})
    taxonomy = taxon_metadata(taxids)
    by_taxid: dict[str, dict] = {}
    for record in product_records:
        tax_id = str(record.get("taxId") or "")
        if not tax_id:
            continue
        accession = best_protein_accession(record, fasta_sequences, human_length)
        if not accession:
            continue
        sequence = fasta_sequences[accession]
        taxname = str(record.get("taxname") or taxonomy.get(tax_id, {}).get("organism_name") or tax_id)
        item = {
            "gene": gene,
            "gene_id": str(record.get("geneId") or ""),
            "symbol": str(record.get("symbol") or gene),
            "tax_id": tax_id,
            "taxname": taxname,
            "common_name": str(record.get("commonName") or taxonomy.get(tax_id, {}).get("common_name") or ""),
            "species": slugify_species(taxname, tax_id),
            "label": taxname,
            "protein_id": accession,
            "sequence": sequence,
            "sequence_length": len(sequence),
            "clade": classify_clade(tax_id, taxonomy.get(tax_id, {})),
            "node_id": f"{gene}|{tax_id}|{accession}",
        }
        current = by_taxid.get(tax_id)
        score = (abs(len(sequence) - human_length), -len(sequence), item["protein_id"])
        current_score = (
            abs(len(current["sequence"]) - human_length),
            -len(current["sequence"]),
            current["protein_id"],
        ) if current else None
        if current is None or score < current_score:
            by_taxid[tax_id] = item
    return list(by_taxid.values())


def run_mafft(gene: str, records: list[dict]) -> dict[str, str]:
    input_path = NCBI_DIR / f"{gene.lower()}_selected.faa"
    aligned_path = NCBI_DIR / f"{gene.lower()}_selected.aligned.faa"
    write_fasta(records, input_path)
    mafft = shutil.which("mafft")
    if not mafft:
        raise RuntimeError("MAFFT is required to align NCBI ortholog proteins")
    with aligned_path.open("w") as out, (NCBI_DIR / f"{gene.lower()}_mafft.log").open("w") as err:
        subprocess.run([mafft, "--auto", "--quiet", str(input_path)], stdout=out, stderr=err, check=True)
    return read_fasta(aligned_path)


def identity_to_human(human_aligned: str, target_aligned: str) -> tuple[float, float]:
    comparable = 0
    identical = 0
    positive = 0
    for left, right in zip(human_aligned, target_aligned):
        if left == "-" or right == "-":
            continue
        comparable += 1
        if left == right:
            identical += 1
            positive += 1
    if comparable == 0:
        return 0.0, 0.0
    value = 100.0 * identical / comparable
    return value, 100.0 * positive / comparable


def write_homology_payload(gene: str, records: list[dict], aligned: dict[str, str]) -> dict:
    human = next(record for record in records if record["tax_id"] == HUMAN_TAX_ID)
    human_aligned = aligned[human["node_id"]]
    homologies = []
    for record in records:
        if record["tax_id"] == HUMAN_TAX_ID:
            continue
        target_aligned = aligned.get(record["node_id"])
        if not target_aligned:
            continue
        perc_id, perc_pos = identity_to_human(human_aligned, target_aligned)
        homologies.append(
            {
                "type": "ncbi_ortholog",
                "method_link_type": "NCBI_ORTHOLOG",
                "taxonomy_level": record["clade"],
                "source": {
                    "species": "homo_sapiens",
                    "id": GENE_IDS[gene],
                    "protein_id": human["protein_id"],
                    "taxon_id": int(HUMAN_TAX_ID),
                    "align_seq": human_aligned,
                    "perc_id": 100.0,
                    "perc_pos": 100.0,
                },
                "target": {
                    "species": record["species"],
                    "taxon_id": int(record["tax_id"]),
                    "id": record["gene_id"],
                    "protein_id": record["protein_id"],
                    "align_seq": target_aligned,
                    "perc_id": perc_id,
                    "perc_pos": perc_pos,
                    "label": record["label"],
                    "common_name": record["common_name"],
                    "clade": record["clade"],
                    "ncbi_taxname": record["taxname"],
                    "sequence_length": record["sequence_length"],
                },
            }
        )
    payload = {
        "source": "NCBI Datasets CLI ortholog vertebrates + MAFFT protein alignment",
        "source_url": "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genes/download-ortholog-data-package/",
        "gene": gene,
        "human_gene_id": GENE_IDS[gene],
        "data": [{"id": GENE_IDS[gene], "homologies": homologies}],
    }
    out_path = OUT_DIR / f"ncbi_{gene.lower()}_all_homologies.json"
    with out_path.open("w") as fh:
        json.dump(payload, fh, indent=2)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true", help="redownload NCBI packages")
    args = parser.parse_args()

    summaries = {}
    species_sets = []
    for gene in GENE_ORDER:
        data_dir = datasets_download(gene, args.refresh)
        records = choose_records(gene, data_dir)
        aligned = run_mafft(gene, records)
        payload = write_homology_payload(gene, records, aligned)
        species = {
            homology["target"]["species"]
            for homology in payload["data"][0]["homologies"]
            if homology["target"].get("align_seq")
        }
        species_sets.append(species)
        clades = Counter(homology["target"].get("clade", "Species ortholog context") for homology in payload["data"][0]["homologies"])
        summaries[gene] = {
            "ortholog_records": len(payload["data"][0]["homologies"]),
            "species": len(species),
            "clades": dict(sorted(clades.items())),
        }
        print(f"{gene}: NCBI aligned ortholog species={len(species)}")

    common = set.intersection(*species_sets) if species_sets else set()
    audit = {
        "source": "NCBI Datasets CLI --ortholog vertebrates",
        "alignment": "MAFFT --auto on one selected protein per NCBI taxon per gene",
        "genes": summaries,
        "common_species_across_genes": len(common),
        "common_species": sorted(common),
        "output_files": {
            gene: str(OUT_DIR / f"ncbi_{gene.lower()}_all_homologies.json")
            for gene in GENE_ORDER
        },
    }
    audit_path = OUT_DIR / "ncbi_ortholog_fetch_audit.json"
    with audit_path.open("w") as fh:
        json.dump(audit, fh, indent=2)
    print(f"Common NCBI ortholog species across {', '.join(GENE_ORDER)}: {len(common)}")
    print(f"Wrote {audit_path}")


if __name__ == "__main__":
    main()
