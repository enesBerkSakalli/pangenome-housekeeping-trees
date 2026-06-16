"""
Fetch compact mammalian ortholog context for the housekeeping genes.

The output is deliberately metadata-only: it stores orthology type, taxonomy
level, Ensembl IDs, and protein percent identity. These records are used as
comparative-genomics leaves attached to the human variant trees.
"""

from __future__ import annotations

import json
import os
import socket
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT_PATH = os.path.join(ROOT, "outputs_3d", "mammal_ortholog_context.json")

GENE_ORDER = ["GAPDH", "ACTB", "RPLP0"]

MAMMAL_SPECIES = [
    {"species": "pan_troglodytes", "label": "Chimpanzee", "clade": "Great apes"},
    {"species": "pan_paniscus", "label": "Bonobo", "clade": "Great apes"},
    {"species": "gorilla_gorilla", "label": "Gorilla", "clade": "Great apes"},
    {"species": "pongo_abelii", "label": "Sumatran orangutan", "clade": "Great apes"},
    {"species": "pongo_pygmaeus", "label": "Bornean orangutan", "clade": "Great apes"},
    {"species": "nomascus_leucogenys", "label": "Gibbon", "clade": "Lesser apes"},
    {"species": "macaca_mulatta", "label": "Macaque", "clade": "Old World monkeys"},
    {"species": "macaca_fascicularis", "label": "Crab-eating macaque", "clade": "Old World monkeys"},
    {"species": "macaca_nemestrina", "label": "Pig-tailed macaque", "clade": "Old World monkeys"},
    {"species": "papio_anubis", "label": "Baboon", "clade": "Old World monkeys"},
    {"species": "cercocebus_atys", "label": "Sooty mangabey", "clade": "Old World monkeys"},
    {"species": "mandrillus_leucophaeus", "label": "Drill", "clade": "Old World monkeys"},
    {"species": "rhinopithecus_bieti", "label": "Black snub-nosed monkey", "clade": "Old World monkeys"},
    {"species": "rhinopithecus_roxellana", "label": "Golden snub-nosed monkey", "clade": "Old World monkeys"},
    {"species": "chlorocebus_sabaeus", "label": "Vervet", "clade": "Old World monkeys"},
    {"species": "callithrix_jacchus", "label": "Marmoset", "clade": "New World monkeys"},
    {"species": "saimiri_boliviensis_boliviensis", "label": "Squirrel monkey", "clade": "New World monkeys"},
    {"species": "aotus_nancymaae", "label": "Owl monkey", "clade": "New World monkeys"},
    {"species": "cebus_capucinus", "label": "Capuchin", "clade": "New World monkeys"},
    {"species": "cebus_imitator", "label": "Panamanian capuchin", "clade": "New World monkeys"},
    {"species": "sapajus_apella", "label": "Tufted capuchin", "clade": "New World monkeys"},
    {"species": "microcebus_murinus", "label": "Mouse lemur", "clade": "Strepsirrhines"},
    {"species": "propithecus_coquereli", "label": "Sifaka", "clade": "Strepsirrhines"},
    {"species": "otolemur_garnettii", "label": "Bushbaby", "clade": "Strepsirrhines"},
    {"species": "prolemur_simus", "label": "Greater bamboo lemur", "clade": "Strepsirrhines"},
    {"species": "carlito_syrichta", "label": "Tarsier", "clade": "Tarsiiformes"},
    {"species": "tupaia_belangeri", "label": "Common tree shrew", "clade": "Scandentia"},
    {"species": "mus_musculus", "label": "Mouse", "clade": "Glires"},
    {"species": "mus_spicilegus", "label": "Steppe mouse", "clade": "Glires"},
    {"species": "rattus_norvegicus", "label": "Rat", "clade": "Glires"},
    {"species": "oryctolagus_cuniculus", "label": "Rabbit", "clade": "Glires"},
    {"species": "cavia_porcellus", "label": "Guinea pig", "clade": "Glires"},
    {"species": "octodon_degus", "label": "Degu", "clade": "Glires"},
    {"species": "sciurus_vulgaris", "label": "Red squirrel", "clade": "Glires"},
    {"species": "ochotona_princeps", "label": "Pika", "clade": "Glires"},
    {"species": "mesocricetus_auratus", "label": "Golden hamster", "clade": "Glires"},
    {"species": "cricetulus_griseus_chok1gshd", "label": "Chinese hamster", "clade": "Glires"},
    {"species": "cricetulus_griseus_picr", "label": "Chinese hamster PICR", "clade": "Glires"},
    {"species": "jaculus_jaculus", "label": "Jerboa", "clade": "Glires"},
    {"species": "dipodomys_ordii", "label": "Kangaroo rat", "clade": "Glires"},
    {"species": "heterocephalus_glaber_female", "label": "Naked mole-rat", "clade": "Glires"},
    {"species": "peromyscus_maniculatus_bairdii", "label": "Deer mouse", "clade": "Glires"},
    {"species": "marmota_marmota_marmota", "label": "Alpine marmot", "clade": "Glires"},
    {"species": "castor_canadensis", "label": "Beaver", "clade": "Glires"},
    {"species": "ictidomys_tridecemlineatus", "label": "Thirteen-lined ground squirrel", "clade": "Glires"},
    {"species": "urocitellus_parryii", "label": "Arctic ground squirrel", "clade": "Glires"},
    {"species": "canis_lupus_familiaris", "label": "Dog", "clade": "Carnivores"},
    {"species": "canis_lupus_dingo", "label": "Dingo", "clade": "Carnivores"},
    {"species": "felis_catus", "label": "Cat", "clade": "Carnivores"},
    {"species": "mustela_putorius_furo", "label": "Ferret", "clade": "Carnivores"},
    {"species": "neovison_vison", "label": "American mink", "clade": "Carnivores"},
    {"species": "ursus_maritimus", "label": "Polar bear", "clade": "Carnivores"},
    {"species": "ailuropoda_melanoleuca", "label": "Giant panda", "clade": "Carnivores"},
    {"species": "vulpes_vulpes", "label": "Red fox", "clade": "Carnivores"},
    {"species": "panthera_leo", "label": "Lion", "clade": "Carnivores"},
    {"species": "panthera_pardus", "label": "Leopard", "clade": "Carnivores"},
    {"species": "leptonychotes_weddellii", "label": "Weddell seal", "clade": "Carnivores"},
    {"species": "bos_taurus", "label": "Cow", "clade": "Cetartiodactyla"},
    {"species": "bos_indicus_hybrid", "label": "Hybrid cattle", "clade": "Cetartiodactyla"},
    {"species": "bos_mutus", "label": "Wild yak", "clade": "Cetartiodactyla"},
    {"species": "sus_scrofa", "label": "Pig", "clade": "Cetartiodactyla"},
    {"species": "ovis_aries", "label": "Sheep", "clade": "Cetartiodactyla"},
    {"species": "capra_hircus", "label": "Goat", "clade": "Cetartiodactyla"},
    {"species": "vicugna_pacos", "label": "Alpaca", "clade": "Cetartiodactyla"},
    {"species": "camelus_dromedarius", "label": "Arabian camel", "clade": "Cetartiodactyla"},
    {"species": "camelus_bactrianus", "label": "Bactrian camel", "clade": "Cetartiodactyla"},
    {"species": "moschus_moschiferus", "label": "Siberian musk deer", "clade": "Cetartiodactyla"},
    {"species": "tursiops_truncatus", "label": "Dolphin", "clade": "Cetartiodactyla"},
    {"species": "delphinapterus_leucas", "label": "Beluga whale", "clade": "Cetartiodactyla"},
    {"species": "monodon_monoceros", "label": "Narwhal", "clade": "Cetartiodactyla"},
    {"species": "phocoena_sinus", "label": "Vaquita", "clade": "Cetartiodactyla"},
    {"species": "physeter_catodon", "label": "Sperm whale", "clade": "Cetartiodactyla"},
    {"species": "orcinus_orca", "label": "Killer whale", "clade": "Cetartiodactyla"},
    {"species": "balaenoptera_musculus", "label": "Blue whale", "clade": "Cetartiodactyla"},
    {"species": "bubalus_bubalis", "label": "Water buffalo", "clade": "Cetartiodactyla"},
    {"species": "bison_bison_bison", "label": "Bison", "clade": "Cetartiodactyla"},
    {"species": "pantholops_hodgsonii", "label": "Tibetan antelope", "clade": "Cetartiodactyla"},
    {"species": "equus_caballus", "label": "Horse", "clade": "Perissodactyla"},
    {"species": "equus_asinus", "label": "Donkey", "clade": "Perissodactyla"},
    {"species": "equus_przewalskii", "label": "Przewalski horse", "clade": "Perissodactyla"},
    {"species": "ceratotherium_simum_simum", "label": "White rhinoceros", "clade": "Perissodactyla"},
    {"species": "pteropus_vampyrus", "label": "Megabat", "clade": "Bats"},
    {"species": "pteropus_alecto", "label": "Black flying fox", "clade": "Bats"},
    {"species": "myotis_lucifugus", "label": "Microbat", "clade": "Bats"},
    {"species": "myotis_davidii", "label": "David's myotis", "clade": "Bats"},
    {"species": "myotis_brandtii", "label": "Brandt's bat", "clade": "Bats"},
    {"species": "miniopterus_natalensis", "label": "Natal long-fingered bat", "clade": "Bats"},
    {"species": "rhinolophus_ferrumequinum", "label": "Horseshoe bat", "clade": "Bats"},
    {"species": "erinaceus_europaeus", "label": "Hedgehog", "clade": "Eulipotyphla"},
    {"species": "sorex_araneus", "label": "Shrew", "clade": "Eulipotyphla"},
    {"species": "condylura_cristata", "label": "Star-nosed mole", "clade": "Eulipotyphla"},
    {"species": "talpa_occidentalis", "label": "Iberian mole", "clade": "Eulipotyphla"},
    {"species": "loxodonta_africana", "label": "Elephant", "clade": "Afrotheria"},
    {"species": "echinops_telfairi", "label": "Tenrec", "clade": "Afrotheria"},
    {"species": "procavia_capensis", "label": "Hyrax", "clade": "Afrotheria"},
    {"species": "trichechus_manatus_latirostris", "label": "Manatee", "clade": "Afrotheria"},
    {"species": "orycteropus_afer_afer", "label": "Aardvark", "clade": "Afrotheria"},
    {"species": "chrysochloris_asiatica", "label": "Cape golden mole", "clade": "Afrotheria"},
    {"species": "dasypus_novemcinctus", "label": "Armadillo", "clade": "Xenarthra"},
    {"species": "choloepus_hoffmanni", "label": "Sloth", "clade": "Xenarthra"},
    {"species": "tamandua_tetradactyla", "label": "Anteater", "clade": "Xenarthra"},
    {"species": "ornithorhynchus_anatinus", "label": "Platypus", "clade": "Monotremata"},
    {"species": "tachyglossus_aculeatus", "label": "Echidna", "clade": "Monotremata"},
    {"species": "monodelphis_domestica", "label": "Opossum", "clade": "Marsupialia"},
    {"species": "sarcophilus_harrisii", "label": "Tasmanian devil", "clade": "Marsupialia"},
    {"species": "dasyurus_viverrinus", "label": "Eastern quoll", "clade": "Marsupialia"},
    {"species": "phascolarctos_cinereus", "label": "Koala", "clade": "Marsupialia"},
    {"species": "vombatus_ursinus", "label": "Wombat", "clade": "Marsupialia"},
    {"species": "notamacropus_eugenii", "label": "Wallaby", "clade": "Marsupialia"},
    {"species": "anser_brachyrhynchus", "label": "Pink-footed goose", "clade": "Birds"},
    {"species": "aquila_chrysaetos_chrysaetos", "label": "Golden eagle", "clade": "Birds"},
    {"species": "coturnix_japonica", "label": "Japanese quail", "clade": "Birds"},
    {"species": "ficedula_albicollis", "label": "Collared flycatcher", "clade": "Birds"},
    {"species": "geospiza_fortis", "label": "Medium ground finch", "clade": "Birds"},
    {"species": "strigops_habroptila", "label": "Kakapo", "clade": "Birds"},
    {"species": "struthio_camelus_australis", "label": "Ostrich", "clade": "Birds"},
    {"species": "taeniopygia_guttata", "label": "Zebra finch", "clade": "Birds"},
    {"species": "anolis_carolinensis", "label": "Green anole", "clade": "Reptiles"},
    {"species": "chelonoidis_abingdonii", "label": "Pinta giant tortoise", "clade": "Reptiles"},
    {"species": "chrysemys_picta_bellii", "label": "Painted turtle", "clade": "Reptiles"},
    {"species": "crocodylus_porosus", "label": "Saltwater crocodile", "clade": "Reptiles"},
    {"species": "gopherus_evgoodei", "label": "Goode's thornscrub tortoise", "clade": "Reptiles"},
    {"species": "laticauda_laticaudata", "label": "Blue-banded sea krait", "clade": "Reptiles"},
    {"species": "naja_naja", "label": "Indian cobra", "clade": "Reptiles"},
    {"species": "sphenodon_punctatus", "label": "Tuatara", "clade": "Reptiles"},
    {"species": "xenopus_tropicalis", "label": "Western clawed frog", "clade": "Amphibians"},
    {"species": "leptobrachium_leishanense", "label": "Leishan spiny toad", "clade": "Amphibians"},
    {"species": "danio_rerio", "label": "Zebrafish", "clade": "Ray-finned fish"},
    {"species": "oryzias_latipes", "label": "Medaka", "clade": "Ray-finned fish"},
    {"species": "gasterosteus_aculeatus", "label": "Stickleback", "clade": "Ray-finned fish"},
    {"species": "salmo_salar", "label": "Atlantic salmon", "clade": "Ray-finned fish"},
    {"species": "tetraodon_nigroviridis", "label": "Green spotted puffer", "clade": "Ray-finned fish"},
    {"species": "xiphophorus_maculatus", "label": "Platyfish", "clade": "Ray-finned fish"},
    {"species": "poecilia_reticulata", "label": "Guppy", "clade": "Ray-finned fish"},
    {"species": "astyanax_mexicanus", "label": "Mexican tetra", "clade": "Ray-finned fish"},
    {"species": "ictalurus_punctatus", "label": "Channel catfish", "clade": "Ray-finned fish"},
    {"species": "lepisosteus_oculatus", "label": "Spotted gar", "clade": "Ray-finned fish"},
    {"species": "erpetoichthys_calabaricus", "label": "Reedfish", "clade": "Ray-finned fish"},
    {"species": "latimeria_chalumnae", "label": "Coelacanth", "clade": "Lobe-finned fish"},
    {"species": "petromyzon_marinus", "label": "Sea lamprey", "clade": "Jawless fish"},
]


def fetch_homology(gene: str, species: str) -> dict | None:
    query = urllib.parse.urlencode(
        {
            "type": "orthologues",
            "target_species": species,
            "content-type": "application/json",
        },
        safe=";",
    )
    url = f"https://rest.ensembl.org/homology/symbol/human/{gene}?{query}"
    request = urllib.request.Request(
        url,
        headers={
            "Content-Type": "application/json",
            "User-Agent": "PangenomeCodex/1.0",
        },
    )
    last_error = ""
    for attempt in range(2):
        try:
            with urllib.request.urlopen(request, timeout=14) as response:
                payload = json.load(response)
            break
        except (urllib.error.URLError, TimeoutError, socket.timeout, json.JSONDecodeError, OSError) as exc:
            last_error = str(exc)
            time.sleep(0.4 + attempt)
    else:
        return {"error": last_error, "species": species}

    records = payload.get("data", [])
    if not records:
        return None
    homologies = records[0].get("homologies", [])
    if not homologies:
        return None

    def rank(record: dict) -> tuple[int, float]:
        target = record.get("target") or {}
        one_to_one = 0 if record.get("type") == "ortholog_one2one" else 1
        identity = float(target.get("perc_id") or 0.0)
        return one_to_one, -identity

    return sorted(homologies, key=rank)[0]


def build_context() -> dict:
    prior = {}
    if os.path.exists(OUT_PATH):
        with open(OUT_PATH) as fh:
            prior_payload = json.load(fh)
        prior = {
            gene: {record.get("species"): record for record in records}
            for gene, records in prior_payload.get("genes", {}).items()
        }

    context = {
        "source": "Ensembl REST GET /homology/symbol/:species/:symbol",
        "source_url": "https://rest.ensembl.org/documentation/info/homology_symbol",
        "human_species": "homo_sapiens",
        "genes": {},
    }
    fetch_jobs = []
    for gene in GENE_ORDER:
        gene_records = []
        for species_meta in MAMMAL_SPECIES:
            cached = prior.get(gene, {}).get(species_meta["species"])
            if cached:
                gene_records.append({**cached, **species_meta})
                continue
            record_index = len(gene_records)
            gene_records.append(None)
            fetch_jobs.append((gene, record_index, species_meta))
        context["genes"][gene] = gene_records

    def build_record(gene: str, species_meta: dict) -> tuple[str, str, dict]:
        homology = fetch_homology(gene, species_meta["species"])
        if not homology or homology.get("error"):
            error = homology.get("error") if homology else ""
            return gene, species_meta["species"], {**species_meta, "available": False, "error": error}
        target = homology.get("target") or {}
        return (
            gene,
            species_meta["species"],
            {
                **species_meta,
                "available": True,
                "orthology_type": homology.get("type", ""),
                "taxonomy_level": homology.get("taxonomy_level", ""),
                "ensembl_gene_id": target.get("id", ""),
                "protein_id": target.get("protein_id", ""),
                "taxon_id": target.get("taxon_id", ""),
                "percent_identity": float(target.get("perc_id") or 0.0),
                "percent_positive": float(target.get("perc_pos") or 0.0),
            },
        )

    if fetch_jobs:
        with ThreadPoolExecutor(max_workers=6) as executor:
            future_map = {
                executor.submit(build_record, gene, species_meta): (gene, record_index)
                for gene, record_index, species_meta in fetch_jobs
            }
            for future in as_completed(future_map):
                gene, record_index = future_map[future]
                _gene, _species, record = future.result()
                context["genes"][gene][record_index] = record
    return context


def main() -> None:
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    context = build_context()
    with open(OUT_PATH, "w") as fh:
        json.dump(context, fh, indent=2)
    for gene, records in context["genes"].items():
        available = sum(1 for record in records if record.get("available"))
        print(f"{gene}: species orthologs={available}/{len(records)}")
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
