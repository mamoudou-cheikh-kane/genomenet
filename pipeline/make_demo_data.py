"""Génère un jeu de données SIMULÉ (data/demo/) pour tester l'application
sans lancer le pipeline. Les valeurs n'ont aucune valeur biologique."""
import csv
import random
from pathlib import Path

rng = random.Random(7)
OUT = Path(__file__).resolve().parent.parent / "data" / "demo"

# Lignées fictives regroupées en complexes clonaux (distances inter-lignées plus faibles)
LINEAGES = {
    "ST258": "CG258", "ST11": "CG258", "ST512": "CG258",
    "ST147": "CG147", "ST392": "CG147",
    "ST307": "CG307", "ST15": "CG15", "ST14": "CG14",
    "ST23": "CG23", "ST101": "CG101", "ST45": "CG45", "ST17": "CG17",
}
COUNTRIES = ["Sénégal", "France", "Italie", "États-Unis", "Inde", "Chine",
             "Brésil", "Nigeria", "Grèce", "Afrique du Sud", "Thaïlande", "Espagne"]
CARBA = ["blaKPC-2", "blaKPC-3", "blaNDM-1", "blaNDM-5", "blaOXA-48", "blaOXA-232"]
OTHER_AMR = ["blaSHV-11", "blaTEM-1", "aac(6')-Ib-cr", "aph(3'')-Ib", "sul1", "sul2",
             "dfrA12", "qnrS1", "tet(A)", "catA1", "fosA", "oqxA", "oqxB", "mph(A)"]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    nodes = []
    idx = 0
    for st, cg in LINEAGES.items():
        n_outbreaks = rng.randint(1, 3)
        for outbreak in range(n_outbreaks):
            country = rng.choice(COUNTRIES)
            base_year = rng.randint(2008, 2023)
            carba_gene = rng.choice(CARBA) if rng.random() < 0.5 else None
            for _ in range(rng.randint(4, 14)):
                idx += 1
                genes = rng.sample(OTHER_AMR, rng.randint(2, 9)) + ["fosA", "oqxA"]
                if carba_gene and rng.random() < 0.85:
                    genes.append(carba_gene)
                esbl = rng.random() < 0.55
                if esbl:
                    genes.append(rng.choice(["blaCTX-M-15", "blaCTX-M-14", "blaCTX-M-27"]))
                genes = sorted(set(genes))
                nodes.append({
                    "id": f"DEMO_{idx:04d}",
                    "strain": f"souche-simulée-{idx}",
                    "country": country if rng.random() < 0.9 else "",
                    "year": base_year + rng.randint(0, 2),
                    "host": "Homo sapiens",
                    "isolation_source": rng.choice(["sang", "urine", "écouvillon rectal", "expectoration"]),
                    "gc_percent": round(rng.gauss(57.3, 0.15), 2),
                    "length_mb": round(rng.gauss(5.55, 0.15), 3),
                    "st": st,
                    "n_amr": len(genes),
                    "carbapenemase": int(any(g in CARBA for g in genes)),
                    "esbl": int(esbl),
                    "amr_genes": ";".join(genes),
                    "_cg": cg, "_outbreak": f"{st}-{outbreak}",
                })

    edges = []
    for i, a in enumerate(nodes):
        for b in nodes[i + 1:]:
            if a["_outbreak"] == b["_outbreak"]:
                d = rng.uniform(0.00005, 0.0012)
            elif a["st"] == b["st"]:
                d = rng.uniform(0.0012, 0.003)
            elif a["_cg"] == b["_cg"]:
                d = rng.uniform(0.003, 0.0065)
            else:
                d = rng.uniform(0.008, 0.02)
            if d <= 0.01:
                edges.append({"source": a["id"], "target": b["id"], "distance": round(d, 6)})

    fields = [k for k in nodes[0] if not k.startswith("_")]
    with open(OUT / "nodes.csv", "w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(nodes)
    with open(OUT / "edges.csv", "w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=["source", "target", "distance"])
        w.writeheader()
        w.writerows(edges)
    print(f"Données simulées : {len(nodes)} génomes, {len(edges)} arêtes dans {OUT}")


if __name__ == "__main__":
    main()
