"""Génère docs/index.html : version statique de GenomeNet (Cytoscape.js),
avec les données intégrées, publiable gratuitement sur GitHub Pages.

Usage : python pipeline/build_static.py
"""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = ROOT / "docs" / "index.html"
TEMPLATE = Path(__file__).resolve().parent / "static_template.html"
SPECIES = "Klebsiella pneumoniae"


def to_number(value, kind=float):
    try:
        return kind(float(value))
    except (TypeError, ValueError):
        return None


def main():
    folder, demo = (DATA, False) if (DATA / "nodes.csv").exists() else (DATA / "demo", True)
    nodes = []
    with open(folder / "nodes.csv", encoding="utf-8") as handle:
        for r in csv.DictReader(handle):
            nodes.append({
                "id": r["id"], "strain": r.get("strain", ""), "country": r.get("country", ""),
                "year": to_number(r.get("year"), int), "host": r.get("host", ""),
                "source": r.get("isolation_source", ""), "gc": to_number(r.get("gc_percent")),
                "len": to_number(r.get("length_mb")), "st": r.get("st", ""),
                "n_amr": to_number(r.get("n_amr"), int),
                "carba": to_number(r.get("carbapenemase"), int),
                "esbl": to_number(r.get("esbl"), int),
                "genes": [g for g in (r.get("amr_genes") or "").split(";") if g],
            })
    edges = []
    with open(folder / "edges.csv", encoding="utf-8") as handle:
        for r in csv.DictReader(handle):
            edges.append([r["source"], r["target"], float(r["distance"])])

    payload = json.dumps({"nodes": nodes, "edges": edges, "demo": demo, "species": SPECIES},
                         ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    html = TEMPLATE.read_text(encoding="utf-8").replace("__GENOMENET_DATA__", payload)
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(html, encoding="utf-8")
    origin = "simulées" if demo else "réelles"
    print(f"{OUT.relative_to(ROOT)} écrit : {len(nodes)} génomes, {len(edges)} liens (données {origin}).")


if __name__ == "__main__":
    main()
