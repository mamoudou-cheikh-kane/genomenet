"""Assemble métadonnées NCBI, distances Mash, MLST et AMRFinderPlus en deux CSV.

Sorties :
  nodes.csv : un génome par ligne
  edges.csv : paires de génomes dont la distance Mash est <= --max-dist
"""
import argparse
import csv
from pathlib import Path

from ncbi_meta import parse_record, read_jsonl

NODE_FIELDS = [
    "id", "strain", "country", "year", "host", "isolation_source",
    "gc_percent", "length_mb", "st", "n_amr", "carbapenemase", "esbl", "amr_genes",
]


def accession_from_path(path):
    return Path(path).name.removesuffix(".fna")


def read_mlst(path):
    """Sortie de `mlst` (T. Seemann) : fichier, schéma, ST, allèles..."""
    st = {}
    if not path or not Path(path).exists():
        return st
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            value = parts[2].strip()
            st[accession_from_path(parts[0])] = f"ST{value}" if value.isdigit() else "ST inconnu"
    return st


def fasta_stats(path):
    """Taille (Mb) et taux de GC (%) calculés directement à partir de la séquence."""
    length = gc = 0
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            seq = line.strip().upper()
            length += len(seq)
            gc += seq.count("G") + seq.count("C")
    if not length:
        return None, None
    return round(length / 1e6, 3), round(100 * gc / length, 2)


def read_amr(amr_dir, accession):
    """Sortie tabulée d'AMRFinderPlus (colonnes v3 ou v4)."""
    path = Path(amr_dir) / f"{accession}.tsv" if amr_dir else None
    if not path or not path.exists():
        return None
    genes, carbapenemase, esbl = [], False, False
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            symbol = row.get("Element symbol") or row.get("Gene symbol") or ""
            element_type = row.get("Type") or row.get("Element type") or "AMR"
            if not symbol or element_type != "AMR":
                continue
            genes.append(symbol)
            subclass = (row.get("Subclass") or "").upper()
            if "CARBAPENEM" in subclass:
                carbapenemase = True
            if symbol.startswith("blaCTX-M"):
                esbl = True
    return {
        "n_amr": len(genes),
        "carbapenemase": int(carbapenemase),
        "esbl": int(esbl),
        "amr_genes": ";".join(sorted(set(genes))),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--accessions", required=True)
    parser.add_argument("--mash", required=True, help="Sortie de `mash dist`")
    parser.add_argument("--mlst", default=None)
    parser.add_argument("--amr-dir", default=None)
    parser.add_argument("--fna-dir", default=None,
                        help="Dossier des génomes, pour calculer taille et GC absents chez NCBI")
    parser.add_argument("--max-dist", type=float, default=0.01)
    parser.add_argument("--outdir", default="data")
    args = parser.parse_args()

    wanted = {line.strip() for line in open(args.accessions, encoding="utf-8") if line.strip()}
    meta = {}
    for record in read_jsonl(args.summary):
        parsed = parse_record(record)
        if parsed["id"] in wanted:
            meta[parsed["id"]] = parsed

    if args.fna_dir:
        computed = 0
        for acc, row in meta.items():
            fna = Path(args.fna_dir) / f"{acc}.fna"
            if (row["gc_percent"] is None or row["length_mb"] is None) and fna.exists():
                length_mb, gc = fasta_stats(fna)
                row["length_mb"] = row["length_mb"] if row["length_mb"] is not None else length_mb
                row["gc_percent"] = row["gc_percent"] if row["gc_percent"] is not None else gc
                computed += 1
        print(f"Taille et GC calculés depuis les séquences pour {computed} génomes")

    st = read_mlst(args.mlst)
    for acc, row in meta.items():
        row["st"] = st.get(acc, "")
        amr = read_amr(args.amr_dir, acc)
        row.update(amr or {"n_amr": "", "carbapenemase": "", "esbl": "", "amr_genes": ""})

    edges, seen = [], set()
    with open(args.mash, encoding="utf-8") as handle:
        for line in handle:
            ref, query, dist = line.split("\t")[:3]
            a, b = accession_from_path(ref), accession_from_path(query)
            if a == b or a not in meta or b not in meta:
                continue
            key = tuple(sorted((a, b)))
            if key in seen:
                continue
            seen.add(key)
            d = float(dist)
            if d <= args.max_dist:
                edges.append({"source": key[0], "target": key[1], "distance": round(d, 6)})

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "nodes.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=NODE_FIELDS)
        writer.writeheader()
        for acc in sorted(meta):
            writer.writerow({k: meta[acc].get(k, "") for k in NODE_FIELDS})
    with open(outdir / "edges.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source", "target", "distance"])
        writer.writeheader()
        writer.writerows(edges)

    print(f"{len(meta)} génomes et {len(edges)} arêtes (distance <= {args.max_dist}) écrits dans {outdir}/")


if __name__ == "__main__":
    main()
