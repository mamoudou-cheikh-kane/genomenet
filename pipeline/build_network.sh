#!/usr/bin/env bash
# Construit data/nodes.csv et data/edges.csv à partir de génomes NCBI.
# Usage : bash pipeline/build_network.sh
# Variables modifiables : SPECIES, N, THREADS, MAX_DIST, SEED
set -euo pipefail

SPECIES="${SPECIES:-Klebsiella pneumoniae}"
AMR_ORGANISM="${AMR_ORGANISM:-Klebsiella_pneumoniae}"
N="${N:-250}"
THREADS="${THREADS:-4}"
MAX_DIST="${MAX_DIST:-0.01}"
SEED="${SEED:-42}"
WORK="work"

cd "$(dirname "$0")/.."
mkdir -p "$WORK/fna" "$WORK/amr" data

echo "[1/6] Liste des assemblages complets RefSeq de ${SPECIES}"
datasets summary genome taxon "$SPECIES" \
  --assembly-level complete --assembly-source refseq --as-json-lines \
  > "$WORK/summary.jsonl"

echo "[2/6] Tirage de ${N} génomes (graine ${SEED})"
python pipeline/select_genomes.py "$WORK/summary.jsonl" --n "$N" --seed "$SEED" \
  --out "$WORK/accessions.txt"

echo "[3/6] Téléchargement des séquences"
datasets download genome accession --inputfile "$WORK/accessions.txt" \
  --include genome --filename "$WORK/genomes.zip"
unzip -o -q "$WORK/genomes.zip" -d "$WORK"
for f in "$WORK"/ncbi_dataset/data/*/*.fna; do
  acc="$(basename "$(dirname "$f")")"
  cp "$f" "$WORK/fna/${acc}.fna"
done

echo "[4/6] Distances Mash"
mash sketch -p "$THREADS" -s 10000 -o "$WORK/genomes" "$WORK"/fna/*.fna
mash dist -p "$THREADS" "$WORK/genomes.msh" "$WORK/genomes.msh" > "$WORK/mash_dist.tsv"

echo "[5/6] Typage MLST et gènes de résistance (optionnels)"
MLST_ARG=()
if command -v mlst >/dev/null 2>&1; then
  mlst --threads "$THREADS" "$WORK"/fna/*.fna > "$WORK/mlst.tsv"
  MLST_ARG=(--mlst "$WORK/mlst.tsv")
else
  echo "  mlst absent : étape ignorée"
fi

AMR_ARG=()
if command -v amrfinder >/dev/null 2>&1; then
  for f in "$WORK"/fna/*.fna; do
    acc="$(basename "$f" .fna)"
    [ -s "$WORK/amr/${acc}.tsv" ] && continue
    amrfinder -n "$f" --organism "$AMR_ORGANISM" --threads "$THREADS" \
      -o "$WORK/amr/${acc}.tsv" >/dev/null 2>&1
  done
  AMR_ARG=(--amr-dir "$WORK/amr")
else
  echo "  amrfinder absent : étape ignorée"
fi

echo "[6/6] Export CSV"
python pipeline/to_csv.py --summary "$WORK/summary.jsonl" \
  --accessions "$WORK/accessions.txt" --mash "$WORK/mash_dist.tsv" \
  ${MLST_ARG[@]+"${MLST_ARG[@]}"} ${AMR_ARG[@]+"${AMR_ARG[@]}"} --max-dist "$MAX_DIST" --outdir data

echo "Terminé. Lancez : python app.py"
