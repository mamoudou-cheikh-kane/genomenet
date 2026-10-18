#!/usr/bin/env bash
# Construit data/nodes.csv et data/edges.csv à partir de génomes NCBI.
# Usage : bash pipeline/build_network.sh
# Variables modifiables : SPECIES, N, THREADS, MAX_DIST, SEED, MLST_SCHEME, MLST_CMD
set -euo pipefail

SPECIES="${SPECIES:-Klebsiella pneumoniae}"
AMR_ORGANISM="${AMR_ORGANISM:-Klebsiella_pneumoniae}"
# Schéma MLST imposé : le choix automatique de mlst peut se tromper d'espèce
MLST_SCHEME="${MLST_SCHEME:-klebsiella}"
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

# La liste des génomes est figée dans data/accessions.txt pour que l'analyse soit
# reproductible : NCBI ajoute des génomes chaque jour, ce qui changerait le tirage.
FROZEN="data/accessions.txt"
if [ -s "$FROZEN" ]; then
  echo "[2/6] Liste figée : $(grep -c . "$FROZEN") génomes lus dans ${FROZEN}"
  echo "      (supprimez ce fichier pour refaire un tirage au sort)"
  cp "$FROZEN" "$WORK/accessions.txt"
else
  echo "[2/6] Tirage de ${N} génomes (graine ${SEED})"
  python pipeline/select_genomes.py "$WORK/summary.jsonl" --n "$N" --seed "$SEED" \
    --out "$WORK/accessions.txt"
  cp "$WORK/accessions.txt" "$FROZEN"
  echo "      Liste enregistrée dans ${FROZEN} (à versionner avec git)"
fi

echo "[3/6] Téléchargement des séquences (reprise automatique en cas de coupure)"
# Un dossier par liste de génomes : relancer le script reprend là où il s'était arrêté
PKG="$WORK/pkg_$(md5sum "$WORK/accessions.txt" | cut -c1-8)"
mkdir -p "$PKG"
datasets download genome accession --inputfile "$WORK/accessions.txt" \
  --include genome --dehydrated --filename "$WORK/dehydrated.zip"
unzip -n -q "$WORK/dehydrated.zip" -d "$PKG"
for attempt in 1 2 3 4 5 6; do
  if datasets rehydrate --directory "$PKG"; then
    break
  fi
  if [ "$attempt" -eq 6 ]; then
    echo "Le téléchargement échoue toujours. Relancez le script plus tard : il reprendra où il s'est arrêté."
    exit 1
  fi
  echo "  Connexion interrompue, nouvelle tentative dans 15 s ($attempt/5)..."
  sleep 15
done

rm -f "$WORK"/fna/*.fna
missing=0
while read -r acc; do
  [ -z "$acc" ] && continue
  f="$(ls "$PKG"/ncbi_dataset/data/"$acc"/*.fna 2>/dev/null | head -n 1 || true)"
  if [ -n "$f" ]; then
    cp "$f" "$WORK/fna/${acc}.fna"
  else
    missing=$((missing + 1))
  fi
done < "$WORK/accessions.txt"
echo "  $(ls "$WORK"/fna/*.fna | wc -l) génomes prêts, ${missing} manquant(s)"

echo "[4/6] Distances Mash"
mash sketch -p "$THREADS" -s 10000 -o "$WORK/genomes" "$WORK"/fna/*.fna 2>/dev/null
mash dist -p "$THREADS" "$WORK/genomes.msh" "$WORK/genomes.msh" > "$WORK/mash_dist.tsv"

echo "[5/6] Typage MLST et gènes de résistance (optionnels)"
# mlst tourne dans son propre environnement conda (voir environment-mlst.yml)
CONDA_BIN="${CONDA_EXE:-$HOME/miniforge3/bin/conda}"
if [ -z "${MLST_CMD:-}" ]; then
  if "$CONDA_BIN" run -n genomenet-mlst mlst --version >/dev/null 2>&1; then
    MLST_CMD="$CONDA_BIN run --no-capture-output -n genomenet-mlst mlst"
  else
    MLST_CMD="mlst"
  fi
fi

MLST_ARG=()
if $MLST_CMD --version >/dev/null 2>&1; then
  if $MLST_CMD --quiet --scheme "$MLST_SCHEME" --threads "$THREADS" "$WORK"/fna/*.fna > "$WORK/mlst.tsv" 2>"$WORK/mlst.log"; then
    MLST_ARG=(--mlst "$WORK/mlst.tsv")
    echo "  mlst : terminé (schéma ${MLST_SCHEME})"
  else
    echo "  mlst a échoué (détails dans $WORK/mlst.log) : étape ignorée"
  fi
else
  echo "  mlst indisponible : étape ignorée (créez l'environnement genomenet-mlst)"
fi

AMR_ARG=()
if command -v amrfinder >/dev/null 2>&1; then
  total=$(ls "$WORK"/fna/*.fna | wc -l)
  i=0
  for f in "$WORK"/fna/*.fna; do
    i=$((i + 1))
    acc="$(basename "$f" .fna)"
    [ -s "$WORK/amr/${acc}.tsv" ] && continue
    echo "  AMRFinderPlus ${i}/${total} : ${acc}"
    amrfinder -n "$f" --organism "$AMR_ORGANISM" --threads "$THREADS" \
      -o "$WORK/amr/${acc}.tsv" >/dev/null 2>&1 || echo "    échec pour ${acc}, ignoré"
  done
  AMR_ARG=(--amr-dir "$WORK/amr")
else
  echo "  amrfinder absent : étape ignorée"
fi

echo "[6/6] Export CSV"
python pipeline/to_csv.py --summary "$WORK/summary.jsonl" \
  --accessions "$WORK/accessions.txt" --mash "$WORK/mash_dist.tsv" --fna-dir "$WORK/fna" \
  ${MLST_ARG[@]+"${MLST_ARG[@]}"} ${AMR_ARG[@]+"${AMR_ARG[@]}"} \
  --max-dist "$MAX_DIST" --outdir data

echo "Terminé. Lancez : python app.py"
