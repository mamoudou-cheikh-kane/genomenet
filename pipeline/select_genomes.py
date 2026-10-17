"""Tire au sort N assemblages, en privilégiant ceux dont le pays d'origine est renseigné."""
import argparse
import random

from ncbi_meta import parse_record, read_jsonl


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", help="Fichier JSON Lines de `datasets summary genome`")
    parser.add_argument("--n", type=int, default=250)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with_country, without_country = [], []
    seen = set()
    for record in read_jsonl(args.summary):
        meta = parse_record(record)
        acc = meta["id"]
        if not acc or acc in seen:
            continue
        seen.add(acc)
        (with_country if meta["country"] else without_country).append(acc)

    rng = random.Random(args.seed)
    rng.shuffle(with_country)
    rng.shuffle(without_country)
    chosen = sorted((with_country + without_country)[: args.n])

    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write("\n".join(chosen) + "\n")

    print(f"{len(seen)} assemblages disponibles, {len(chosen)} retenus "
          f"({min(len(with_country), args.n)} avec pays renseigné).")


if __name__ == "__main__":
    main()
