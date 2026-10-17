"""Lecture des rapports JSON Lines produits par `datasets summary genome`."""
import json
import re

MISSING = {
    "", "missing", "not collected", "not applicable", "not provided",
    "unknown", "na", "n/a", "none", "-", "restricted access",
}


def clean(value):
    value = (value or "").strip()
    return "" if value.lower() in MISSING else value


def read_jsonl(path):
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                yield json.loads(line)


def biosample_attributes(record):
    biosample = (record.get("assemblyInfo") or {}).get("biosample") or {}
    attrs = {}
    for attr in biosample.get("attributes") or []:
        name = attr.get("name")
        if name:
            attrs[name] = attr.get("value", "")
    return attrs


def parse_record(record):
    """Extrait les champs utiles d'un rapport d'assemblage NCBI."""
    attrs = biosample_attributes(record)
    geo = clean(attrs.get("geo_loc_name") or attrs.get("country"))
    country = geo.split(":")[0].strip() if geo else ""

    date = clean(attrs.get("collection_date"))
    match = re.search(r"(19|20)\d{2}", date)
    year = int(match.group(0)) if match else None

    organism = record.get("organism") or {}
    strain = (organism.get("infraspecificNames") or {}).get("strain") or attrs.get("strain") or ""

    stats = record.get("assemblyStats") or {}
    gc = stats.get("gcPercent")
    length = stats.get("totalSequenceLength")

    return {
        "id": record.get("accession"),
        "strain": clean(strain),
        "country": country,
        "year": year,
        "host": clean(attrs.get("host")),
        "isolation_source": clean(attrs.get("isolation_source")),
        "gc_percent": float(gc) if gc not in (None, "") else None,
        "length_mb": round(int(length) / 1e6, 3) if length not in (None, "") else None,
    }
