"""Lecture des rapports JSON Lines produits par `datasets summary genome`."""
import json
import re

MISSING = {
    "", "missing", "not collected", "not applicable", "not provided",
    "unknown", "na", "n/a", "none", "-", "restricted access",
}

# Clés au format camelCase parfois utilisées par NCBI à la place des attributs name/value
CAMEL_KEYS = {
    "geoLocName": "geo_loc_name",
    "collectionDate": "collection_date",
    "isolationSource": "isolation_source",
    "host": "host",
    "strain": "strain",
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


def collect_attributes(obj, attrs):
    """Parcourt tout le rapport et récupère les attributs BioSample où qu'ils soient."""
    if isinstance(obj, dict):
        name, value = obj.get("name"), obj.get("value")
        if isinstance(name, str) and isinstance(value, (str, int, float)):
            if clean(str(value)) or name not in attrs:
                attrs[name] = str(value)
        for key, val in obj.items():
            if key in CAMEL_KEYS and isinstance(val, (str, int, float)):
                attrs.setdefault(CAMEL_KEYS[key], str(val))
            collect_attributes(val, attrs)
    elif isinstance(obj, list):
        for item in obj:
            collect_attributes(item, attrs)
    return attrs


def parse_record(record):
    """Extrait les champs utiles d'un rapport d'assemblage NCBI."""
    attrs = collect_attributes(record, {})
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
