"""GenomeNet : réseau interactif de similarité de génomes bactériens (Dash + Cytoscape)."""
from pathlib import Path

import dash_cytoscape as cyto
import pandas as pd
from dash import Dash, Input, Output, dcc, html

cyto.load_extra_layouts()

DATA_DIR = Path(__file__).resolve().parent / "data"
SPECIES = "Klebsiella pneumoniae"

# --------------------------------------------------------------------------- données

def load_data():
    real = DATA_DIR / "nodes.csv", DATA_DIR / "edges.csv"
    if all(p.exists() for p in real):
        folder, demo = DATA_DIR, False
    else:
        folder, demo = DATA_DIR / "demo", True
    nodes = pd.read_csv(folder / "nodes.csv", dtype={"id": str})
    edges = pd.read_csv(folder / "edges.csv", dtype={"source": str, "target": str})

    defaults = {"strain": "", "country": "", "year": None, "host": "", "isolation_source": "",
                "gc_percent": None, "length_mb": None, "st": "", "n_amr": None,
                "carbapenemase": None, "esbl": None, "amr_genes": ""}
    for col, value in defaults.items():
        if col not in nodes:
            nodes[col] = value
    for col in ["strain", "country", "host", "isolation_source", "st", "amr_genes"]:
        nodes[col] = nodes[col].fillna("").astype(str)
    return nodes, edges, demo


NODES, EDGES, IS_DEMO = load_data()
NODE_IDS = NODES["id"].tolist()
MAX_DIST = float(EDGES["distance"].max()) if len(EDGES) else 0.01
DEFAULT_THRESHOLD = min(0.0025, MAX_DIST)
HAS_AMR = NODES["n_amr"].notna().any()
HAS_ST = (NODES["st"] != "").any()

PALETTE = ["#1f6f78", "#e0a526", "#b5446e", "#3f7cc4", "#6a9a3a", "#8a5cc2",
           "#d9623b", "#2a9d8f", "#a0522d", "#5c6f91", "#c2a83e", "#7b3f61"]
GREY = "#b9c0c4"

# --------------------------------------------------------------------------- calculs

def connected_components(ids, edge_df):
    parent = {i: i for i in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in zip(edge_df["source"], edge_df["target"]):
        if a in parent and b in parent:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    groups = {}
    for i in ids:
        groups.setdefault(find(i), []).append(i)
    ordered = sorted(groups.values(), key=len, reverse=True)
    label = {}
    rank = 0
    for members in ordered:
        if len(members) == 1:
            label[members[0]] = "Isolé"
        else:
            rank += 1
            for m in members:
                label[m] = f"Groupe {rank}"
    return label, ordered


def categorical_colors(values, top=12):
    counts = pd.Series(values).value_counts()
    keep = [v for v in counts.index if v not in ("", "Inconnu", "Isolé", "ST inconnu")][:top]
    mapping = {v: PALETTE[i % len(PALETTE)] for i, v in enumerate(keep)}
    colors = [mapping.get(v, GREY) for v in values]
    legend = [(v, mapping[v]) for v in keep]
    if any(c == GREY for c in colors):
        legend.append(("Autres ou non renseigné", GREY))
    return colors, legend


def year_colors(years):
    known = [y for y in years if pd.notna(y)]
    if not known:
        return [GREY] * len(years), [("Année non renseignée", GREY)]
    lo, hi = min(known), max(known)
    start, end = (225, 236, 238), (18, 52, 59)

    def mix(y):
        t = 0 if hi == lo else (y - lo) / (hi - lo)
        return "#%02x%02x%02x" % tuple(int(s + (e - s) * t) for s, e in zip(start, end))

    colors = [mix(y) if pd.notna(y) else GREY for y in years]
    return colors, [(str(int(lo)), mix(lo)), (str(int(hi)), mix(hi)), ("Non renseignée", GREY)]


def yes_no(value):
    if pd.isna(value) or value == "":
        return "Non renseigné"
    return "Oui" if int(value) else "Non"

# --------------------------------------------------------------------------- interface

COLOR_OPTIONS = [{"label": "Groupe au seuil choisi", "value": "cluster"}]
if HAS_ST:
    COLOR_OPTIONS.append({"label": "Lignée (ST)", "value": "st"})
COLOR_OPTIONS += [{"label": "Pays", "value": "country"},
                  {"label": "Année de prélèvement", "value": "year"}]
if HAS_AMR:
    COLOR_OPTIONS.append({"label": "Carbapénémase", "value": "carbapenemase"})

LAYOUTS = [{"label": "Force (fcose)", "value": "fcose"},
           {"label": "Force (cose)", "value": "cose"},
           {"label": "Concentrique", "value": "concentric"},
           {"label": "Cercle", "value": "circle"},
           {"label": "Grille", "value": "grid"}]

SEARCH_OPTIONS = [{"label": f"{r.id}  {r.strain}".strip(), "value": r.id}
                  for r in NODES.itertuples()]

STYLESHEET = [
    {"selector": "node", "style": {
        "background-color": "data(color)", "width": "data(size)", "height": "data(size)",
        "border-width": 1, "border-color": "#ffffff", "label": ""}},
    {"selector": "edge", "style": {
        "width": "mapData(sim, 0, 1, 0.6, 3.5)", "line-color": "#8fa3a8",
        "opacity": 0.45, "curve-style": "haystack"}},
    {"selector": "node:selected", "style": {
        "border-width": 3, "border-color": "#12343b", "label": "data(label)",
        "font-size": 9, "text-valign": "top", "text-margin-y": -4}},
    {"selector": ".focus", "style": {
        "border-width": 3, "border-color": "#12343b", "label": "data(label)", "font-size": 9}},
    {"selector": ".neighbor", "style": {"label": "data(label)", "font-size": 7}},
    {"selector": "edge.neighbor", "style": {"line-color": "#12343b", "opacity": 0.9}},
    {"selector": ".faded", "style": {"opacity": 0.12}},
]

app = Dash(__name__, title="GenomeNet", external_stylesheets=[
    "https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600"
    "&family=Source+Serif+4:opsz,wght@8..60,600&display=swap"])
server = app.server

demo_banner = html.Div(
    "Données simulées pour la démonstration. Lancez pipeline/build_network.sh "
    "pour construire le réseau à partir de vrais génomes NCBI.",
    className="banner") if IS_DEMO else None

app.layout = html.Div(className="page", children=[
    html.Header(className="top", children=[
        html.H1("GenomeNet"),
        html.P([f"Réseau de similarité de {len(NODES)} génomes de ",
                html.I(SPECIES),
                ". Deux génomes sont reliés quand leur distance Mash est sous le seuil choisi."]),
        demo_banner,
    ]),
    html.Main(className="grid", children=[
        html.Aside(className="panel controls", children=[
            html.Label("Seuil de distance Mash", htmlFor="threshold"),
            dcc.Slider(id="threshold", min=0, max=round(MAX_DIST, 4),
                       step=round(MAX_DIST / 200, 6), value=DEFAULT_THRESHOLD,
                       marks={round(MAX_DIST * k / 4, 4): f"{MAX_DIST * k / 4:.4f}" for k in range(5)},
                       tooltip={"placement": "top"}),
            html.P("Environ 0,001 correspond à des souches très proches, "
                   "typiques d'une transmission récente.", className="hint"),
            html.Label("Couleur des nœuds", htmlFor="color-by"),
            dcc.Dropdown(id="color-by", options=COLOR_OPTIONS,
                         value="st" if HAS_ST else "cluster", clearable=False),
            html.Label("Disposition", htmlFor="layout"),
            dcc.Dropdown(id="layout", options=LAYOUTS, value="fcose", clearable=False),
            html.Label("Rechercher un génome", htmlFor="search"),
            dcc.Dropdown(id="search", options=SEARCH_OPTIONS, placeholder="Accession ou souche"),
            dcc.Checklist(id="hide-isolated", className="check",
                          options=[{"label": " Masquer les génomes isolés", "value": "hide"}],
                          value=[]),
            html.Div(id="stats", className="stats"),
            html.Div(id="legend", className="legend"),
        ]),
        html.Section(className="panel graph-wrap", children=[
            cyto.Cytoscape(id="graph", elements=[], stylesheet=STYLESHEET,
                           layout={"name": "fcose"}, minZoom=0.2, maxZoom=4,
                           style={"width": "100%", "height": "100%"}),
        ]),
        html.Aside(id="details", className="panel details",
                   children=html.P("Cliquez sur un génome pour voir sa fiche.", className="hint")),
    ]),
    html.Footer(className="foot", children=[
        "Données : NCBI RefSeq. Distances : Mash. Typage : mlst, AMRFinderPlus. ",
        html.A("Code source", href="https://github.com/mamoudou-cheikh-kane/genomenet", target="_blank"),
    ]),
])

# --------------------------------------------------------------------------- callbacks

@app.callback(
    Output("graph", "elements"), Output("graph", "layout"),
    Output("legend", "children"), Output("stats", "children"),
    Input("threshold", "value"), Input("color-by", "value"), Input("layout", "value"),
    Input("hide-isolated", "value"), Input("search", "value"),
)
def update_graph(threshold, color_by, layout_name, hide_isolated, focus):
    threshold = threshold or 0
    kept = EDGES[EDGES["distance"] <= threshold]
    cluster_of, groups = connected_components(NODE_IDS, kept)

    degree = pd.concat([kept["source"], kept["target"]]).value_counts()
    nodes = NODES.copy()
    nodes["degree"] = nodes["id"].map(degree).fillna(0).astype(int)
    nodes["cluster"] = nodes["id"].map(cluster_of)
    if "hide" in (hide_isolated or []):
        nodes = nodes[nodes["degree"] > 0]

    if color_by == "year":
        colors, legend = year_colors(nodes["year"].tolist())
    elif color_by == "carbapenemase":
        values = [yes_no(v) for v in nodes["carbapenemase"]]
        mapping = {"Oui": "#b5446e", "Non": "#1f6f78", "Non renseigné": GREY}
        colors = [mapping[v] for v in values]
        legend = [(k, c) for k, c in mapping.items() if k in values]
    else:
        colors, legend = categorical_colors(nodes[color_by].replace("", "Inconnu").tolist())

    neighbors = set()
    if focus:
        mask = (kept["source"] == focus) | (kept["target"] == focus)
        neighbors = set(kept.loc[mask, "source"]) | set(kept.loc[mask, "target"])

    elements = []
    shown = set(nodes["id"])
    for row, color in zip(nodes.itertuples(), colors):
        n_amr = 0 if pd.isna(row.n_amr) else float(row.n_amr)
        classes = ""
        if focus:
            classes = "focus" if row.id == focus else ("neighbor" if row.id in neighbors else "faded")
        elements.append({"data": {
            "id": row.id, "label": row.st or row.id, "color": color,
            "size": round(14 + min(n_amr, 25) * 0.8, 1)}, "classes": classes})
    for row in kept.itertuples():
        if row.source in shown and row.target in shown:
            classes = ""
            if focus:
                classes = "neighbor" if focus in (row.source, row.target) else "faded"
            elements.append({"data": {
                "id": f"{row.source}__{row.target}", "source": row.source, "target": row.target,
                "distance": row.distance,
                "sim": 1 - row.distance / threshold if threshold else 1}, "classes": classes})

    layout = {"name": layout_name, "animate": False, "fit": True, "padding": 20}
    if layout_name == "fcose":
        layout.update({"quality": "default", "randomize": True, "nodeRepulsion": 6000,
                       "idealEdgeLength": 50, "nodeSeparation": 60})

    multi = [g for g in groups if len(g) > 1]
    stats = [
        html.Div([html.Strong(len(nodes)), " génomes affichés"]),
        html.Div([html.Strong(int(len(elements) - len(nodes))), " liens"]),
        html.Div([html.Strong(len(multi)), " groupes de plus d'un génome"]),
        html.Div([html.Strong(len(multi[0]) if multi else 0), " génomes dans le plus grand groupe"]),
    ]
    legend_items = [html.Div(className="legend-item", children=[
        html.Span(className="swatch", style={"background": c}), html.Span(label)])
        for label, c in legend]
    return elements, layout, legend_items, stats


@app.callback(Output("details", "children"), Input("graph", "tapNodeData"))
def show_details(data):
    if not data:
        return html.P("Cliquez sur un génome pour voir sa fiche.", className="hint")
    row = NODES.set_index("id").loc[data["id"]]

    def fmt(value, suffix=""):
        if value is None or (isinstance(value, float) and pd.isna(value)) or value == "":
            return "Non renseigné"
        if isinstance(value, float) and value.is_integer():
            value = int(value)
        return f"{value}{suffix}"

    fields = [
        ("Souche", fmt(row["strain"])), ("Lignée", fmt(row["st"])),
        ("Pays", fmt(row["country"])), ("Année", fmt(row["year"])),
        ("Hôte", fmt(row["host"])), ("Prélèvement", fmt(row["isolation_source"])),
        ("Taille", fmt(row["length_mb"], " Mb")), ("GC", fmt(row["gc_percent"], " %")),
        ("Carbapénémase", yes_no(row["carbapenemase"])), ("BLSE (CTX-M)", yes_no(row["esbl"])),
    ]
    genes = [g for g in row["amr_genes"].split(";") if g]
    children = [
        html.H2(data["id"]),
        html.Dl([item for k, v in fields for item in (html.Dt(k), html.Dd(v))]),
    ]
    if genes:
        children += [html.H3(f"Gènes de résistance ({len(genes)})"),
                     html.Div([html.Span(g, className="gene") for g in genes], className="genes")]
    if not IS_DEMO:
        children.append(html.A("Voir la fiche NCBI", target="_blank",
                               href=f"https://www.ncbi.nlm.nih.gov/datasets/genome/{data['id']}/"))
    return children


if __name__ == "__main__":
    app.run(debug=True)
