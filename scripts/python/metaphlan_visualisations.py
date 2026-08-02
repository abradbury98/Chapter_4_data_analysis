"""
MetaPhlAn4 visualisations for barcode01_profile.txt
Covers:
  1. Phylum-level bar chart
  2. Species-level bar chart
  3. Phylum pie chart
  4. GraPhlAn pipeline (shell commands printed at end)

Requirements:
  pip install pandas matplotlib seaborn
  conda install -c biobakery graphlan export2graphlan   # for GraPhlAn
"""

import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

# ── 1. Parse MetaPhlAn profile ────────────────────────────────────────────────

PROFILE = Path("barcode01_profile.txt")

def parse_metaphlan(path: Path) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            clade, _, abundance, *_ = parts
            rows.append({"clade": clade, "abundance": float(abundance)})
    return pd.DataFrame(rows)

df = parse_metaphlan(PROFILE)

def level_df(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """Extract rows at a specific taxonomic level (k, p, c, o, f, g, s, t)."""
    pat = re.compile(rf"^(k__[^|]+\|)*{prefix}__[^|]+$")
    sub = df[df["clade"].apply(lambda c: bool(pat.match(c)))].copy()
    sub["name"] = sub["clade"].apply(lambda c: c.split("|")[-1].replace(f"{prefix}__", "").replace("_", " "))
    return sub.sort_values("abundance", ascending=False).reset_index(drop=True)

phylum_df  = level_df(df, "p")
species_df = level_df(df, "s")

# ── 2. Colour palette ─────────────────────────────────────────────────────────

PHYLUM_COLOURS = {
    "Firmicutes":       "#4C72B0",
    "Campylobacterota": "#DD8452",
    "Pseudomonadota":   "#55A868",
    "Bacteroidota":     "#C44E52",
}

def phylum_colour(species_clade: str) -> str:
    for phylum, colour in PHYLUM_COLOURS.items():
        if phylum in species_clade:
            return colour
    return "#8C8C8C"

# ── 3. Phylum bar chart ───────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(7, 4))
colours = [PHYLUM_COLOURS.get(n, "#8C8C8C") for n in phylum_df["name"]]
bars = ax.barh(phylum_df["name"], phylum_df["abundance"], color=colours, edgecolor="white")
ax.bar_label(bars, fmt="%.2f%%", padding=3, fontsize=9)
ax.set_xlabel("Relative abundance (%)", fontsize=11)
ax.set_title("Barcode01 – Phylum-level relative abundance", fontsize=12, pad=10)
ax.invert_yaxis()
ax.set_xlim(0, phylum_df["abundance"].max() * 1.15)
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig("barcode01_phylum_bar.png", dpi=200)
plt.close()
print("Saved: barcode01_phylum_bar.png")

# ── 4. Species bar chart ──────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(8, 6))
sp_colours = [phylum_colour(c) for c in species_df["clade"]]
bars = ax.barh(species_df["name"], species_df["abundance"], color=sp_colours, edgecolor="white")
ax.bar_label(bars, fmt="%.2f%%", padding=3, fontsize=8)
ax.set_xlabel("Relative abundance (%)", fontsize=11)
ax.set_title("Barcode01 – Species-level relative abundance", fontsize=12, pad=10)
ax.invert_yaxis()
ax.set_xlim(0, species_df["abundance"].max() * 1.18)
ax.spines[["top", "right"]].set_visible(False)
legend_patches = [mpatches.Patch(color=c, label=p) for p, c in PHYLUM_COLOURS.items()]
ax.legend(handles=legend_patches, title="Phylum", fontsize=8, title_fontsize=9,
          loc="lower right")
plt.tight_layout()
plt.savefig("barcode01_species_bar.png", dpi=200)
plt.close()
print("Saved: barcode01_species_bar.png")

# ── 5. Phylum pie chart ───────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(6, 6))
wedge_colours = [PHYLUM_COLOURS.get(n, "#8C8C8C") for n in phylum_df["name"]]
wedges, texts, autotexts = ax.pie(
    phylum_df["abundance"],
    labels=phylum_df["name"],
    colors=wedge_colours,
    autopct="%1.1f%%",
    startangle=140,
    pctdistance=0.75,
    wedgeprops=dict(edgecolor="white", linewidth=1.5),
)
for t in autotexts:
    t.set_fontsize(9)
ax.set_title("Barcode01 – Phylum composition", fontsize=12, pad=12)
plt.tight_layout()
plt.savefig("barcode01_phylum_pie.png", dpi=200)
plt.close()
print("Saved: barcode01_phylum_pie.png")

# ── 6. GraPhlAn pipeline ─────────────────────────────────────────────────────

GRAPHLAN_CMDS = """
# ── GraPhlAn pipeline (run in bash) ──────────────────────────────────────────
# Install once:
#   conda install -c biobakery graphlan export2graphlan

# Step 1 – convert MetaPhlAn profile to tree + annotation files
export2graphlan.py \\
    -i barcode01_profile.txt \\
    -t barcode01.tree \\
    -a barcode01.annot \\
    --title "Barcode01" \\
    --annotations 2,3,4,5,6,7 \\
    --external_annotations 7 \\
    --min_clade_size 1 \\
    --def_clade_size 10 \\
    --def_font_size 10 \\
    --most_abundant 25 \\
    --abundance_threshold 1

# Step 2 – annotate the tree
graphlan_annotate.py \\
    --annot barcode01.annot \\
    barcode01.tree \\
    barcode01_annotated.xml

# Step 3 – render to PNG
graphlan.py \\
    --dpi 150 \\
    --size 14 \\
    barcode01_annotated.xml \\
    barcode01_cladogram.png
"""

print(GRAPHLAN_CMDS)
