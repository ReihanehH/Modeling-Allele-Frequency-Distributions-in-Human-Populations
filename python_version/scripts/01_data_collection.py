import pandas as pd
import os

# ── Paths ────────────────────────────────────────────────────────────────────
RAW_DIR  = "data/raw"
OUT_DIR  = "data/processed"
os.makedirs(RAW_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

# ── Gene files ───────────────────────────────────────────────────────────────
FILES = {
    "APOE":   "data/raw/gnomad_APOE.csv",
    "CYP2C19": "data/raw/gnomad_CYP2C19.csv",
    "HLA-B":  "data/raw/gnomad_HLA-B.csv",
}

# ── Populations present in gnomAD export ─────────────────────────────────────
POPULATIONS = [
    "African/African American",
    "Admixed American",
    "Ashkenazi Jewish",
    "East Asian",
    "European (Finnish)",
    "Middle Eastern",
    "European (non-Finnish)",
    "Amish",
    "South Asian",
    "Remaining",
]

# Columns we keep from the original file for each variant
VARIANT_COLS = [
    "gnomAD ID", "Chromosome", "Position", "Reference", "Alternate",
    "VEP Annotation", "Allele Frequency",          # overall AF
    "Allele Count", "Allele Number",                # overall AC / AN
]

def tidy_gene(gene, filepath):
    """Read one gnomAD CSV and return a tidy long-format DataFrame."""
    df = pd.read_csv(filepath, low_memory=False)
    df["Gene"] = gene

    rows = []

    for _, row in df.iterrows():
        base = {col: row.get(col) for col in VARIANT_COLS if col in df.columns}
        base["Gene"] = gene

        # Overall population row
        rows.append({
            **base,
            "Population": "Overall",
            "AC": row.get("Allele Count"),
            "AN": row.get("Allele Number"),
            "AF": row.get("Allele Frequency"),
            "Homozygote Count": row.get("Homozygote Count"),
        })

        # Per-population rows
        for pop in POPULATIONS:
            ac = row.get(f"Allele Count {pop}")
            an = row.get(f"Allele Number {pop}")
            hom = row.get(f"Homozygote Count {pop}")

            # Compute AF from AC/AN (avoid division by zero)
            if pd.notna(an) and an > 0 and pd.notna(ac):
                af = ac / an
            else:
                af = None

            rows.append({
                **base,
                "Population": pop,
                "AC": ac,
                "AN": an,
                "AF": af,
                "Homozygote Count": hom,
            })

    tidy = pd.DataFrame(rows)

    # Drop rows where both AC and AN are missing
    tidy = tidy.dropna(subset=["AC", "AN"])

    return tidy


def main():
    all_dfs = []

    for gene, filepath in FILES.items():
        print(f"Processing {gene}...")
        tidy = tidy_gene(gene, filepath)

        out_path = f"{OUT_DIR}/{gene}_tidy.csv"
        tidy.to_csv(out_path, index=False)
        print(f"  {len(tidy)} rows → {out_path}")

        all_dfs.append(tidy)

    # Combined file
    combined = pd.concat(all_dfs, ignore_index=True)
    combined.to_csv(f"{OUT_DIR}/all_genes_tidy.csv", index=False)
    print(f"\nCombined: {len(combined)} rows → {OUT_DIR}/all_genes_tidy.csv")

    # Summary
    print("\n── Variant counts by gene and population ──")
    summary = (
        combined[combined["Population"] != "Overall"]
        .groupby(["Gene", "Population"])
        .size()
        .unstack(fill_value=0)
    )
    print(summary.to_string())


if __name__ == "__main__":
    main()