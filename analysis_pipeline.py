"""End-to-end bioactivity analysis for PHD finger protein 23 (ChEMBL CHEMBL2424508).

The pipeline downloads and cleans IC50 records, computes pIC50 and Lipinski descriptors,
creates exploratory plots, writes a data-driven conclusions report, and only fits a
regression model when the dataset is sufficiently informative.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from chembl_webresource_client.new_client import new_client

TARGET_CHEMBL_ID = "CHEMBL2424508"
TARGET_NAME = "PHD finger protein 23"
OUTDIR = Path("outputs")
PLOTDIR = OUTDIR / "plots"
OUTDIR.mkdir(exist_ok=True)
PLOTDIR.mkdir(exist_ok=True)


def fetch_bioactivity() -> pd.DataFrame:
    records = new_client.activity.filter(
        target_chembl_id=TARGET_CHEMBL_ID, standard_type="IC50"
    ).only(
        "molecule_chembl_id", "canonical_smiles", "standard_value",
        "standard_units", "standard_relation", "pchembl_value",
    )
    return pd.DataFrame.from_dict(records)


def classify_activity(value_nm: float) -> str:
    if value_nm <= 1000:
        return "active"
    if value_nm >= 10000:
        return "inactive"
    return "intermediate"


def clean_data(df: pd.DataFrame) -> pd.DataFrame:
    keep = ["molecule_chembl_id", "canonical_smiles", "standard_value",
            "standard_units", "standard_relation", "pchembl_value"]
    df = df[[c for c in keep if c in df.columns]].copy()
    df = df.dropna(subset=["molecule_chembl_id", "canonical_smiles", "standard_value"])
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df.dropna(subset=["standard_value"])
    df = df[df["standard_value"] > 0]
    if "standard_units" in df.columns:
        df = df[(df["standard_units"].isna()) | (df["standard_units"] == "nM")]
    df = df.drop_duplicates(subset=["molecule_chembl_id", "canonical_smiles", "standard_value"])
    df["bioactivity_class"] = df["standard_value"].map(classify_activity)
    capped = df["standard_value"].clip(upper=100_000_000)
    df["pIC50"] = 9 - np.log10(capped)
    return df.reset_index(drop=True)


def add_lipinski(df: pd.DataFrame) -> pd.DataFrame:
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski
    except ImportError:
        print("RDKit is not installed; skipping Lipinski descriptors.")
        return df

    rows = []
    for smi in df["canonical_smiles"]:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rows.append((np.nan, np.nan, np.nan, np.nan))
        else:
            rows.append((Descriptors.MolWt(mol), Descriptors.MolLogP(mol),
                         Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)))
    lip = pd.DataFrame(rows, columns=["MW", "LogP", "NumHDonors", "NumHAcceptors"])
    return pd.concat([df.reset_index(drop=True), lip], axis=1)


def save_summaries(df: pd.DataFrame) -> None:
    df.to_csv(OUTDIR / "phf23_bioactivity_processed.csv", index=False)
    df["bioactivity_class"].value_counts(dropna=False).rename_axis("class").reset_index(name="count").to_csv(
        OUTDIR / "activity_class_counts.csv", index=False)
    numeric = df.select_dtypes(include=[np.number])
    if not numeric.empty:
        numeric.describe().T.to_csv(OUTDIR / "numeric_summary.csv")


def make_plots(df: pd.DataFrame) -> None:
    """Create a dedicated plots section for biological and chemical interpretation."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; plots skipped.")
        return

    counts = df["bioactivity_class"].value_counts()
    if not counts.empty:
        ax = counts.plot(kind="bar", title=f"Bioactivity classes: {TARGET_NAME}")
        ax.set_xlabel("Bioactivity class")
        ax.set_ylabel("Number of compounds")
        plt.tight_layout()
        plt.savefig(PLOTDIR / "01_bioactivity_class_distribution.png", dpi=250)
        plt.close()

    if "pIC50" in df.columns and not df.empty:
        ax = df["pIC50"].plot(kind="hist", bins=min(15, max(5, len(df) // 5)),
                              title="Distribution of pIC50 values")
        ax.set_xlabel("pIC50")
        ax.set_ylabel("Frequency")
        plt.tight_layout()
        plt.savefig(PLOTDIR / "02_pIC50_distribution.png", dpi=250)
        plt.close()

    if {"MW", "pIC50"}.issubset(df.columns):
        plot_df = df[["MW", "pIC50"]].dropna()
        if not plot_df.empty:
            ax = plot_df.plot.scatter(x="MW", y="pIC50", title="Molecular weight vs pIC50")
            ax.set_xlabel("Molecular weight (Da)")
            ax.set_ylabel("pIC50")
            plt.tight_layout()
            plt.savefig(PLOTDIR / "03_MW_vs_pIC50.png", dpi=250)
            plt.close()

    if {"LogP", "pIC50"}.issubset(df.columns):
        plot_df = df[["LogP", "pIC50"]].dropna()
        if not plot_df.empty:
            ax = plot_df.plot.scatter(x="LogP", y="pIC50", title="LogP vs pIC50")
            ax.set_xlabel("LogP")
            ax.set_ylabel("pIC50")
            plt.tight_layout()
            plt.savefig(PLOTDIR / "04_LogP_vs_pIC50.png", dpi=250)
            plt.close()

    descriptors = [c for c in ["MW", "LogP", "NumHDonors", "NumHAcceptors", "pIC50"] if c in df.columns]
    if len(descriptors) >= 2:
        corr = df[descriptors].corr(numeric_only=True)
        fig, ax = plt.subplots(figsize=(7, 6))
        image = ax.imshow(corr.values, aspect="auto")
        ax.set_xticks(range(len(corr.columns)), corr.columns, rotation=45, ha="right")
        ax.set_yticks(range(len(corr.index)), corr.index)
        ax.set_title("Correlation matrix of molecular descriptors")
        for i in range(len(corr.index)):
            for j in range(len(corr.columns)):
                ax.text(j, i, f"{corr.iloc[i, j]:.2f}", ha="center", va="center")
        fig.colorbar(image, ax=ax, label="Pearson correlation")
        plt.tight_layout()
        plt.savefig(PLOTDIR / "05_descriptor_correlation_matrix.png", dpi=250)
        plt.close()


def write_conclusions(df: pd.DataFrame) -> None:
    """Write conclusions from the actual processed data instead of hard-coding claims."""
    counts = df["bioactivity_class"].value_counts()
    total = len(df)
    active = int(counts.get("active", 0))
    intermediate = int(counts.get("intermediate", 0))
    inactive = int(counts.get("inactive", 0))

    lines = [
        "# Conclusions",
        "",
        f"This analysis examined **{total} cleaned IC50 bioactivity records** for "
        f"**{TARGET_NAME} ({TARGET_CHEMBL_ID})** retrieved from ChEMBL.",
        "",
        "## Bioactivity profile",
        "",
        f"- Active compounds (IC50 <= 1,000 nM): **{active}**",
        f"- Intermediate compounds (1,000 < IC50 < 10,000 nM): **{intermediate}**",
        f"- Inactive compounds (IC50 >= 10,000 nM): **{inactive}**",
        "",
    ]

    if total:
        dominant = counts.idxmax()
        dominant_n = int(counts.max())
        dominant_pct = 100 * dominant_n / total
        lines += [
            "## Interpretation",
            "",
            f"The most common bioactivity class is **{dominant}**, representing "
            f"**{dominant_n}/{total} ({dominant_pct:.1f}%)** of the processed observations.",
            "",
        ]
        if "pIC50" in df.columns:
            lines.append(f"The observed pIC50 values range from **{df['pIC50'].min():.2f}** to "
                         f"**{df['pIC50'].max():.2f}**, with a median of **{df['pIC50'].median():.2f}**.")
            lines.append("")

    if df["bioactivity_class"].nunique() < 2:
        lines += [
            "## Main limitation",
            "",
            "Only one bioactivity class is represented after preprocessing. Therefore, an "
            "active-vs-inactive classification model would not be scientifically meaningful for "
            "this data slice. This is a limitation of the available target data rather than a "
            "reason to manufacture or relabel observations.",
            "",
        ]
    elif min(active, inactive) < 5:
        lines += [
            "## Main limitation",
            "",
            "The active and inactive classes are strongly imbalanced. Any classification result "
            "should therefore be interpreted cautiously and evaluated with methods appropriate "
            "for imbalanced data.",
            "",
        ]

    lines += [
        "## Overall conclusion",
        "",
        "The workflow successfully converts raw ChEMBL measurements into a reproducible dataset "
        "containing activity classes, pIC50 values and molecular descriptors, together with plots "
        "that summarize the chemical and bioactivity space. The amount and diversity of available "
        "PHF23 bioactivity data should determine whether predictive modeling is justified. Future "
        "work could expand the dataset with additional compatible activity measurements or apply "
        "the same workflow to a target with a larger and better-balanced ChEMBL dataset.",
        "",
        "## Generated plots",
        "",
        "Plots are saved in `outputs/plots/` and include the bioactivity-class distribution, "
        "pIC50 distribution, molecular-weight/pIC50 relationship, LogP/pIC50 relationship and "
        "descriptor correlation matrix.",
    ]
    (OUTDIR / "CONCLUSIONS.md").write_text("\n".join(lines), encoding="utf-8")


def optional_regression(df: pd.DataFrame) -> None:
    features = ["MW", "LogP", "NumHDonors", "NumHAcceptors"]
    if not all(c in df.columns for c in features):
        print("Regression skipped: Lipinski descriptor columns are unavailable.")
        return
    model_df = df[features + ["pIC50"]].dropna()
    if len(model_df) < 20 or model_df["pIC50"].nunique() < 5:
        print(f"Regression skipped: insufficient informative observations (n={len(model_df)}, "
              f"unique pIC50={model_df['pIC50'].nunique()}).")
        return
    try:
        from sklearn.ensemble import RandomForestRegressor
        from sklearn.metrics import mean_squared_error, r2_score
        from sklearn.model_selection import train_test_split
    except ImportError:
        print("scikit-learn is not installed; regression skipped.")
        return
    X, y = model_df[features], model_df["pIC50"]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    model = RandomForestRegressor(n_estimators=300, random_state=42)
    model.fit(X_train, y_train)
    pred = model.predict(X_test)
    pd.DataFrame({"metric": ["R2", "RMSE", "n_train", "n_test"],
                  "value": [r2_score(y_test, pred), math.sqrt(mean_squared_error(y_test, pred)),
                            len(X_train), len(X_test)]}).to_csv(OUTDIR / "random_forest_metrics.csv", index=False)
    pd.DataFrame({"observed_pIC50": y_test, "predicted_pIC50": pred}).to_csv(
        OUTDIR / "random_forest_predictions.csv", index=False)


def main() -> None:
    print(f"Fetching ChEMBL IC50 data for {TARGET_NAME} ({TARGET_CHEMBL_ID})...")
    raw = fetch_bioactivity()
    raw.to_csv(OUTDIR / "phf23_bioactivity_raw.csv", index=False)
    print(f"Retrieved {len(raw)} raw records.")
    processed = add_lipinski(clean_data(raw))
    save_summaries(processed)
    make_plots(processed)
    write_conclusions(processed)
    print("\nActivity-class counts:")
    print(processed["bioactivity_class"].value_counts(dropna=False))
    optional_regression(processed)
    print(f"\nFinished. Results are in: {OUTDIR.resolve()}")
    print("Plots: outputs/plots/")
    print("Conclusions: outputs/CONCLUSIONS.md")


if __name__ == "__main__":
    main()
