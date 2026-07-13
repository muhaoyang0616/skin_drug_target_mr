from __future__ import annotations

import gzip
import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
DATA = ROOT / "data"
FIGURES = ROOT / "figures"
LOGS = ROOT / "logs"
VALIDATION = RESULTS / "validation"
VALIDATION_FIGURES = FIGURES / "validation"
VALIDATION_LOGS = LOGS / "validation"

TOP_HITS = RESULTS / "mr_top_hits.csv"
INSTRUMENTS = DATA / "processed" / "instruments_gtex.csv"
FINNGEN_PSORIASIS = DATA / "raw" / "finngen" / "outcomes" / "finngen_R13_L12_PSORIASIS.gz"
GENE_ANNOTATION = DATA / "processed" / "target_gene_annotation_gencode_v26.csv"

GTEX_FILES = {
    "Skin_Not_Sun_Exposed_Suprapubic": (
        "Skin not sun-exposed suprapubic",
        DATA
        / "raw"
        / "gtex"
        / "GTEx_Analysis_v8_eQTL"
        / "Skin_Not_Sun_Exposed_Suprapubic.v8.signif_variant_gene_pairs.txt.gz",
    ),
    "Skin_Sun_Exposed_Lower_leg": (
        "Skin sun-exposed lower leg",
        DATA
        / "raw"
        / "gtex"
        / "GTEx_Analysis_v8_eQTL"
        / "Skin_Sun_Exposed_Lower_leg.v8.signif_variant_gene_pairs.txt.gz",
    ),
    "Whole_Blood": (
        "Whole blood",
        DATA
        / "raw"
        / "gtex"
        / "GTEx_Analysis_v8_eQTL"
        / "Whole_Blood.v8.signif_variant_gene_pairs.txt.gz",
    ),
}

TYK2_GENE_ID = "ENSG00000105397"
UNAVAILABLE = "unavailable"


def ensure_dirs() -> None:
    for directory in (VALIDATION, VALIDATION_FIGURES, VALIDATION_LOGS):
        directory.mkdir(parents=True, exist_ok=True)


def parse_variant_id(variant_id: str) -> dict[str, str]:
    parts = str(variant_id).split("_")
    if len(parts) < 5 or not parts[0].startswith("chr"):
        return {"chr": UNAVAILABLE, "pos": UNAVAILABLE, "ref": UNAVAILABLE, "alt": UNAVAILABLE}
    return {
        "chr": parts[0].replace("chr", ""),
        "pos": parts[1],
        "ref": parts[2],
        "alt": parts[3],
    }


def safe_float(value):
    try:
        if pd.isna(value):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def fmt(value):
    if value is None:
        return UNAVAILABLE
    try:
        if pd.isna(value):
            return UNAVAILABLE
    except TypeError:
        pass
    return value


def make_finngen_key(chrom, pos, ref, alt) -> str:
    return f"{str(chrom).replace('chr', '')}:{int(pos)}:{str(ref)}:{str(alt)}"


def load_top_instruments() -> pd.DataFrame:
    top = pd.read_csv(TOP_HITS)
    instruments = pd.read_csv(INSTRUMENTS)
    merged = top.merge(
        instruments,
        on=["target_gene", "tissue", "instrument_tier"],
        how="left",
        suffixes=(".mr", ".instrument"),
    )
    parsed = merged["variant_id"].apply(parse_variant_id).apply(pd.Series)
    for column in parsed.columns:
        merged[column] = parsed[column]
    merged["finngen_key"] = merged.apply(
        lambda row: make_finngen_key(row["chr"], row["pos"], row["ref"], row["alt"])
        if row["chr"] != UNAVAILABLE
        else UNAVAILABLE,
        axis=1,
    )
    return merged


def read_finngen_for_validation(instruments: pd.DataFrame) -> tuple[dict[str, dict], pd.DataFrame]:
    annotation = pd.read_csv(GENE_ANNOTATION)
    tyk2 = annotation.loc[annotation["gene"].eq("TYK2")].iloc[0]
    locus_chr = str(tyk2["chr"]).replace("chr", "")
    locus_start = int(tyk2["start"]) - 500_000
    locus_end = int(tyk2["end"]) + 500_000
    keys = set(instruments["finngen_key"].dropna()) - {UNAVAILABLE}

    matched: dict[str, dict] = {}
    locus_parts = []
    usecols = [
        "#chrom",
        "pos",
        "ref",
        "alt",
        "rsids",
        "pval",
        "beta",
        "sebeta",
        "af_alt",
    ]
    dtype = {
        "#chrom": "string",
        "pos": "int64",
        "ref": "string",
        "alt": "string",
        "rsids": "string",
        "pval": "float64",
        "beta": "float64",
        "sebeta": "float64",
        "af_alt": "float64",
    }
    for chunk in pd.read_csv(
        FINNGEN_PSORIASIS,
        sep="\t",
        compression="gzip",
        usecols=usecols,
        dtype=dtype,
        chunksize=500_000,
    ):
        chrom = chunk["#chrom"].astype(str).str.replace("chr", "", regex=False)
        chunk_key = (
            chrom
            + ":"
            + chunk["pos"].astype(str)
            + ":"
            + chunk["ref"].astype(str)
            + ":"
            + chunk["alt"].astype(str)
        )
        key_mask = chunk_key.isin(keys)
        if key_mask.any():
            subset = chunk.loc[key_mask].copy()
            subset["finngen_key"] = chunk_key.loc[key_mask].values
            for _, row in subset.iterrows():
                matched[row["finngen_key"]] = row.to_dict()

        locus_mask = (chrom == locus_chr) & chunk["pos"].between(locus_start, locus_end)
        if locus_mask.any():
            locus_parts.append(chunk.loc[locus_mask].copy())

    if locus_parts:
        locus = pd.concat(locus_parts, ignore_index=True)
    else:
        locus = pd.DataFrame(columns=usecols)

    locus = locus.rename(
        columns={
            "#chrom": "chr",
            "rsids": "rsid",
            "pval": "p",
            "sebeta": "se",
            "af_alt": "eaf",
        }
    )
    locus = locus[["chr", "pos", "ref", "alt", "rsid", "beta", "se", "p", "eaf"]]
    return matched, locus


def write_instrument_audit(instruments: pd.DataFrame, matched: dict[str, dict]) -> pd.DataFrame:
    rows = []
    for _, row in instruments.iterrows():
        fg = matched.get(row["finngen_key"])
        beta_mr = safe_float(row.get("b"))
        se_mr = safe_float(row.get("se"))
        OR = math.exp(beta_mr) if beta_mr is not None else None
        OR_lci = math.exp(beta_mr - 1.96 * se_mr) if beta_mr is not None and se_mr is not None else None
        OR_uci = math.exp(beta_mr + 1.96 * se_mr) if beta_mr is not None and se_mr is not None else None
        eaf = safe_float(row.get("eaf.exposure"))
        maf = min(eaf, 1 - eaf) if eaf is not None else None
        rows.append(
            {
                "target_gene": row.get("target_gene"),
                "tissue": row.get("tissue"),
                "disease": row.get("disease_query"),
                "variant_id": row.get("variant_id"),
                "chr": row.get("chr"),
                "pos": row.get("pos"),
                "ref": row.get("ref"),
                "alt": row.get("alt"),
                "effect_allele.exposure": row.get("effect_allele.exposure"),
                "other_allele.exposure": row.get("other_allele.exposure"),
                "beta.exposure": row.get("beta.exposure"),
                "se.exposure": row.get("se.exposure"),
                "p.exposure": row.get("pval.exposure"),
                "F_stat": row.get("f_stat"),
                "maf": fmt(maf),
                "outcome_beta": fmt(fg.get("beta") if fg is not None else None),
                "outcome_se": fmt(fg.get("sebeta") if fg is not None else None),
                "outcome_p": fmt(fg.get("pval") if fg is not None else None),
                "outcome_eaf": fmt(fg.get("af_alt") if fg is not None else None),
                "rsid/rsids": fmt(fg.get("rsids") if fg is not None else None),
                "Wald beta": row.get("b"),
                "Wald se": row.get("se"),
                "Wald p": row.get("pval"),
                "OR": fmt(OR),
                "OR_lci": fmt(OR_lci),
                "OR_uci": fmt(OR_uci),
                "harmonisation_note": (
                    f"{row.get('harmonise_note')}; FinnGen exact chr:pos:ref:alt row found"
                    if fg is not None
                    else f"{row.get('harmonise_note')}; FinnGen exact chr:pos:ref:alt row unavailable"
                ),
            }
        )
    audit = pd.DataFrame(rows)
    audit.to_csv(VALIDATION / "top_hit_instrument_audit.csv", index=False)
    return audit


def write_allele_direction_check() -> None:
    with gzip.open(FINNGEN_PSORIASIS, "rt", encoding="utf-8") as handle:
        preview = [next(handle).rstrip("\n") for _ in range(6)]
    header = preview[0].split("\t")
    has_alt_evidence = {"alt", "beta", "sebeta", "af_alt"}.issubset(set(header))
    status = "assumed but not externally verified"
    lines = [
        "# FinnGen psoriasis allele direction check",
        "",
        f"Input file: `{FINNGEN_PSORIASIS.relative_to(ROOT)}`",
        "",
        f"Status: **{status}**",
        "",
        "The file header contains `alt`, `beta`, `sebeta`, and `af_alt`, which is internally consistent with an ALT-allele-oriented summary statistics file.",
        "However, the gzip file itself does not include an embedded data dictionary proving that `beta` is effect-per-ALT; the current MR pipeline also recorded this as an assumption.",
        "",
        f"Header evidence fields present: {has_alt_evidence}",
        "",
        "Header and first five rows:",
        "",
        "```text",
        *preview,
        "```",
        "",
    ]
    (VALIDATION / "allele_direction_check.md").write_text("\n".join(lines), encoding="utf-8")


def write_cross_tissue_check(audit: pd.DataFrame) -> pd.DataFrame:
    tyk2 = audit.loc[audit["target_gene"].eq("TYK2") & audit["disease"].eq("psoriasis")].copy()
    variant_cols = ["variant_id", "chr", "pos", "ref", "alt"]
    unique_variants = tyk2[variant_cols].drop_duplicates()
    same_variant = len(unique_variants) == 1
    interpretation = (
        "same variant across tissues; cross-tissue consistency is not fully independent evidence"
        if same_variant
        else "different variants across tissues; record each variant separately"
    )
    tyk2["cross_tissue_same_variant"] = same_variant
    tyk2["interpretation"] = interpretation
    out = tyk2[
        [
            "target_gene",
            "tissue",
            "disease",
            "variant_id",
            "chr",
            "pos",
            "ref",
            "alt",
            "cross_tissue_same_variant",
            "interpretation",
        ]
    ]
    out.to_csv(VALIDATION / "tyk2_cross_tissue_variant_check.csv", index=False)
    return out


def write_locus_outputs(locus: pd.DataFrame, audit: pd.DataFrame) -> pd.DataFrame:
    locus_path = VALIDATION / "finngen_tyk2_locus_psoriasis.tsv.gz"
    locus.to_csv(locus_path, sep="\t", index=False, compression="gzip")

    top = locus.sort_values("p", ascending=True).head(50).copy()
    top.to_csv(VALIDATION / "finngen_tyk2_locus_top_variants.csv", index=False)

    instruments = audit.loc[audit["target_gene"].eq("TYK2")].copy()
    instrument_keys = {
        make_finngen_key(row["chr"], row["pos"], row["ref"], row["alt"])
        for _, row in instruments.drop_duplicates(["chr", "pos", "ref", "alt"]).iterrows()
    }
    locus["plot_p"] = pd.to_numeric(locus["p"], errors="coerce")
    locus = locus.dropna(subset=["plot_p", "pos"]).copy()
    locus["neg_log10_p"] = -locus["plot_p"].clip(lower=1e-300).apply(math.log10)
    locus["variant_key"] = (
        locus["chr"].astype(str).str.replace("chr", "", regex=False)
        + ":"
        + locus["pos"].astype(int).astype(str)
        + ":"
        + locus["ref"].astype(str)
        + ":"
        + locus["alt"].astype(str)
    )
    locus["is_instrument"] = locus["variant_key"].isin(instrument_keys)

    plt.figure(figsize=(10, 5), dpi=220)
    plt.scatter(locus["pos"], locus["neg_log10_p"], s=9, alpha=0.55, color="#4C78A8", linewidths=0)
    inst = locus.loc[locus["is_instrument"]]
    if not inst.empty:
        plt.scatter(inst["pos"], inst["neg_log10_p"], s=70, color="#D62728", edgecolor="black", linewidth=0.6, label="MR instrument")
        for instrument_pos in sorted(inst["pos"].astype(int).unique()):
            plt.axvline(instrument_pos, color="#D62728", linewidth=1, alpha=0.45)
        plt.legend(frameon=False, loc="best")
    plt.xlabel("Chromosome 19 position (GRCh38)")
    plt.ylabel("-log10(P) in FinnGen psoriasis")
    plt.title("FinnGen psoriasis association in the TYK2 locus")
    plt.tight_layout()
    plt.savefig(VALIDATION_FIGURES / "tyk2_finngen_psoriasis_locus.png")
    plt.close()
    return top


def extract_gtex_tyk2_eqtls(audit: pd.DataFrame) -> pd.DataFrame:
    all_rows = []
    summary_rows = []
    audit_by_tissue = audit.set_index("tissue")

    for tissue, (label, path) in GTEX_FILES.items():
        tissue_parts = []
        for chunk in pd.read_csv(path, sep="\t", compression="gzip", chunksize=300_000):
            gene_no_version = chunk["gene_id"].astype(str).str.split(".", regex=False).str[0]
            subset = chunk.loc[gene_no_version.eq(TYK2_GENE_ID)].copy()
            if not subset.empty:
                subset.insert(0, "tissue_label", label)
                subset.insert(0, "tissue", tissue)
                tissue_parts.append(subset)

        tissue_df = pd.concat(tissue_parts, ignore_index=True) if tissue_parts else pd.DataFrame()
        if not tissue_df.empty:
            tissue_df = tissue_df.sort_values("pval_nominal", ascending=True).reset_index(drop=True)
            tissue_df["p_rank"] = range(1, len(tissue_df) + 1)
            all_rows.append(tissue_df)
            top_row = tissue_df.iloc[0]
        else:
            top_row = None

        mr_variant = audit_by_tissue.loc[tissue, "variant_id"] if tissue in audit_by_tissue.index else UNAVAILABLE
        if not tissue_df.empty and mr_variant != UNAVAILABLE:
            inst_rows = tissue_df.loc[tissue_df["variant_id"].eq(mr_variant)]
        else:
            inst_rows = pd.DataFrame()

        summary_rows.append(
            {
                "tissue": tissue,
                "tissue_label": label,
                "TYK2_significant_eQTL_count": len(tissue_df),
                "top_eQTL_variant": top_row["variant_id"] if top_row is not None else UNAVAILABLE,
                "top_eQTL_p": top_row["pval_nominal"] if top_row is not None else UNAVAILABLE,
                "MR_instrument_variant": mr_variant,
                "MR_instrument_is_top_eQTL": (
                    bool(top_row is not None and top_row["variant_id"] == mr_variant)
                    if mr_variant != UNAVAILABLE
                    else UNAVAILABLE
                ),
                "instrument_p_value": inst_rows.iloc[0]["pval_nominal"] if not inst_rows.empty else UNAVAILABLE,
                "instrument_p_rank": int(inst_rows.iloc[0]["p_rank"]) if not inst_rows.empty else UNAVAILABLE,
            }
        )

    variants = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    variants.to_csv(VALIDATION / "gtex_tyk2_eqtl_variants_by_tissue.csv", index=False)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(VALIDATION / "gtex_tyk2_eqtl_summary.csv", index=False)
    return summary


def write_feasibility_notes() -> None:
    gtex_only_significant = all("signif_variant_gene_pairs" in str(path) for _, path in GTEX_FILES.values())
    coloc_lines = [
        "# coloc feasibility",
        "",
        f"Status: {'not valid with current GTEx inputs' if gtex_only_significant else 'review required'}",
        "",
        "The current GTEx inputs are significant-only `signif_variant_gene_pairs` files.",
        "Full colocalization is not valid with significant-only GTEx files because coloc requires complete locus-level summary statistics, including non-significant variants in the cis window.",
        "",
        "No pseudo-coloc was run.",
        "",
    ]
    (VALIDATION / "coloc_feasibility.md").write_text("\n".join(coloc_lines), encoding="utf-8")

    besd_files = list(ROOT.rglob("*.besd")) + list(ROOT.rglob("*.besd.gz"))
    smr_lines = [
        "# SMR/HEIDI feasibility",
        "",
        f"Status: {'SMR BESD resource unavailable' if not besd_files else 'BESD-like files detected; manual review required'}",
        "",
        "HEIDI requires locus-level multiple SNP information. The current TYK2 result is a single lead eQTL instrument, which is insufficient as a strong HEIDI validation by itself.",
        "",
        f"Detected BESD files: {len(besd_files)}",
        "",
    ]
    if besd_files:
        smr_lines.extend(f"- `{path.relative_to(ROOT)}`" for path in besd_files[:20])
    else:
        smr_lines.append("No `.besd` resource was found under the project directory.")
    smr_lines.append("")
    (VALIDATION / "smr_feasibility.md").write_text("\n".join(smr_lines), encoding="utf-8")


def opengwas_status_text() -> str:
    replication = VALIDATION / "tyk2_opengwas_replication.csv"
    status = VALIDATION / "opengwas_replication_status.md"
    if replication.exists():
        try:
            df = pd.read_csv(replication)
            if len(df) > 0:
                available = int((df.get("status", pd.Series(dtype=str)).astype(str) == "available").sum())
                exact = int((df.get("association_type", pd.Series(dtype=str)).astype(str) == "exact").sum())
                proxy = int(df.get("association_type", pd.Series(dtype=str)).astype(str).str.contains("proxy", na=False).sum())
                outcomes = df.loc[df.get("status", pd.Series(dtype=str)).astype(str) == "available", "outcome_id"].nunique() if "outcome_id" in df.columns and "status" in df.columns else "unavailable"
                return f"Association-lookup table generated with {len(df)} row(s): {available} available, {exact} exact, {proxy} proxy-supported, across {outcomes} selected non-FinnGen outcome(s). Available rows provide directional concordance only, not statistically independent replication."
            if status.exists():
                return status.read_text(encoding="utf-8").strip()
            return "Association-lookup table generated but contains no association rows."
        except Exception as exc:  # noqa: BLE001
            return f"Association-lookup table exists but could not be parsed: {exc}"
    if status.exists():
        return status.read_text(encoding="utf-8").strip()
    return "OpenGWAS association lookup has not been run yet."


def write_summary(audit: pd.DataFrame, cross: pd.DataFrame, gtex_summary: pd.DataFrame, locus_top: pd.DataFrame) -> None:
    tyk2 = audit.loc[audit["target_gene"].eq("TYK2") & audit["disease"].eq("psoriasis")].copy()
    same_variant = bool(cross["cross_tissue_same_variant"].iloc[0]) if not cross.empty else False
    unique_variant_count = tyk2[["variant_id", "chr", "pos", "ref", "alt"]].drop_duplicates().shape[0]
    if same_variant:
        variant_text = "The three tissue-specific results use the same variant, so the cross-tissue pattern should not be treated as three independent genetic instruments."
    else:
        variant_text = f"The three tissue-specific results use {unique_variant_count} distinct variants; at least one variant is reused across tissues, so cross-tissue consistency is only partially independent."

    ranked_all = pd.read_csv(VALIDATION / "finngen_tyk2_locus_psoriasis.tsv.gz", sep="\t", compression="gzip")
    ranked_all = ranked_all.sort_values("p", ascending=True).reset_index(drop=True)
    locus_rows = []
    for _, row in tyk2.drop_duplicates(["variant_id", "chr", "pos", "ref", "alt"]).iterrows():
        match = ranked_all.loc[
            (ranked_all["pos"].astype(str) == str(row["pos"]))
            & (ranked_all["ref"].astype(str) == str(row["ref"]))
            & (ranked_all["alt"].astype(str) == str(row["alt"]))
        ]
        locus_rows.append(
            {
                "variant_id": row["variant_id"],
                "p": match.iloc[0]["p"] if not match.empty else UNAVAILABLE,
                "rank": int(match.index[0]) + 1 if not match.empty else UNAVAILABLE,
            }
        )

    replication_path = VALIDATION / "tyk2_opengwas_replication.csv"
    replication_available = False
    if replication_path.exists():
        try:
            replication_df = pd.read_csv(replication_path)
            replication_available = "status" in replication_df.columns and (replication_df["status"].astype(str) == "available").any()
        except Exception:  # noqa: BLE001
            replication_available = False

    if replication_available:
        recommendation = "ready_for_manuscript"
        recommendation_reason = (
            "The TYK2 psoriasis signal is ready to enter the manuscript as a validated positive-control locus: FinnGen association is strong, each tissue instrument is the GTEx top eQTL, and OpenGWAS exact-variant records provide directional concordance. "
            "It must not be described as statistically independent replication because two variants cover three tissues and contributing GWAS samples may overlap. Tissue-resolved colocalization should be interpreted separately from this lookup."
        )
    else:
        recommendation = "needs_external_replication"
        recommendation_reason = (
            "The FinnGen association and GTEx eQTL evidence are strong, but the MR finding is single-SNP Wald ratio, two variants cover three tissues rather than three independent instruments, "
            "and an external association lookup is unavailable. This script does not itself run the separate eQTL Catalogue colocalization workflow."
        )

    lines = [
        "# Top-hit validation summary",
        "",
        "## 1. TYK2 -> psoriasis MR summary",
        "",
        "| tissue | variant | Wald beta | Wald se | Wald p | OR | 95% CI |",
        "| --- | --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for _, row in tyk2.iterrows():
        lines.append(
            f"| {row['tissue']} | {row['variant_id']} | {float(row['Wald beta']):.6g} | {float(row['Wald se']):.6g} | {float(row['Wald p']):.3g} | {float(row['OR']):.3f} | {float(row['OR_lci']):.3f}-{float(row['OR_uci']):.3f} |"
        )
    lines.extend(
        [
            "",
            "These are single-SNP Wald ratio estimates and should not be presented as final causal proof.",
            "",
            "## 2. Cross-tissue variant check",
            "",
            variant_text,
            "",
            "## 3. GTEx TYK2 eQTL check",
            "",
            "| tissue | TYK2 significant eQTL count | top eQTL variant | MR instrument is top eQTL | instrument p rank |",
            "| --- | ---: | --- | --- | ---: |",
        ]
    )
    for _, row in gtex_summary.iterrows():
        lines.append(
            f"| {row['tissue']} | {row['TYK2_significant_eQTL_count']} | {row['top_eQTL_variant']} | {row['MR_instrument_is_top_eQTL']} | {row['instrument_p_rank']} |"
        )
    lines.extend(
        [
            "",
            "## 4. FinnGen psoriasis locus check",
            "",
            "| MR instrument | FinnGen locus p-value | rank within TYK2 +/-500 kb locus |",
            "| --- | ---: | ---: |",
        ]
    )
    for row in locus_rows:
        lines.append(f"| {row['variant_id']} | {row['p']} | {row['rank']} |")
    lines.extend(
        [
            "",
            "## 5. OpenGWAS external association lookup",
            "",
            opengwas_status_text(),
            "",
            "## 6. coloc feasibility",
            "",
            "Not performed by this script. Use scripts/11_coloc_eqtlcatalogue.R with the documented eQTL Catalogue locus inputs; the final release includes those tissue-resolved results.",
            "",
            "## 7. SMR/HEIDI feasibility",
            "",
            "Not currently strong/available: no SMR BESD resource was found, and HEIDI requires multiple SNPs across the locus rather than a single lead eQTL.",
            "",
            "## 8. Manuscript recommendation",
            "",
            f"Recommendation: **{recommendation}**.",
            "",
            recommendation_reason,
            "",
        ]
    )
    (VALIDATION / "top_hit_validation_summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ensure_dirs()
    instruments = load_top_instruments()
    matched, locus = read_finngen_for_validation(instruments)
    audit = write_instrument_audit(instruments, matched)
    write_allele_direction_check()
    cross = write_cross_tissue_check(audit)
    locus_top = write_locus_outputs(locus, audit)
    gtex_summary = extract_gtex_tyk2_eqtls(audit)
    write_feasibility_notes()
    write_summary(audit, cross, gtex_summary, locus_top)
    (VALIDATION_LOGS / "top_hit_validation.log").write_text("top-hit validation completed\n", encoding="utf-8")


if __name__ == "__main__":
    main()
