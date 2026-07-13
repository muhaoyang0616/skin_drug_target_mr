#!/usr/bin/env python3
"""Reproduce the primary GTEx-FinnGen MR analysis from local raw files.

Run from any directory by providing explicit ``--source`` and ``--output``
paths. The source directory must be a copy of the skin_drug_target_mr project
containing the seven compressed GTEx and FinnGen inputs used by the analysis.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


SOURCE = Path.cwd()
OUT = Path(__file__).resolve().parent / "reanalysis"

PRIMARY_P = 5e-8
RELAXED_P = 1e-5
MIN_F = 10.0
CHUNK_SIZE = 1_000_000

TISSUES: dict[str, tuple[str, Path]] = {}
OUTCOMES: dict[str, tuple[str, Path]] = {}


def configure_paths(source: Path, output: Path) -> None:
    global SOURCE, OUT, TISSUES, OUTCOMES
    SOURCE = source.expanduser().resolve()
    OUT = output.expanduser().resolve()
    TISSUES = {
        "Skin_Not_Sun_Exposed_Suprapubic": (
            "Skin not sun-exposed suprapubic",
            SOURCE
            / "data/raw/gtex/GTEx_Analysis_v8_eQTL"
            / "Skin_Not_Sun_Exposed_Suprapubic.v8.signif_variant_gene_pairs.txt.gz",
        ),
        "Skin_Sun_Exposed_Lower_leg": (
            "Skin sun-exposed lower leg",
            SOURCE
            / "data/raw/gtex/GTEx_Analysis_v8_eQTL"
            / "Skin_Sun_Exposed_Lower_leg.v8.signif_variant_gene_pairs.txt.gz",
        ),
        "Whole_Blood": (
            "Whole blood",
            SOURCE
            / "data/raw/gtex/GTEx_Analysis_v8_eQTL"
            / "Whole_Blood.v8.signif_variant_gene_pairs.txt.gz",
        ),
    }
    OUTCOMES = {
        "alopecia_areata": (
            "L12_ALOPECAREATA",
            SOURCE / "data/raw/finngen/outcomes/finngen_R13_L12_ALOPECAREATA.gz",
        ),
        "atopic_dermatitis": (
            "L12_DERMATITISECZEMA",
            SOURCE / "data/raw/finngen/outcomes/finngen_R13_L12_DERMATITISECZEMA.gz",
        ),
        "psoriasis": (
            "L12_PSORIASIS",
            SOURCE / "data/raw/finngen/outcomes/finngen_R13_L12_PSORIASIS.gz",
        ),
        "vitiligo": (
            "L12_VITILIGO",
            SOURCE / "data/raw/finngen/outcomes/finngen_R13_L12_VITILIGO.gz",
        ),
    }


def normal_two_sided_p(z_value: float) -> float:
    return math.erfc(abs(z_value) / math.sqrt(2.0))


def bh_adjust(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    adjusted = [math.nan] * len(values)
    running = 1.0
    total = len(values)
    for rank_index in range(total - 1, -1, -1):
        original_index, value = indexed[rank_index]
        rank = rank_index + 1
        running = min(running, value * total / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def parse_variant_id(value: str) -> tuple[str, int, str, str]:
    parts = value.split("_")
    if len(parts) < 4:
        raise ValueError(f"Unexpected GTEx variant_id: {value}")
    chrom = parts[0].removeprefix("chr").upper()
    return chrom, int(parts[1]), parts[2].upper(), parts[3].upper()


def variant_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{str(chrom).removeprefix('chr').upper()}:{int(pos)}:{ref.upper()}:{alt.upper()}"


def input_metadata(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {
        "path": path.relative_to(SOURCE).as_posix(),
        "bytes": stat.st_size,
        "modified_ns": stat.st_mtime_ns,
    }


def read_target_map() -> tuple[dict[str, str], list[str]]:
    annotation_path = SOURCE / "data/processed/target_gene_annotation_gencode_v26.csv"
    annotation = pd.read_csv(annotation_path, dtype=str)
    target_path = SOURCE / "data/target_genes/skin_druggable_targets.csv"
    targets = pd.read_csv(target_path, dtype=str)["gene"].str.upper().tolist()
    filtered = annotation[annotation["gene"].str.upper().isin(targets)].copy()
    mapping = dict(zip(filtered["gene_id"].str.replace(r"\.[0-9]+$", "", regex=True), filtered["gene"].str.upper()))
    if len(mapping) != len(targets):
        raise RuntimeError(f"Mapped {len(mapping)} of {len(targets)} targets")
    return mapping, targets


def scan_gtex(mapping: dict[str, str]) -> tuple[pd.DataFrame, dict[str, int]]:
    usecols = [
        "variant_id",
        "gene_id",
        "ma_samples",
        "maf",
        "pval_nominal",
        "slope",
        "slope_se",
    ]
    selected_chunks: list[pd.DataFrame] = []
    rows_scanned: dict[str, int] = {}
    wanted_ids = set(mapping)

    for tissue, (label, path) in TISSUES.items():
        scanned = 0
        tissue_chunks: list[pd.DataFrame] = []
        for chunk in pd.read_csv(
            path,
            sep="\t",
            compression="gzip",
            usecols=usecols,
            chunksize=CHUNK_SIZE,
            dtype={"variant_id": "string", "gene_id": "string"},
        ):
            scanned += len(chunk)
            stripped = chunk["gene_id"].str.replace(r"\.[0-9]+$", "", regex=True)
            keep = stripped.isin(wanted_ids)
            if not keep.any():
                continue
            subset = chunk.loc[keep].copy()
            subset["gene_id_stripped"] = stripped.loc[keep]
            subset["target_gene"] = subset["gene_id_stripped"].map(mapping)
            subset["f_stat"] = (subset["slope"] / subset["slope_se"]) ** 2
            subset = subset[
                subset["pval_nominal"].notna()
                & subset["f_stat"].notna()
                & (subset["pval_nominal"] < RELAXED_P)
                & (subset["f_stat"] > MIN_F)
            ]
            if not subset.empty:
                tissue_chunks.append(subset)

        rows_scanned[tissue] = scanned
        if tissue_chunks:
            tissue_df = pd.concat(tissue_chunks, ignore_index=True)
            tissue_df["tissue"] = tissue
            tissue_df["tissue_label"] = label
            selected_chunks.append(tissue_df)

    if not selected_chunks:
        raise RuntimeError("No eligible GTEx rows found")
    return pd.concat(selected_chunks, ignore_index=True), rows_scanned


def select_instruments(
    candidates: pd.DataFrame, targets: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    selected: list[dict[str, object]] = []
    availability: list[dict[str, object]] = []

    for gene in targets:
        for tissue, (label, _) in TISSUES.items():
            cell = candidates[
                (candidates["target_gene"] == gene) & (candidates["tissue"] == tissue)
            ].sort_values(["pval_nominal", "variant_id"], kind="mergesort")
            primary = cell[cell["pval_nominal"] < PRIMARY_P]
            if not primary.empty:
                lead = primary.iloc[0]
                tier = "primary"
                status = "primary_p_lt_5e-8"
            elif not cell.empty:
                lead = cell.iloc[0]
                tier = "exploratory_relaxed"
                status = "exploratory_5e-8_to_1e-5"
            else:
                lead = None
                tier = "none"
                status = "no_eligible_row_in_significant_pair_file"

            row: dict[str, object] = {
                "target_gene": gene,
                "tissue": tissue,
                "tissue_label": label,
                "instrument_status": status,
                "variant_id": "",
                "chrom": "",
                "pos": np.nan,
                "ref": "",
                "alt": "",
                "beta_exposure": np.nan,
                "se_exposure": np.nan,
                "p_exposure": np.nan,
                "maf": np.nan,
                "ma_samples": np.nan,
                "f_stat": np.nan,
            }
            if lead is not None:
                chrom, pos, ref, alt = parse_variant_id(str(lead["variant_id"]))
                row.update(
                    {
                        "variant_id": str(lead["variant_id"]),
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "beta_exposure": float(lead["slope"]),
                        "se_exposure": float(lead["slope_se"]),
                        "p_exposure": float(lead["pval_nominal"]),
                        "maf": float(lead["maf"]),
                        "ma_samples": float(lead["ma_samples"]),
                        "f_stat": float(lead["f_stat"]),
                        "tier": tier,
                    }
                )
                selected.append(row.copy())
            availability.append(row)

    selected_df = pd.DataFrame(selected)
    primary_df = selected_df[selected_df["tier"] == "primary"].reset_index(drop=True)
    relaxed_df = selected_df[selected_df["tier"] == "exploratory_relaxed"].reset_index(drop=True)
    return primary_df, relaxed_df, pd.DataFrame(availability)


def scan_finngen(
    instruments: pd.DataFrame,
) -> tuple[dict[str, dict[str, object]], dict[str, int]]:
    positions = {
        (str(row.chrom).upper(), int(row.pos))
        for row in instruments.itertuples(index=False)
    }
    position_values = {pos for _, pos in positions}
    matches: dict[str, dict[str, object]] = {}
    rows_scanned: dict[str, int] = {}
    usecols = ["#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"]

    for outcome_name, (endpoint, path) in OUTCOMES.items():
        scanned = 0
        found = 0
        for chunk in pd.read_csv(
            path,
            sep="\t",
            compression="gzip",
            usecols=usecols,
            chunksize=CHUNK_SIZE,
            dtype={"#chrom": "string", "ref": "string", "alt": "string", "rsids": "string"},
        ):
            scanned += len(chunk)
            subset = chunk[chunk["pos"].isin(position_values)].copy()
            if subset.empty:
                continue
            subset["chrom_norm"] = subset["#chrom"].str.removeprefix("chr").str.upper()
            subset = subset[
                [
                    (chrom, int(pos)) in positions
                    for chrom, pos in zip(subset["chrom_norm"], subset["pos"])
                ]
            ]
            for row in subset.itertuples(index=False):
                key = variant_key(row.chrom_norm, row.pos, row.ref, row.alt)
                matches[f"{outcome_name}|{key}"] = {
                    "endpoint": endpoint,
                    "rsid": "" if pd.isna(row.rsids) else str(row.rsids),
                    "beta_outcome": float(row.beta),
                    "se_outcome": float(row.sebeta),
                    "p_outcome": float(row.pval),
                    "eaf_outcome": float(row.af_alt),
                    "match_orientation": "forward",
                }
                found += 1
        rows_scanned[outcome_name] = scanned
        print(f"FinnGen {outcome_name}: scanned {scanned:,} rows; retained {found}", flush=True)
    return matches, rows_scanned


def build_mr_results(
    instruments: pd.DataFrame, matches: dict[str, dict[str, object]]
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    missing: list[dict[str, object]] = []

    for instrument in instruments.itertuples(index=False):
        forward_key = variant_key(instrument.chrom, instrument.pos, instrument.ref, instrument.alt)
        reverse_key = variant_key(instrument.chrom, instrument.pos, instrument.alt, instrument.ref)
        for outcome_name, (endpoint, _) in OUTCOMES.items():
            match = matches.get(f"{outcome_name}|{forward_key}")
            orientation = "forward"
            if match is None:
                match = matches.get(f"{outcome_name}|{reverse_key}")
                orientation = "reversed"
            if match is None:
                missing.append(
                    {
                        "target_gene": instrument.target_gene,
                        "tissue": instrument.tissue,
                        "variant_id": instrument.variant_id,
                        "outcome": outcome_name,
                        "endpoint": endpoint,
                    }
                )
                continue

            beta_outcome = float(match["beta_outcome"])
            if orientation == "reversed":
                beta_outcome = -beta_outcome
            beta = beta_outcome / float(instrument.beta_exposure)
            se_first = abs(float(match["se_outcome"]) / float(instrument.beta_exposure))
            se_delta = math.sqrt(
                (float(match["se_outcome"]) ** 2 / float(instrument.beta_exposure) ** 2)
                + (
                    beta_outcome**2
                    * float(instrument.se_exposure) ** 2
                    / float(instrument.beta_exposure) ** 4
                )
            )
            row = {
                "target_gene": instrument.target_gene,
                "tissue": instrument.tissue,
                "tissue_label": instrument.tissue_label,
                "outcome": outcome_name,
                "endpoint": endpoint,
                "tier": instrument.tier,
                "method": "Wald ratio",
                "nsnp": 1,
                "variant_id": instrument.variant_id,
                "rsid": match["rsid"],
                "chrom": instrument.chrom,
                "pos": instrument.pos,
                "effect_allele_exposure": instrument.alt,
                "other_allele_exposure": instrument.ref,
                "beta_exposure": instrument.beta_exposure,
                "se_exposure": instrument.se_exposure,
                "p_exposure": instrument.p_exposure,
                "f_stat": instrument.f_stat,
                "beta_outcome": beta_outcome,
                "se_outcome": match["se_outcome"],
                "p_outcome": match["p_outcome"],
                "eaf_outcome": match["eaf_outcome"],
                "harmonisation": orientation,
                "beta_mr": beta,
                "se_mr": se_first,
                "p_value": normal_two_sided_p(beta / se_first),
                "OR": math.exp(beta),
                "CI_lower": math.exp(beta - 1.96 * se_first),
                "CI_upper": math.exp(beta + 1.96 * se_first),
                "se_full_delta": se_delta,
                "p_full_delta": normal_two_sided_p(beta / se_delta),
            }
            rows.append(row)

    result = pd.DataFrame(rows)
    if result.empty:
        return result, missing
    result["FDR"] = bh_adjust(result["p_value"].tolist())
    result["FDR_full_delta"] = bh_adjust(result["p_full_delta"].tolist())
    result["significant_category"] = np.where(
        result["FDR"] < 0.05,
        "FDR_significant",
        np.where(result["p_value"] < 0.05, "nominal", "null"),
    )
    result["full_delta_category"] = np.where(
        result["FDR_full_delta"] < 0.05,
        "FDR_significant",
        np.where(result["p_full_delta"] < 0.05, "nominal", "null"),
    )
    return result, missing


def compare_saved(primary: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    public_path = SOURCE / "results/table_s3_primary_mr_68.csv"
    legacy_path = SOURCE / "results/mr_all_results.csv"
    saved_path = public_path if public_path.exists() else legacy_path
    if not saved_path.exists():
        raise FileNotFoundError(
            "Expected results/table_s3_primary_mr_68.csv or results/mr_all_results.csv"
        )
    saved = pd.read_csv(saved_path)
    if saved_path == public_path:
        saved_columns = {
            "beta_mr": "saved_beta",
            "se_mr": "saved_se",
            "p_value": "saved_p",
            "FDR": "saved_fdr",
        }
    else:
        saved_columns = {
            "b": "saved_beta",
            "se": "saved_se",
            "pval": "saved_p",
            "fdr": "saved_fdr",
        }
    old = saved.rename(columns=saved_columns)[
        ["target_gene", "tissue", "endpoint", "saved_beta", "saved_se", "saved_p", "saved_fdr"]
    ]
    merged = primary.merge(old, on=["target_gene", "tissue", "endpoint"], how="outer", indicator=True)
    for new_col, old_col, label in (
        ("beta_mr", "saved_beta", "beta_abs_error"),
        ("se_mr", "saved_se", "se_abs_error"),
        ("p_value", "saved_p", "p_abs_error"),
        ("FDR", "saved_fdr", "fdr_abs_error"),
    ):
        merged[label] = (merged[new_col] - merged[old_col]).abs()
    summary = {
        "saved_source": str(saved_path.relative_to(SOURCE)),
        "saved_rows": len(saved),
        "reproduced_rows": len(primary),
        "matched_rows": int((merged["_merge"] == "both").sum()),
        "only_reproduced": int((merged["_merge"] == "left_only").sum()),
        "only_saved": int((merged["_merge"] == "right_only").sum()),
        "max_abs_errors": {
            column: (None if merged[column].dropna().empty else float(merged[column].max()))
            for column in ("beta_abs_error", "se_abs_error", "p_abs_error", "fdr_abs_error")
        },
    }
    return merged, summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=Path.cwd(),
        help="skin_drug_target_mr project root containing data/raw and results",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "reanalysis",
        help="directory for regenerated CSV and JSON audit outputs",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configure_paths(args.source, args.output)
    if not (SOURCE / "analysis_config.yaml").exists():
        raise FileNotFoundError(
            f"--source does not look like the project root: {SOURCE}"
        )
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"Source project: {SOURCE}", flush=True)
    print(f"Audit output: {OUT}", flush=True)
    for _, path in TISSUES.values():
        if not path.exists():
            raise FileNotFoundError(path)
    for _, path in OUTCOMES.values():
        if not path.exists():
            raise FileNotFoundError(path)

    mapping, targets = read_target_map()
    print("Scanning GTEx raw files", flush=True)
    candidates, gtex_rows = scan_gtex(mapping)
    primary_instruments, relaxed_instruments, availability = select_instruments(candidates, targets)
    all_instruments = pd.concat([primary_instruments, relaxed_instruments], ignore_index=True)

    primary_instruments.to_csv(OUT / "primary_instruments_from_raw.csv", index=False)
    relaxed_instruments.to_csv(OUT / "exploratory_relaxed_instruments_from_raw.csv", index=False)
    availability.to_csv(OUT / "instrument_availability_57_from_raw.csv", index=False)
    candidates.sort_values(["target_gene", "tissue", "pval_nominal"]).to_csv(
        OUT / "eligible_gtex_candidates_from_raw.csv", index=False
    )
    print(
        f"GTEx: {len(primary_instruments)} primary cells; "
        f"{len(relaxed_instruments)} relaxed-only cells",
        flush=True,
    )

    matches, finngen_rows = scan_finngen(all_instruments)
    primary_results, primary_missing = build_mr_results(primary_instruments, matches)
    relaxed_results, relaxed_missing = build_mr_results(relaxed_instruments, matches)
    primary_results.to_csv(OUT / "primary_mr_results_from_raw.csv", index=False)
    relaxed_results.to_csv(OUT / "exploratory_relaxed_mr_results_from_raw.csv", index=False)

    comparison, comparison_summary = compare_saved(primary_results)
    comparison.to_csv(OUT / "comparison_with_saved_68.csv", index=False)

    primary_hits = primary_results[primary_results["FDR"] < 0.05]
    primary_nominal = primary_results[
        (primary_results["p_value"] < 0.05) & (primary_results["FDR"] >= 0.05)
    ]
    delta_hits = primary_results[primary_results["FDR_full_delta"] < 0.05]
    relaxed_hits = relaxed_results[relaxed_results["FDR"] < 0.05]

    audit = {
        "analysis": "single-SNP lead cis-eQTL Wald-ratio MR",
        "thresholds": {
            "primary_p": PRIMARY_P,
            "exploratory_relaxed_p": RELAXED_P,
            "minimum_f_statistic": MIN_F,
            "lead_variant_tie_break": "lowest variant_id after p-value ordering",
            "multiple_testing": "BH separately within primary and exploratory families",
        },
        "input_files": [
            input_metadata(path)
            for _, path in list(TISSUES.values()) + list(OUTCOMES.values())
        ],
        "gtex_rows_scanned": gtex_rows,
        "finngen_rows_scanned": finngen_rows,
        "target_gene_tissue_cells": len(targets) * len(TISSUES),
        "primary_instrument_cells": len(primary_instruments),
        "exploratory_relaxed_only_cells": len(relaxed_instruments),
        "cells_with_any_p_lt_1e5_instrument": len(all_instruments),
        "primary_mr_rows": len(primary_results),
        "primary_missing_outcome_matches": primary_missing,
        "exploratory_mr_rows": len(relaxed_results),
        "exploratory_missing_outcome_matches": relaxed_missing,
        "saved_68_comparison": comparison_summary,
        "primary_fdr_hits": primary_hits[
            ["target_gene", "tissue", "outcome", "OR", "CI_lower", "CI_upper", "p_value", "FDR"]
        ].to_dict(orient="records"),
        "primary_nominal_hits": primary_nominal[
            ["target_gene", "tissue", "outcome", "OR", "CI_lower", "CI_upper", "p_value", "FDR"]
        ].to_dict(orient="records"),
        "full_delta_fdr_hits": delta_hits[
            ["target_gene", "tissue", "outcome", "OR", "p_full_delta", "FDR_full_delta"]
        ].to_dict(orient="records"),
        "exploratory_relaxed_fdr_hits": relaxed_hits[
            ["target_gene", "tissue", "outcome", "OR", "CI_lower", "CI_upper", "p_value", "FDR"]
        ].to_dict(orient="records"),
    }
    (OUT / "reanalysis_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(audit, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
