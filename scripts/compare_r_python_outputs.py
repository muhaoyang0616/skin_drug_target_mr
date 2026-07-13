#!/usr/bin/env python3
"""Compare corrected R MR outputs with the independent Python reproduction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


KEYS = ["target_gene", "tissue", "outcome"]
METRICS = {
    "beta": ("beta_mr", "b"),
    "se": ("se_mr", "se"),
    "p": ("p_value", "pval"),
    "fdr": ("FDR", "fdr"),
}


def compare_family(
    python_path: Path,
    r_path: Path,
    output_path: Path,
) -> dict[str, object]:
    python_df = pd.read_csv(python_path)
    r_df = pd.read_csv(r_path).rename(columns={"disease_query": "outcome"})
    r_df = r_df[(r_df["method"] == "Wald ratio") & r_df["method_is_primary"]]

    python_cols = KEYS + [pair[0] for pair in METRICS.values()]
    r_cols = KEYS + [pair[1] for pair in METRICS.values()]
    merged = python_df[python_cols].merge(
        r_df[r_cols],
        on=KEYS,
        how="outer",
        indicator=True,
        validate="one_to_one",
        suffixes=("_python", "_r"),
    )

    max_errors: dict[str, float | None] = {}
    for label, (python_col, r_col) in METRICS.items():
        error_col = f"{label}_abs_error"
        merged[error_col] = (merged[python_col] - merged[r_col]).abs()
        max_errors[error_col] = (
            float(merged[error_col].max()) if merged[error_col].notna().any() else None
        )

    merged.to_csv(output_path, index=False)
    return {
        "python_rows": int(len(python_df)),
        "r_rows": int(len(r_df)),
        "matched_rows": int((merged["_merge"] == "both").sum()),
        "python_only_rows": int((merged["_merge"] == "left_only").sum()),
        "r_only_rows": int((merged["_merge"] == "right_only").sum()),
        "max_abs_errors": max_errors,
    }


def compare_instruments(
    python_dir: Path,
    r_project: Path,
    output_path: Path,
) -> dict[str, object]:
    python_df = pd.concat(
        [
            pd.read_csv(python_dir / "primary_instruments_from_raw.csv"),
            pd.read_csv(python_dir / "exploratory_relaxed_instruments_from_raw.csv"),
        ],
        ignore_index=True,
    )
    r_df = pd.read_csv(r_project / "data/processed/instruments_gtex.csv")
    r_df = r_df.rename(
        columns={
            "pval.exposure": "p_exposure_r",
            "beta.exposure": "beta_exposure_r",
            "se.exposure": "se_exposure_r",
            "f_stat": "f_stat_r",
            "instrument_tier": "tier_r",
        }
    )
    python_df = python_df.rename(
        columns={
            "p_exposure": "p_exposure_python",
            "beta_exposure": "beta_exposure_python",
            "se_exposure": "se_exposure_python",
            "f_stat": "f_stat_python",
            "tier": "tier_python",
        }
    )
    keys = ["target_gene", "tissue", "variant_id"]
    merged = python_df[
        keys
        + [
            "tier_python",
            "beta_exposure_python",
            "se_exposure_python",
            "p_exposure_python",
            "f_stat_python",
        ]
    ].merge(
        r_df[
            keys
            + [
                "tier_r",
                "beta_exposure_r",
                "se_exposure_r",
                "p_exposure_r",
                "f_stat_r",
            ]
        ],
        on=keys,
        how="outer",
        indicator=True,
        validate="one_to_one",
    )
    for label in ("beta_exposure", "se_exposure", "p_exposure", "f_stat"):
        merged[f"{label}_abs_error"] = (
            merged[f"{label}_python"] - merged[f"{label}_r"]
        ).abs()
    merged["tier_match"] = merged["tier_python"] == merged["tier_r"]
    merged.to_csv(output_path, index=False)
    error_columns = [
        "beta_exposure_abs_error",
        "se_exposure_abs_error",
        "p_exposure_abs_error",
        "f_stat_abs_error",
    ]
    return {
        "python_rows": int(len(python_df)),
        "r_rows": int(len(r_df)),
        "matched_rows": int((merged["_merge"] == "both").sum()),
        "python_only_rows": int((merged["_merge"] == "left_only").sum()),
        "r_only_rows": int((merged["_merge"] == "right_only").sum()),
        "tier_mismatches": int((~merged["tier_match"].fillna(False)).sum()),
        "max_abs_errors": {
            column: float(merged[column].max()) for column in error_columns
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python-dir", type=Path, required=True)
    parser.add_argument("--r-project", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    summary = {
        "instruments": compare_instruments(
            args.python_dir,
            args.r_project,
            args.output / "r_python_instrument_comparison.csv",
        ),
        "primary": compare_family(
            args.python_dir / "primary_mr_results_from_raw.csv",
            args.r_project / "results/mr_all_results.csv",
            args.output / "r_python_primary_comparison.csv",
        ),
        "exploratory_relaxed": compare_family(
            args.python_dir / "exploratory_relaxed_mr_results_from_raw.csv",
            args.r_project / "results/mr_exploratory_relaxed_results.csv",
            args.output / "r_python_exploratory_comparison.csv",
        ),
    }
    output_json = args.output / "r_python_comparison_summary.json"
    output_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
