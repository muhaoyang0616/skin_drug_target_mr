#!/usr/bin/env python3
"""Build manuscript figures from the verified reanalysis tables."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PUBLIC_ROOT = SCRIPT_DIR.parent

if (PUBLIC_ROOT / "results" / "table_s3_primary_mr_68.csv").is_file():
    RESULTS = PUBLIC_ROOT / "results"
    PRIMARY_PATH = RESULTS / "table_s3_primary_mr_68.csv"
    AVAILABILITY_PATH = RESULTS / "table_s2_instrument_availability_57.csv"
    OPENGWAS_PATH = RESULTS / "table_s6_opengwas_validation_23.csv"
    COLOC_PATH = RESULTS / "table_s7_colocalization.csv"
    OUT = PUBLIC_ROOT / "manuscript" / "figures"
else:
    reanalysis = SCRIPT_DIR / "reanalysis_portable_20260711"
    source_project = SCRIPT_DIR.parents[1] / "skin_drug_target_mr"
    PRIMARY_PATH = reanalysis / "primary_mr_results_from_raw.csv"
    AVAILABILITY_PATH = reanalysis / "instrument_availability_57_from_raw.csv"
    OPENGWAS_PATH = source_project / "results" / "validation" / "tyk2_opengwas_replication.csv"
    COLOC_PATH = source_project / "manuscript_assets" / "table_eqtlcatalogue_coloc.csv"
    OUT = SCRIPT_DIR.parent / "manuscript_assets" / "figures"

COLORS = {
    "ink": "#20262E",
    "muted": "#66717E",
    "grid": "#D6DCE2",
    "primary": "#007C83",
    "exploratory": "#D89B2B",
    "risk": "#B84A3A",
    "protective": "#2878B5",
    "green": "#3B8C6E",
    "panel": "#F5F7F8",
}

TARGET_ORDER = [
    "JAK1",
    "JAK2",
    "JAK3",
    "TYK2",
    "IL4R",
    "IL13",
    "IL17A",
    "IL17RA",
    "IL23A",
    "IL23R",
    "TNF",
    "PDE4A",
    "PDE4B",
    "PDE4C",
    "PDE4D",
    "IL36RN",
    "TSLP",
    "FLG",
    "CARD14",
]

TISSUE_SHORT = {
    "Skin_Sun_Exposed_Lower_leg": "Sun skin",
    "Skin_Not_Sun_Exposed_Suprapubic": "Non-sun skin",
    "Whole_Blood": "Blood",
}

OUTCOME_LABEL = {
    "psoriasis": "Psoriasis",
    "atopic_dermatitis": "Dermatitis/eczema",
    "vitiligo": "Vitiligo",
    "alopecia_areata": "Alopecia areata",
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 11,
            "axes.titlesize": 14,
            "axes.labelsize": 11,
            "axes.edgecolor": COLORS["ink"],
            "axes.linewidth": 0.8,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def add_box(ax, xy, width, height, title, lines, color, title_color="white"):
    x, y = xy
    box = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.015",
        linewidth=1.2,
        edgecolor=color,
        facecolor="white",
    )
    ax.add_patch(box)
    header = FancyBboxPatch(
        (x, y + height - 0.10),
        width,
        0.10,
        boxstyle="round,pad=0.012,rounding_size=0.015",
        linewidth=0,
        facecolor=color,
    )
    ax.add_patch(header)
    ax.text(x + width / 2, y + height - 0.05, title, ha="center", va="center", color=title_color, weight="bold")
    ax.text(
        x + 0.025,
        y + height - 0.135,
        "\n".join(lines),
        ha="left",
        va="top",
        color=COLORS["ink"],
        fontsize=10,
        linespacing=1.18,
    )


def build_figure1() -> None:
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(
        0.5,
        0.97,
        "Tissue-aware cis-eQTL MR workflow",
        ha="center",
        va="top",
        fontsize=18,
        weight="bold",
        color=COLORS["ink"],
    )

    add_box(
        ax,
        (0.05, 0.66),
        0.25,
        0.23,
        "Target panel",
        ["19 pathway genes", "3 GTEx tissues", "57 gene-tissue cells"],
        COLORS["ink"],
    )
    add_box(
        ax,
        (0.375, 0.66),
        0.25,
        0.23,
        "Primary instruments",
        ["eQTL P < 5e-8; F > 10", "18 instrumented cells", "17 matched to FinnGen"],
        COLORS["primary"],
    )
    add_box(
        ax,
        (0.70, 0.66),
        0.25,
        0.23,
        "Primary MR family",
        ["4 FinnGen outcomes", "68 Wald-ratio tests", "BH-FDR across 68"],
        COLORS["protective"],
    )

    add_box(
        ax,
        (0.14, 0.31),
        0.31,
        0.24,
        "Primary result",
        ["3 FDR-significant rows", "All TYK2 to psoriasis", "Delta sensitivity unchanged"],
        COLORS["green"],
    )
    add_box(
        ax,
        (0.55, 0.31),
        0.31,
        0.24,
        "TYK2 validation",
        ["GTEx eQTL rank and FinnGen locus", "8 external GWAS records", "Tissue-resolved colocalization"],
        COLORS["risk"],
    )

    add_box(
        ax,
        (0.23, 0.035),
        0.54,
        0.18,
        "Separate exploratory family",
        [
            "5 relaxed-only cells (P < 1e-5; F > 10) and 20 tests",
            "PDE4A locus-level analysis non-supportive",
        ],
        COLORS["exploratory"],
        title_color=COLORS["ink"],
    )

    arrow = dict(arrowstyle="-|>", color=COLORS["muted"], lw=1.8, mutation_scale=14)
    ax.annotate("", xy=(0.365, 0.775), xytext=(0.31, 0.775), arrowprops=arrow)
    ax.annotate("", xy=(0.69, 0.775), xytext=(0.635, 0.775), arrowprops=arrow)
    ax.annotate("", xy=(0.30, 0.56), xytext=(0.76, 0.655), arrowprops=arrow)
    ax.annotate("", xy=(0.54, 0.44), xytext=(0.46, 0.44), arrowprops=arrow)
    ax.plot([0.175, 0.025, 0.025, 0.215], [0.655, 0.655, 0.125, 0.125], color=COLORS["muted"], lw=1.8)
    ax.annotate("", xy=(0.23, 0.125), xytext=(0.21, 0.125), arrowprops=arrow)
    save(fig, "figure1_study_design_flowchart")


def build_figure2(primary: pd.DataFrame) -> None:
    cells = (
        primary[["target_gene", "tissue"]]
        .drop_duplicates()
        .assign(
            gene_order=lambda d: d["target_gene"].map({g: i for i, g in enumerate(TARGET_ORDER)}),
            tissue_order=lambda d: d["tissue"].map(
                {
                    "Skin_Sun_Exposed_Lower_leg": 0,
                    "Skin_Not_Sun_Exposed_Suprapubic": 1,
                    "Whole_Blood": 2,
                }
            ),
        )
        .sort_values(["gene_order", "tissue_order"])
    )
    labels = [f"{row.target_gene} | {TISSUE_SHORT[row.tissue]}" for row in cells.itertuples()]
    cell_y = {(row.target_gene, row.tissue): i for i, row in enumerate(cells.itertuples())}
    outcomes = ["psoriasis", "atopic_dermatitis", "vitiligo", "alopecia_areata"]
    outcome_x = {name: i for i, name in enumerate(outcomes)}

    fig, ax = plt.subplots(figsize=(10.5, 10.5))
    norm = mcolors.TwoSlopeNorm(vmin=-1.2, vcenter=0, vmax=1.2)
    cmap = plt.get_cmap("RdBu_r")

    for category, marker, edge, linewidth in [
        ("null", "o", "#9AA5AF", 0.8),
        ("nominal", "D", COLORS["ink"], 1.3),
        ("FDR_significant", "*", COLORS["ink"], 1.5),
    ]:
        subset = primary[primary["significant_category"] == category]
        if subset.empty:
            continue
        x = [outcome_x[value] for value in subset["outcome"]]
        y = [cell_y[(g, t)] for g, t in zip(subset["target_gene"], subset["tissue"])]
        strength = np.clip(-np.log10(subset["p_value"].astype(float)), 0, 14)
        size = 35 + 20 * strength
        ax.scatter(
            x,
            y,
            c=subset["beta_mr"].astype(float),
            cmap=cmap,
            norm=norm,
            s=size,
            marker=marker,
            edgecolors=edge,
            linewidths=linewidth,
            alpha=1.0,
            zorder=3,
        )

    ax.set_xticks(range(len(outcomes)), [OUTCOME_LABEL[value] for value in outcomes])
    ax.set_yticks(range(len(labels)), labels)
    ax.invert_yaxis()
    ax.set_xlim(-0.5, len(outcomes) - 0.5)
    ax.set_ylim(len(labels) - 0.5, -0.5)
    ax.grid(which="major", color=COLORS["grid"], linewidth=0.7)
    ax.set_axisbelow(True)
    ax.set_title("Primary MR results (68 tests)", loc="left", weight="bold", pad=30)
    ax.text(
        0,
        1.005,
        "Color: log odds ratio   Size: -log10(P)   Shape: FDR / nominal / null",
        transform=ax.transAxes,
        color=COLORS["muted"],
        va="bottom",
    )

    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Wald-ratio log OR (clipped at +/-1.2)")
    handles = [
        Line2D([0], [0], marker="*", color="none", markerfacecolor="#7FA9C5", markeredgecolor=COLORS["ink"], markersize=12, label="FDR < 0.05"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor="#7FA9C5", markeredgecolor=COLORS["ink"], markersize=8, label="Nominal P < 0.05"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#7FA9C5", markeredgecolor="white", markersize=8, label="Null"),
    ]
    ax.legend(
        handles=handles,
        loc="lower right",
        frameon=True,
        framealpha=1.0,
        edgecolor=COLORS["grid"],
    )
    fig.subplots_adjust(left=0.29, right=0.91, top=0.91, bottom=0.09)
    save(fig, "figure2_mr_results")


def build_figure3(availability: pd.DataFrame) -> None:
    tissue_order = [
        "Skin_Sun_Exposed_Lower_leg",
        "Skin_Not_Sun_Exposed_Suprapubic",
        "Whole_Blood",
    ]
    status_value = {
        "no_eligible_row_in_significant_pair_file": 0,
        "exploratory_5e-8_to_1e-5": 1,
        "primary_p_lt_5e-8": 2,
    }
    matrix = np.zeros((len(TARGET_ORDER), len(tissue_order)))
    text_matrix = np.full(matrix.shape, "--", dtype=object)
    for row in availability.itertuples(index=False):
        i = TARGET_ORDER.index(row.target_gene)
        j = tissue_order.index(row.tissue)
        value = status_value[row.instrument_status]
        matrix[i, j] = value
        text_matrix[i, j] = {0: "--", 1: "E", 2: "P"}[value]
    text_matrix[TARGET_ORDER.index("PDE4A"), tissue_order.index("Whole_Blood")] = "P*"

    cmap = mcolors.ListedColormap(["#F2F4F5", COLORS["exploratory"], COLORS["primary"]])
    norm = mcolors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    fig, ax = plt.subplots(figsize=(8.8, 9.2))
    ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks(
        range(3),
        ["Sun-exposed\nskin", "Non-sun-exposed\nskin", "Whole\nblood"],
    )
    ax.set_yticks(range(len(TARGET_ORDER)), TARGET_ORDER)
    ax.tick_params(axis="x", rotation=0, pad=8)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            color = "white" if matrix[i, j] == 2 else COLORS["ink"]
            ax.text(j, i, text_matrix[i, j], ha="center", va="center", color=color, weight="bold")
    ax.set_xticks(np.arange(-0.5, 3, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(TARGET_ORDER), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_title("Cis-eQTL instrument availability (57 cells)", loc="left", weight="bold", pad=30)
    ax.text(
        0,
        1.005,
        "Primary: 18 cells   |   Relaxed-only: 5 cells   |   Any eligible instrument: 23 cells",
        transform=ax.transAxes,
        color=COLORS["muted"],
        va="bottom",
    )
    handles = [
        Line2D([0], [0], marker="s", color="none", markerfacecolor=COLORS["primary"], markersize=11, label="P: primary"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=COLORS["exploratory"], markersize=11, label="E: relaxed-only"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#F2F4F5", markeredgecolor=COLORS["grid"], markersize=11, label="No eligible row"),
    ]
    ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, -0.18), ncol=3, frameon=False)
    fig.subplots_adjust(left=0.22, right=0.98, top=0.91, bottom=0.21)
    save(fig, "figure3_instrument_availability_heatmap")


def build_figure4(primary: pd.DataFrame) -> None:
    tyk2 = primary[(primary["target_gene"] == "TYK2") & (primary["outcome"] == "psoriasis")].copy()
    tissue_order = [
        "Skin_Sun_Exposed_Lower_leg",
        "Skin_Not_Sun_Exposed_Suprapubic",
        "Whole_Blood",
    ]
    tyk2["order"] = tyk2["tissue"].map({v: i for i, v in enumerate(tissue_order)})
    tyk2 = tyk2.sort_values("order")

    coloc = pd.read_csv(COLOC_PATH)
    coloc = coloc[coloc["gene"] == "TYK2"].copy()
    coloc_map = {
        "skin_sun_exposed_lower_leg_skin": "Skin_Sun_Exposed_Lower_leg",
        "skin_not_sun_exposed_suprapubic": "Skin_Not_Sun_Exposed_Suprapubic",
        "whole_blood": "Whole_Blood",
    }
    coloc["tissue_key"] = coloc["tissue_or_cell_type"].map(coloc_map)
    coloc["order"] = coloc["tissue_key"].map({v: i for i, v in enumerate(tissue_order)})
    coloc = coloc.sort_values("order")

    fig = plt.figure(figsize=(11.5, 8.2))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.0, 1.2],
        height_ratios=[1.0, 1.05],
        wspace=0.34,
        hspace=0.48,
    )

    ax1 = fig.add_subplot(gs[0, 0])
    y = np.arange(3)
    ax1.errorbar(
        tyk2["OR"],
        y,
        xerr=[tyk2["OR"] - tyk2["CI_lower"], tyk2["CI_upper"] - tyk2["OR"]],
        fmt="o",
        color=COLORS["protective"],
        ecolor=COLORS["protective"],
        capsize=4,
        markersize=8,
        linewidth=2,
    )
    ax1.axvline(1, color=COLORS["muted"], linestyle="--", linewidth=1)
    ax1.set_xscale("log")
    ax1.set_xlim(0.32, 1.08)
    ax1.set_yticks(y, [TISSUE_SHORT[v] for v in tissue_order])
    ax1.invert_yaxis()
    ax1.set_xlabel("OR for psoriasis (95% CI)")
    ax1.set_title("A  Primary MR", loc="left", weight="bold")
    ax1.grid(axis="x", color=COLORS["grid"])

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.axis("off")
    ax2.set_title("B  eQTL and locus ranks", loc="left", weight="bold")
    table_rows = [
        ["Sun skin", "rs35251378", "1 / 4"],
        ["Non-sun skin", "rs280497", "1 / 11"],
        ["Blood", "rs280497", "1 / 11"],
    ]
    table = ax2.table(
        cellText=table_rows,
        colLabels=["Tissue", "MR SNP", "Ranks"],
        loc="center",
        cellLoc="center",
        colLoc="center",
        colWidths=[0.33, 0.34, 0.33],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 1.8)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor(COLORS["grid"])
        if row == 0:
            cell.set_facecolor(COLORS["ink"])
            cell.set_text_props(color="white", weight="bold")
        else:
            cell.set_facecolor("#F7F8F9" if row % 2 else "white")
    ax2.text(
        0.5,
        0.20,
        "Ranks: eQTL / FinnGen locus",
        ha="center",
        transform=ax2.transAxes,
        color=COLORS["muted"],
    )
    ax2.text(
        0.5,
        0.10,
        "FinnGen locus lead: rs34536443",
        ha="center",
        transform=ax2.transAxes,
        color=COLORS["muted"],
    )

    ax3 = fig.add_subplot(gs[1, :])
    hypotheses = ["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]
    colors = ["#B8BEC5", "#6B9AC4", COLORS["exploratory"], "#C76D4A", COLORS["green"]]
    left = np.zeros(3)
    for hypothesis, color in zip(hypotheses, colors):
        values = coloc[hypothesis.replace(".", "_")].astype(float).to_numpy()
        ax3.barh(y, values, left=left, color=color, edgecolor="white", height=0.58, label=hypothesis)
        left += values
    ax3.set_yticks(y, [TISSUE_SHORT[v] for v in tissue_order])
    ax3.invert_yaxis()
    ax3.set_xlim(0, 1)
    ax3.set_xlabel("Posterior probability")
    ax3.set_title("C  Colocalization", loc="left", weight="bold")
    for i, value in enumerate(coloc["PP_H4"].astype(float)):
        ax3.text(1.01, i, f"H4={value:.3f}", va="center", fontsize=10.5, color=COLORS["ink"])
    ax3.legend(loc="lower center", bbox_to_anchor=(0.5, -0.28), ncol=5, frameon=False)

    fig.suptitle("TYK2 positive-control validation", x=0.02, y=0.98, ha="left", fontsize=18, weight="bold")
    fig.subplots_adjust(left=0.10, right=0.94, top=0.90, bottom=0.13)
    save(fig, "figure4_tyk2_validation_summary")


def build_supplementary_figure(opengwas: pd.DataFrame) -> None:
    tissue_order = ["Skin_Sun_Exposed_Lower_leg", "Skin_Not_Sun_Exposed_Suprapubic", "Whole_Blood"]
    tissue_colors = {
        "Skin_Sun_Exposed_Lower_leg": COLORS["protective"],
        "Skin_Not_Sun_Exposed_Suprapubic": COLORS["exploratory"],
        "Whole_Blood": COLORS["green"],
    }
    opengwas = opengwas.copy()
    opengwas["tissue_key"] = opengwas["tissue"].str.replace("_Suprapubic", "", regex=False).str.replace("_Lower_leg", "", regex=False)
    # The source uses full labels for two tissues and Whole_Blood.
    normalized = {
        "Skin_Not_Sun_Exposed": "Skin_Not_Sun_Exposed_Suprapubic",
        "Skin_Sun_Exposed": "Skin_Sun_Exposed_Lower_leg",
        "Whole_Blood": "Whole_Blood",
    }
    opengwas["tissue_full"] = opengwas["tissue_key"].map(normalized).fillna(opengwas["tissue"])
    opengwas["tissue_order"] = opengwas["tissue_full"].map({v: i for i, v in enumerate(tissue_order)})
    opengwas = opengwas.sort_values(["outcome_id", "tissue_order"]).reset_index(drop=True)
    labels = [
        f"{row.outcome_id} | {TISSUE_SHORT.get(row.tissue_full, row.tissue)}"
        for row in opengwas.itertuples()
    ]
    y = np.arange(len(opengwas))

    fig, ax = plt.subplots(figsize=(11.5, 12.5))
    for tissue in tissue_order:
        subset = opengwas[opengwas["tissue_full"] == tissue]
        yy = subset.index.to_numpy()
        ax.errorbar(
            subset["OR"].astype(float),
            yy,
            xerr=[
                subset["OR"].astype(float) - subset["OR_lci"].astype(float),
                subset["OR_uci"].astype(float) - subset["OR"].astype(float),
            ],
            fmt="o",
            color=tissue_colors[tissue],
            ecolor=tissue_colors[tissue],
            capsize=2,
            markersize=5,
            linewidth=1.3,
            label=TISSUE_SHORT[tissue],
        )
    ax.axvline(1, color=COLORS["muted"], linestyle="--")
    ax.set_xscale("log")
    ax.set_xlim(0.15, 1.15)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("External-outcome Wald-ratio OR (95% CI)")
    ax.grid(axis="x", color=COLORS["grid"])
    ax.legend(loc="upper right", framealpha=1.0)
    fig.subplots_adjust(left=0.39, right=0.97, top=0.98, bottom=0.07)
    save(fig, "figureS1_opengwas_forest_plot")


def main() -> None:
    setup_style()
    primary = pd.read_csv(PRIMARY_PATH, keep_default_na=False)
    availability = pd.read_csv(AVAILABILITY_PATH)
    opengwas = pd.read_csv(OPENGWAS_PATH)
    opengwas = opengwas[
        (opengwas["status"] == "available")
        & (opengwas["association_type"] == "exact")
    ]
    build_figure1()
    build_figure2(primary)
    build_figure3(availability)
    build_figure4(primary)
    build_supplementary_figure(opengwas)
    print(f"Wrote revised figures to {OUT}")


if __name__ == "__main__":
    main()
