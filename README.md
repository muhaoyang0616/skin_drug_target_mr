# Genetic Tractability of Dermatologic Immune Pathway Genes

This repository contains the validated analysis code, machine-readable results,
figures, and manuscript sources for:

> Genetic Tractability of Clinically Relevant Dermatologic Immune Pathway
> Genes: A Tissue-Aware cis-eQTL Mendelian Randomization Screen

## Validated Analysis

- A curated panel of 19 genes was evaluated in three GTEx v8 tissues.
- Primary instruments required cis-eQTL `P < 5e-8` and `F > 10`.
- One deterministic lead variant was retained per gene-tissue cell because no
  matched LD reference panel was available.
- Eighteen of 57 cells had a primary instrument. One PDE4A whole-blood variant
  was absent from all four FinnGen outcomes, leaving 68 primary tests.
- Three primary tests were FDR-significant, all linking TYK2 expression proxies
  to psoriasis. These are treated as recovery of an established locus, not a
  new TYK2 target discovery or a direct proxy for drug inhibition.
- The five relaxed-only cells generated 20 separately corrected exploratory
  tests. The PDE4A-psoriasis signal lacked locus-level colocalization support.

The raw-data audit independently reread seven compressed GTEx and FinnGen files.
Python and R agreed for 23 instruments, 68 primary estimates, and 20 exploratory
estimates. MR beta, standard-error, p-value, and FDR differences were below
`5e-15`; the maximum F-statistic recalculation difference was below `6e-13`.

## Repository Layout

```text
analysis_config.yaml
data/
  target_genes/
  processed/
  raw/finngen/outcome_files.csv
manuscript/
  main.tex
  supplementary.tex
  sections/
  figures/
results/
scripts/
validation/
```

The versioned public release is available at:
`https://github.com/muhaoyang0616/skin_drug_target_mr/releases/tag/v1.0.0`.

Raw GTEx and FinnGen files are not redistributed. Their expected relative paths,
source URLs, byte sizes, and SHA-256 values are documented in
`SOURCE_PROVENANCE.md`.

## Environment

Validated software versions:

- R 4.6.0
- TwoSampleMR 0.7.6
- coloc 5.2.3
- data.table 1.18.4
- dplyr 1.2.1
- yaml 2.3.12
- readr 2.2.0
- Python 3.10.4
- pandas 2.3.3
- NumPy 2.2.6
- Matplotlib 3.10.9
- requests 2.33.1

See `R_PACKAGES.md` and `requirements.txt`.
`requests` is used only by the optional eQTL Catalogue HTTP fallback helper.

## Required Raw Inputs

Place the GTEx files under:

```text
data/raw/gtex/GTEx_Analysis_v8_eQTL/
```

Place the four FinnGen R13 files under:

```text
data/raw/finngen/outcomes/
```

The exact filenames are listed in `data/raw/README.md` and
`data/raw/finngen/outcome_files.csv`.

## Reproduce the Primary R Analysis

From the repository root:

```bash
Rscript scripts/00_setup.R
Rscript scripts/02_prepare_gtex_instruments.R
Rscript scripts/04_run_mr_gtex_finngen.R
```

`04_run_mr_gtex_finngen.R` scans the large FinnGen files in chunks and retains
only positions required by the selected instruments.

## Independent Python Audit

```bash
python scripts/reproduce_primary_mr.py --source . --output validation/python_reanalysis
```

The audit rereads the downloaded raw inputs and compares the regenerated primary
estimates with the published `results/table_s3_primary_mr_68.csv`; it does not
depend on an R-generated temporary result file.

The expected row counts are 18 primary instrument cells, five relaxed-only
cells, 68 primary tests, and 20 exploratory tests.

## Optional Validation Refresh

Run the top-hit and colocalization workflows after downloading the documented
raw inputs:

```bash
python scripts/10_top_hit_validation.py
Rscript scripts/11_coloc_eqtlcatalogue.R
Rscript scripts/12_coloc_pde4a_jak2.R
```

The last script retains its historical filename for compatibility but evaluates
only the exploratory PDE4A locus in the final analysis. An OpenGWAS refresh is
optional and additionally requires `ieugwasr` and an `OPENGWAS_JWT`:

```bash
Rscript scripts/10_opengwas_psoriasis_replication.R
```

## Regenerate Figures

```bash
python scripts/make_figures.py
```

Figures are written to `manuscript/figures/` in PDF and 300-dpi PNG formats.
PDF fonts are embedded Arial subsets.

## Manuscript Build

The manuscript source uses standard LaTeX plus `natbib`, `booktabs`, and
`hyperref`. Run BibTeX between LaTeX passes for `main.tex`. The supplementary
file is self-contained apart from `manuscript/figures/`.

## Interpretation Boundaries

This release does not include or support the invalidated legacy 112-row
multi-SNP analysis, reverse MR, Steiger directionality, formal power percentages,
or pQTL-null claims. Absence of an eligible instrument is reported as
untestability in the current data, not evidence against a therapeutic target.

## License

Source code is released under the MIT License; see `LICENSE`. Manuscript source,
figures, documentation, and machine-readable result tables are released under
the Creative Commons Attribution 4.0 International License; see
`LICENSE-DATA.md`. Third-party source data are not redistributed and remain
subject to the terms of their originating resources.
