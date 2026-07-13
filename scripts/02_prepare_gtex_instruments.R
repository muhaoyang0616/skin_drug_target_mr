args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)) else "scripts"
source(file.path(script_dir, "utils.R"))

cfg <- read_config()
ensure_dir(cfg$project$processed_dir)
ensure_dir(cfg$project$results_dir)
ensure_dir(cfg$data$gtex$raw_dir)

annotation <- read_csv_required(cfg$data$annotations$target_gene_annotation, "target gene GENCODE annotation")
annotation <- annotation[!is.na(annotation$gene_id) & nzchar(annotation$gene_id), , drop = FALSE]

instrument_cols <- c(
  "SNP", "variant_id", "variant_key", "target_gene", "gene_id", "tissue", "tissue_label",
  "exposure", "id.exposure", "effect_allele.exposure", "other_allele.exposure",
  "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure",
  "samplesize.exposure", "chr.exposure", "pos.exposure", "f_stat",
  "instrument_tier", "clump_status", "source"
)

parse_variant_id <- function(x) {
  parts <- strsplit(as.character(x), "_", fixed = TRUE)
  data.frame(
    chr = vapply(parts, function(p) if (length(p) >= 1) p[[1]] else NA_character_, character(1)),
    pos = vapply(parts, function(p) if (length(p) >= 2) p[[2]] else NA_character_, character(1)),
    ref = vapply(parts, function(p) if (length(p) >= 3) p[[3]] else NA_character_, character(1)),
    alt = vapply(parts, function(p) if (length(p) >= 4) p[[4]] else NA_character_, character(1)),
    stringsAsFactors = FALSE
  )
}

ensure_gtex_tissue_file <- function(path, tissue_name) {
  if (file.exists(path) && file.info(path)$size > 1024 * 1024) {
    return(invisible(path))
  }
  tar_file <- cfg$data$gtex$tar_file
  if (!file.exists(tar_file)) {
    stop(
      "Missing GTEx tissue file for ", tissue_name, ": ", path, "\n",
      "GTEx v8 single tissue URLs are not used. Download the full tar archive to: ", tar_file,
      call. = FALSE
    )
  }
  status_message("Extracting ", basename(path), " from GTEx tar archive.")
  utils::untar(tar_file, files = basename(path), exdir = cfg$data$gtex$raw_dir)
  if (!file.exists(path) || file.info(path)$size <= 1024 * 1024) {
    stop("Could not extract required GTEx tissue file from tar archive: ", basename(path), call. = FALSE)
  }
  invisible(path)
}

prepare_one_tissue <- function(tissue_name, tissue_cfg) {
  path <- tissue_cfg$file
  ensure_gtex_tissue_file(path, tissue_name)
  status_message("Reading GTEx tissue file: ", tissue_name)
  x <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  cols <- names(x)
  variant_col <- find_col(cols, c("variant_id", "variant"), "GTEx variant_id")
  gene_col <- find_col(cols, c("gene_id", "gene"), "GTEx gene_id")
  p_col <- find_col(cols, c("pval_nominal", "p", "pval"), "GTEx nominal p-value")
  beta_col <- find_col(cols, c("slope", "beta", "effect"), "GTEx slope/beta")
  se_col <- find_col(cols, c("slope_se", "se", "stderr"), "GTEx slope SE")
  maf_col <- find_col(cols, c("maf", "eaf", "af"), "GTEx MAF", required = FALSE)
  ma_samples_col <- find_col(cols, c("ma_samples", "samplesize", "n"), "GTEx minor allele samples", required = FALSE)

  x$gene_id_stripped <- strip_ensembl_version(x[[gene_col]])
  x <- merge(
    x,
    annotation[, c("gene", "gene_id"), drop = FALSE],
    by.x = "gene_id_stripped",
    by.y = "gene_id",
    all = FALSE
  )
  if (nrow(x) == 0) {
    return(empty_df(instrument_cols))
  }

  beta <- safe_numeric(x[[beta_col]])
  se <- safe_numeric(x[[se_col]])
  pval <- standardize_p(x[[p_col]])
  f_stat <- (beta / se)^2
  keep_sensitivity <- !is.na(pval) & !is.na(f_stat) & pval < cfg$instruments$pval_sensitivity_threshold & f_stat > cfg$instruments$f_stat_min
  x <- x[keep_sensitivity, , drop = FALSE]
  beta <- beta[keep_sensitivity]
  se <- se[keep_sensitivity]
  pval <- pval[keep_sensitivity]
  f_stat <- f_stat[keep_sensitivity]
  if (nrow(x) == 0) {
    return(empty_df(instrument_cols))
  }

  parsed <- parse_variant_id(x[[variant_col]])
  parsed$chr <- normalize_chr(parsed$chr)
  parsed$pos <- safe_numeric(parsed$pos)
  tier <- ifelse(
    pval < cfg$instruments$pval_threshold & f_stat > cfg$instruments$f_stat_min,
    "primary",
    "exploratory_relaxed"
  )
  out <- data.frame(
    SNP = as.character(x[[variant_col]]),
    variant_id = as.character(x[[variant_col]]),
    variant_key = make_variant_key(parsed$chr, parsed$pos, parsed$ref, parsed$alt),
    target_gene = toupper(x$gene),
    gene_id = x$gene_id_stripped,
    tissue = tissue_name,
    tissue_label = tissue_cfg$label,
    exposure = paste(toupper(x$gene), tissue_name, sep = "_"),
    id.exposure = paste(toupper(x$gene), tissue_name, sep = "_"),
    effect_allele.exposure = toupper(parsed$alt),
    other_allele.exposure = toupper(parsed$ref),
    eaf.exposure = if (!is.na(maf_col)) safe_numeric(x[[maf_col]]) else NA_real_,
    beta.exposure = beta,
    se.exposure = se,
    pval.exposure = pval,
    samplesize.exposure = if (!is.na(ma_samples_col)) safe_numeric(x[[ma_samples_col]]) else NA_real_,
    chr.exposure = parsed$chr,
    pos.exposure = parsed$pos,
    f_stat = f_stat,
    instrument_tier = tier,
    clump_status = "lead_variant_only_no_ld_reference",
    source = "GTEx_v8_significant_variant_gene_pairs",
    stringsAsFactors = FALSE
  )
  out <- out[complete.cases(out[, c("variant_key", "target_gene", "beta.exposure", "se.exposure", "pval.exposure")]), , drop = FALSE]
  if (nrow(out) == 0) return(empty_df(instrument_cols))

  do.call(rbind, lapply(split(out, list(out$target_gene, out$tissue), drop = TRUE), function(df) {
    df <- df[order(df$pval.exposure, df$variant_id), , drop = FALSE]
    primary <- df[df$instrument_tier == "primary", , drop = FALSE]
    if (nrow(primary) > 0) primary[1, , drop = FALSE] else df[1, , drop = FALSE]
  }))
}

tissues <- cfg$data$gtex$tissues
all_inst <- list()
for (tissue_name in names(tissues)) {
  all_inst[[tissue_name]] <- prepare_one_tissue(tissue_name, tissues[[tissue_name]])
}

instruments <- data.table::rbindlist(all_inst, fill = TRUE)
if (nrow(instruments) == 0) {
  instruments <- empty_df(instrument_cols)
} else {
  instruments <- as.data.frame(instruments)
  instruments <- instruments[, instrument_cols, drop = FALSE]
  instruments <- unique(instruments)
}

write_csv(instruments, file.path(cfg$project$processed_dir, "instruments_gtex.csv"))

targets <- read_csv_required(cfg$targets$target_gene_file, "target gene file")
gene_counts <- expand.grid(
  gene = unique(toupper(targets$gene)),
  instrument_tier = c("primary", "exploratory_relaxed"),
  stringsAsFactors = FALSE
)
if (nrow(instruments) > 0) {
  obs <- aggregate(variant_key ~ target_gene + instrument_tier, instruments, length)
  names(obs) <- c("gene", "instrument_tier", "n_instruments")
  gene_counts <- merge(gene_counts, obs, by = c("gene", "instrument_tier"), all.x = TRUE)
} else {
  gene_counts$n_instruments <- NA_integer_
}
gene_counts$n_instruments[is.na(gene_counts$n_instruments)] <- 0L
write_csv(gene_counts, file.path(cfg$project$results_dir, "instrument_counts_by_gene_gtex.csv"))

tissue_counts <- expand.grid(
  tissue = names(tissues),
  instrument_tier = c("primary", "exploratory_relaxed"),
  stringsAsFactors = FALSE
)
if (nrow(instruments) > 0) {
  obs_tissue <- aggregate(target_gene ~ tissue + instrument_tier, instruments, function(x) length(unique(x)))
  names(obs_tissue) <- c("tissue", "instrument_tier", "n_target_genes_with_instruments")
  tissue_counts <- merge(tissue_counts, obs_tissue, by = c("tissue", "instrument_tier"), all.x = TRUE)
} else {
  tissue_counts$n_target_genes_with_instruments <- NA_integer_
}
tissue_counts$n_target_genes_with_instruments[is.na(tissue_counts$n_target_genes_with_instruments)] <- 0L
write_csv(tissue_counts, file.path(cfg$project$results_dir, "instrument_counts_by_tissue_gtex.csv"))

append_log("logs/run_log.txt", "GTEx instruments prepared: ", nrow(instruments), " rows; clump_status=lead_variant_only_no_ld_reference")
status_message("Wrote GTEx instruments: ", nrow(instruments), " rows")
