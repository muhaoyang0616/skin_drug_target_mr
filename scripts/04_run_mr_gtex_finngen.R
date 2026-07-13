args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)) else "scripts"
source(file.path(script_dir, "utils.R"))

cfg <- read_config()
inst_path <- file.path(cfg$project$processed_dir, "instruments_gtex.csv")
outcome_index_path <- cfg$data$finngen$outcome_index_path
ensure_dir(cfg$project$results_dir)

mr_cols <- c(
  "target_gene", "tissue", "tissue_label", "disease_query", "endpoint", "method",
  "instrument_tier", "nsnp", "b", "se", "pval", "fdr", "method_is_primary",
  "n_outcome_variants", "n_harmonised", "harmonise_note"
)

write_empty_mr <- function(reason) {
  empty <- empty_df(mr_cols)
  write_csv(empty, file.path(cfg$project$results_dir, "mr_all_results.csv"))
  write_csv(empty, file.path(cfg$project$results_dir, "all_mr_results.csv"))
  write_csv(empty, file.path(cfg$project$results_dir, "mr_significant_results.csv"))
  write_csv(empty, file.path(cfg$project$results_dir, "mr_suggestive_results.csv"))
  write_csv(empty, file.path(cfg$project$results_dir, "mr_top_hits.csv"))
  write_csv(empty, file.path(cfg$project$results_dir, "mr_exploratory_relaxed_results.csv"))
  write_csv(empty, file.path(cfg$project$results_dir, "sensitivity_results.csv"))
  write_unavailable(file.path(cfg$project$results_dir, "mr_unavailable.csv"), "GTEx-FinnGen MR", reason)
  append_log("logs/run_log.txt", "GTEx-FinnGen MR unavailable: ", reason)
  status_message("GTEx-FinnGen MR unavailable: ", reason)
}

if (!file.exists(inst_path)) {
  write_empty_mr("Missing data/processed/instruments_gtex.csv. Run scripts/02_prepare_gtex_instruments.R first.")
  quit(status = 0)
}
if (!file.exists(outcome_index_path)) {
  write_empty_mr("Missing FinnGen outcome index. Run scripts/01_find_finngen_endpoints.R first.")
  quit(status = 0)
}

instruments <- read_csv_required(inst_path, "GTEx instruments")
outcome_index <- read_csv_required(outcome_index_path, "FinnGen outcome file index")
if (nrow(instruments) == 0) {
  write_empty_mr("No GTEx instruments are available for the target gene panel.")
  quit(status = 0)
}
if (nrow(outcome_index) == 0) {
  write_empty_mr("No FinnGen outcome files were downloaded for selected endpoints.")
  quit(status = 0)
}

format_outcome <- function(path, endpoint, disease_query) {
  check_file(path, paste0("FinnGen outcome ", endpoint))
  header <- data.table::fread(path, nrows = 0, data.table = FALSE, showProgress = FALSE)
  cols <- names(header)
  chr_col <- find_col(cols, c("chr", "chrom", "chromosome", "#chrom", "x.chrom"), "outcome chromosome")
  pos_col <- find_col(cols, c("pos", "position", "bp"), "outcome position")
  ref_col <- find_col(cols, c("ref", "reference_allele", "a2", "nea"), "outcome reference allele")
  alt_col <- find_col(cols, c("alt", "effect_allele", "a1", "ea"), "outcome alternate/effect allele")
  beta_col <- find_col(cols, c("beta", "b", "log_or", "effect"), "outcome beta/logOR")
  se_col <- find_col(cols, c("sebeta", "se", "stderr", "standard_error"), "outcome SE")
  p_col <- find_col(cols, c("pval", "p", "pvalue", "p_value"), "outcome p-value")
  eaf_col <- find_col(cols, c("af_alt", "eaf", "effect_allele_frequency", "af"), "outcome alt/effect allele frequency", required = FALSE)
  rsid_col <- find_col(cols, c("rsids", "rsid", "snp"), "outcome rsid", required = FALSE)
  select_cols <- unique(na.omit(c(chr_col, pos_col, ref_col, alt_col, beta_col, se_col, p_col, eaf_col, rsid_col)))
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required for memory-safe FinnGen chunk scanning.", call. = FALSE)
  }
  target_positions <- unique(paste(
    normalize_chr(instruments$chr.exposure),
    safe_numeric(instruments$pos.exposure),
    sep = ":"
  ))
  scanned_rows <- 0L
  callback <- readr::DataFrameCallback$new(function(chunk, pos) {
    scanned_rows <<- scanned_rows + nrow(chunk)
    chunk_chr <- normalize_chr(chunk[[chr_col]])
    chunk_pos <- safe_numeric(chunk[[pos_col]])
    keep <- paste(chunk_chr, chunk_pos, sep = ":") %in% target_positions
    as.data.frame(chunk[keep, select_cols, drop = FALSE])
  })
  df <- readr::read_tsv_chunked(
    path,
    callback = callback,
    chunk_size = 250000,
    col_types = readr::cols(.default = readr::col_character()),
    progress = FALSE,
    show_col_types = FALSE
  )
  chr <- normalize_chr(df[[chr_col]])
  pos <- safe_numeric(df[[pos_col]])
  ref <- toupper(as.character(df[[ref_col]]))
  alt <- toupper(as.character(df[[alt_col]]))
  out <- data.frame(
    variant_key = make_variant_key(chr, pos, ref, alt),
    reverse_variant_key = make_variant_key(chr, pos, alt, ref),
    chr.outcome = chr,
    pos.outcome = pos,
    ref.outcome = ref,
    alt.outcome = alt,
    beta.outcome = safe_numeric(df[[beta_col]]),
    se.outcome = safe_numeric(df[[se_col]]),
    pval.outcome = standardize_p(df[[p_col]]),
    eaf.outcome = if (!is.na(eaf_col)) safe_numeric(df[[eaf_col]]) else NA_real_,
    rsid.outcome = if (!is.na(rsid_col)) as.character(df[[rsid_col]]) else NA_character_,
    outcome = endpoint,
    endpoint = endpoint,
    disease_query = disease_query,
    stringsAsFactors = FALSE
  )
  attr(out, "n_rows_scanned") <- scanned_rows
  out
}

harmonise_by_variant_key <- function(exp_dat, outcome) {
  forward <- merge(exp_dat, outcome, by = "variant_key", all = FALSE)
  if (nrow(forward) > 0) {
    forward$harmonise_note <- "matched_chr_pos_ref_alt; FinnGen beta assumed alt allele"
  }
  rev_out <- outcome
  rev_out$variant_key <- rev_out$reverse_variant_key
  rev_out$beta.outcome <- -rev_out$beta.outcome
  if (nrow(rev_out) > 0) {
    rev_out$harmonise_note <- "matched_ref_alt_reversed; FinnGen beta flipped"
  }
  reversed <- merge(exp_dat, rev_out, by = "variant_key", all = FALSE)
  out <- as.data.frame(data.table::rbindlist(list(forward, reversed), fill = TRUE))
  out <- out[complete.cases(out[, c("beta.exposure", "se.exposure", "beta.outcome", "se.outcome")]), , drop = FALSE]
  out <- out[is.finite(out$beta.exposure) & out$beta.exposure != 0 & is.finite(out$se.outcome) & out$se.outcome > 0, , drop = FALSE]
  out <- out[!duplicated(out$variant_key), , drop = FALSE]
  out
}

wald_ratio <- function(h) {
  b <- h$beta.outcome / h$beta.exposure
  se <- abs(h$se.outcome / h$beta.exposure)
  p <- 2 * stats::pnorm(abs(b / se), lower.tail = FALSE)
  c(b = b, se = se, pval = p)
}

ivw_fixed <- function(h) {
  theta <- h$beta.outcome / h$beta.exposure
  theta_se <- sqrt((h$se.outcome^2 / h$beta.exposure^2) + ((h$beta.outcome^2 * h$se.exposure^2) / h$beta.exposure^4))
  w <- 1 / theta_se^2
  keep <- is.finite(theta) & is.finite(w) & w > 0
  theta <- theta[keep]
  w <- w[keep]
  b <- sum(w * theta) / sum(w)
  se <- sqrt(1 / sum(w))
  p <- 2 * stats::pnorm(abs(b / se), lower.tail = FALSE)
  c(b = b, se = se, pval = p)
}

weighted_median <- function(h) {
  theta <- h$beta.outcome / h$beta.exposure
  theta_se <- sqrt((h$se.outcome^2 / h$beta.exposure^2) + ((h$beta.outcome^2 * h$se.exposure^2) / h$beta.exposure^4))
  w <- 1 / theta_se^2
  keep <- is.finite(theta) & is.finite(w) & w > 0
  theta <- theta[keep]
  w <- w[keep]
  ord <- order(theta)
  theta <- theta[ord]
  w <- w[ord] / sum(w[ord])
  b <- theta[which(cumsum(w) >= 0.5)[[1]]]
  c(b = b, se = NA_real_, pval = NA_real_)
}

egger <- function(h) {
  fit <- tryCatch(stats::lm(beta.outcome ~ beta.exposure, data = h, weights = 1 / (se.outcome^2)), error = function(e) NULL)
  if (is.null(fit)) return(c(b = NA_real_, se = NA_real_, pval = NA_real_))
  sm <- summary(fit)$coefficients
  if (!("beta.exposure" %in% rownames(sm))) return(c(b = NA_real_, se = NA_real_, pval = NA_real_))
  c(b = sm["beta.exposure", "Estimate"], se = sm["beta.exposure", "Std. Error"], pval = sm["beta.exposure", "Pr(>|t|)"])
}

run_methods <- function(h, tier, n_outcome_variants) {
  n <- nrow(h)
  rows <- list()
  if (n == 1) {
    est <- wald_ratio(h)
    rows[["Wald ratio"]] <- c(est, method_is_primary = TRUE)
  } else if (n >= 2) {
    est <- ivw_fixed(h)
    rows[["Inverse variance weighted"]] <- c(est, method_is_primary = TRUE)
    if (n >= 3) {
      rows[["Weighted median"]] <- c(weighted_median(h), method_is_primary = FALSE)
      rows[["MR-Egger"]] <- c(egger(h), method_is_primary = FALSE)
    }
  }
  if (length(rows) == 0) return(empty_df(mr_cols))
  do.call(rbind, lapply(names(rows), function(method) {
    est <- rows[[method]]
    data.frame(
      target_gene = h$target_gene[[1]],
      tissue = h$tissue[[1]],
      tissue_label = h$tissue_label[[1]],
      disease_query = h$disease_query[[1]],
      endpoint = h$endpoint[[1]],
      method = method,
      instrument_tier = tier,
      nsnp = n,
      b = as.numeric(est[["b"]]),
      se = as.numeric(est[["se"]]),
      pval = as.numeric(est[["pval"]]),
      fdr = NA_real_,
      method_is_primary = as.logical(est[["method_is_primary"]]),
      n_outcome_variants = n_outcome_variants,
      n_harmonised = n,
      harmonise_note = paste(unique(h$harmonise_note), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))
}

all_results <- list()
exploratory_results <- list()

for (i in seq_len(nrow(outcome_index))) {
  path <- as.character(outcome_index$file[[i]])
  endpoint <- as.character(outcome_index$endpoint[[i]])
  disease_query <- as.character(outcome_index$disease_query[[i]])
  status_message("Reading FinnGen outcome for MR: ", endpoint)
  outcome <- tryCatch(format_outcome(path, endpoint, disease_query), error = function(e) {
    append_log("logs/run_log.txt", "Could not format FinnGen outcome ", endpoint, ": ", conditionMessage(e))
    NULL
  })
  if (is.null(outcome) || nrow(outcome) == 0) next

  for (tier in c("primary", "exploratory_relaxed")) {
    exp_tier <- instruments[instruments$instrument_tier == tier, , drop = FALSE]
    if (nrow(exp_tier) == 0) next
    groups <- split(exp_tier, list(exp_tier$target_gene, exp_tier$tissue), drop = TRUE)
    for (group_name in names(groups)) {
      h <- harmonise_by_variant_key(groups[[group_name]], outcome)
      if (nrow(h) == 0) next
      res <- run_methods(h, tier, attr(outcome, "n_rows_scanned"))
      if (nrow(res) == 0) next
      if (tier == "primary") {
        all_results[[length(all_results) + 1]] <- res
      } else {
        exploratory_results[[length(exploratory_results) + 1]] <- res
      }
    }
  }
}

mr_out <- if (length(all_results) > 0) data.table::rbindlist(all_results, fill = TRUE) else empty_df(mr_cols)
mr_out <- as.data.frame(mr_out)
if (nrow(mr_out) > 0) {
  primary_idx <- which(mr_out$method_is_primary & is.finite(mr_out$pval))
  mr_out$fdr[primary_idx] <- stats::p.adjust(mr_out$pval[primary_idx], method = cfg$mr$fdr_method)
}

exploratory_out <- if (length(exploratory_results) > 0) data.table::rbindlist(exploratory_results, fill = TRUE) else empty_df(mr_cols)
exploratory_out <- as.data.frame(exploratory_out)
if (nrow(exploratory_out) > 0) {
  exploratory_idx <- which(exploratory_out$method_is_primary & is.finite(exploratory_out$pval))
  exploratory_out$fdr[exploratory_idx] <- stats::p.adjust(
    exploratory_out$pval[exploratory_idx],
    method = cfg$mr$fdr_method
  )
}

sig <- if (nrow(mr_out) > 0) {
  mr_out[mr_out$method_is_primary & is.finite(mr_out$fdr) & mr_out$fdr <= cfg$mr$top_hit_fdr_threshold, , drop = FALSE]
} else {
  mr_out
}
suggestive <- if (nrow(mr_out) > 0) {
  mr_out[
    mr_out$method_is_primary & is.finite(mr_out$pval) &
      mr_out$pval <= cfg$mr$nominal_p_threshold &
      (is.na(mr_out$fdr) | mr_out$fdr > cfg$mr$top_hit_fdr_threshold),
    , drop = FALSE
  ]
} else {
  mr_out
}

write_csv(mr_out, file.path(cfg$project$results_dir, "mr_all_results.csv"))
write_csv(mr_out, file.path(cfg$project$results_dir, "all_mr_results.csv"))
write_csv(sig, file.path(cfg$project$results_dir, "mr_significant_results.csv"))
write_csv(suggestive, file.path(cfg$project$results_dir, "mr_suggestive_results.csv"))
write_csv(sig, file.path(cfg$project$results_dir, "mr_top_hits.csv"))
write_csv(exploratory_out, file.path(cfg$project$results_dir, "mr_exploratory_relaxed_results.csv"))
write_csv(exploratory_out, file.path(cfg$project$results_dir, "sensitivity_results.csv"))

if (nrow(mr_out) == 0) {
  write_unavailable(
    file.path(cfg$project$results_dir, "mr_unavailable.csv"),
    "GTEx-FinnGen MR",
    "No harmonised GTEx lead instruments overlapped downloaded FinnGen outcomes."
  )
} else if (file.exists(file.path(cfg$project$results_dir, "mr_unavailable.csv"))) {
  unlink(file.path(cfg$project$results_dir, "mr_unavailable.csv"))
}

append_log("logs/run_log.txt", "GTEx-FinnGen MR complete: ", nrow(mr_out), " primary rows; ", nrow(sig), " FDR-significant rows.")
status_message("GTEx-FinnGen MR complete. Rows: ", nrow(mr_out), "; top hits: ", nrow(sig), "; exploratory relaxed rows: ", nrow(exploratory_out))
