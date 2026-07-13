args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  "scripts"
}
source(file.path(script_dir, "utils.R"))

required <- c("jsonlite", "data.table", "coloc", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

cfg <- read_config()
validation_dir <- file.path(cfg$project$results_dir, "validation")
fig_validation_dir <- file.path(cfg$project$figures_dir, "validation")
log_dir <- file.path("logs", "validation")
ensure_dir(validation_dir)
ensure_dir(fig_validation_dir)
ensure_dir(log_dir)

status_path <- file.path(validation_dir, "coloc_eqtlcatalogue_status.md")
summary_path <- file.path(validation_dir, "coloc_tyk2_psoriasis_summary.csv")
run_log_path <- file.path(log_dir, "coloc_eqtlcatalogue_api.log")

api_base_v2 <- "https://www.ebi.ac.uk/eqtl/api/v2"
api_base_v3 <- "https://www.ebi.ac.uk/eqtl/api/v3"
tyk2_gene_id <- "ENSG00000105397"
chrom <- "19"
locus_start <- 9835400L
locus_end <- 10854000L
page_size <- 1000L
chunk_size <- 100000L
min_overlap <- 50L

targets <- data.frame(
  tissue = c(
    "skin_sun_exposed_lower_leg_skin",
    "skin_not_sun_exposed_suprapubic",
    "whole_blood"
  ),
  requested_dataset_id = c(
    "GTEx_V8.skin_sun_exposed_lower_leg_skin",
    "GTEx_V8.skin_not_sun_exposed_suprapubic",
    "GTEx_V8.whole_blood"
  ),
  dataset_id_hint = c(
    "QTD000316",
    "QTD000311",
    "QTD000356"
  ),
  sample_group_pattern = c(
    "skin_sun_exposed",
    "skin_not_sun_exposed",
    "^blood$"
  ),
  tissue_pattern = c(
    "skin",
    "skin",
    "blood"
  ),
  stringsAsFactors = FALSE
)

status_lines <- c(
  "# eQTL Catalogue TYK2 coloc status",
  "",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Principle: no coloc result is fabricated. REST API data are used when available.",
  ""
)

add_status <- function(...) {
  status_lines <<- c(status_lines, paste0(...))
}

write_status <- function() {
  writeLines(status_lines, status_path, useBytes = TRUE)
}

fail_run <- function(reason) {
  add_status("Status: **unavailable**")
  add_status("")
  add_status("Reason: ", reason)
  write_status()
  stop(reason, call. = FALSE)
}

if (length(missing) > 0) {
  fail_run(paste0("Missing required R packages: ", paste(missing, collapse = ", ")))
}

append_log(run_log_path, "Starting eQTL Catalogue REST coloc run")

api_url <- function(base, path, query = list()) {
  query <- query[!vapply(query, is.null, logical(1))]
  if (length(query) == 0) {
    return(paste0(base, path))
  }
  vals <- vapply(query, function(x) utils::URLencode(as.character(x), reserved = TRUE), character(1))
  paste0(base, path, "?", paste(paste(names(vals), vals, sep = "="), collapse = "&"))
}

api_get <- function(url, label, timeout_sec = 120, retries = 3) {
  append_log(run_log_path, label, " | GET ", url)
  old_timeout <- getOption("timeout")
  options(timeout = max(timeout_sec, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  last_error <- NA_character_
  for (attempt in seq_len(retries)) {
    res <- tryCatch(
      jsonlite::fromJSON(url, simplifyDataFrame = TRUE),
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(res)) {
      return(res)
    }
    append_log(run_log_path, label, " | ERROR attempt ", attempt, "/", retries, " | ", last_error)
    if (attempt < retries) {
      Sys.sleep(min(10, attempt * 2))
    }
  }
  structure(list(error = last_error), class = "api_error")
}

is_api_error <- function(x) inherits(x, "api_error")

as_records <- function(x) {
  if (is_api_error(x) || is.null(x)) {
    return(data.frame())
  }
  if (is.data.frame(x)) {
    return(x)
  }
  if (is.list(x) && "value" %in% names(x) && is.data.frame(x$value)) {
    return(x$value)
  }
  data.frame()
}

datasets_url <- api_url(api_base_v2, "/datasets", list(study_label = "GTEx", start = 0, size = 1000))
datasets_raw <- api_get(datasets_url, "datasets")
datasets <- as_records(datasets_raw)
if (is_api_error(datasets_raw) || nrow(datasets) == 0) {
  fail_run("GET https://www.ebi.ac.uk/eqtl/api/v2/datasets did not return usable dataset metadata.")
}

add_status("## API availability")
add_status("")
add_status("- `GET https://www.ebi.ac.uk/eqtl/api/v2/datasets`: available")
add_status("- GTEx-filtered dataset rows retrieved: ", nrow(datasets))
add_status("")

select_dataset <- function(target_row) {
  exact_hint <- datasets[
    datasets$dataset_id == target_row$dataset_id_hint &
      datasets$quant_method == "ge" &
      grepl("GTEx", datasets$study_label, ignore.case = TRUE) &
      grepl(target_row$sample_group_pattern, datasets$sample_group, ignore.case = TRUE),
    ,
    drop = FALSE
  ]
  if (nrow(exact_hint) > 0) {
    exact_hint$selection_rule <- "exact_dataset_id_hint_validated"
    return(exact_hint[order(exact_hint$dataset_id), , drop = FALSE][1, , drop = FALSE])
  }

  exact_label <- datasets[
    datasets$dataset_id == target_row$requested_dataset_id &
      datasets$quant_method == "ge" &
      grepl("GTEx", datasets$study_label, ignore.case = TRUE) &
      grepl(target_row$sample_group_pattern, datasets$sample_group, ignore.case = TRUE),
    ,
    drop = FALSE
  ]
  if (nrow(exact_label) > 0) {
    exact_label$selection_rule <- "exact_requested_dataset_id_validated"
    return(exact_label[order(exact_label$dataset_id), , drop = FALSE][1, , drop = FALSE])
  }

  ds <- datasets[
    datasets$quant_method == "ge" &
      grepl("GTEx", datasets$study_label, ignore.case = TRUE),
    ,
    drop = FALSE
  ]
  fuzzy <- ds[
    grepl(target_row$sample_group_pattern, ds$sample_group, ignore.case = TRUE),
    ,
    drop = FALSE
  ]
  if (nrow(fuzzy) == 0) {
    fuzzy <- ds[
      grepl(target_row$tissue_pattern, ds$tissue_label, ignore.case = TRUE) |
        grepl(target_row$tissue_pattern, ds$sample_group, ignore.case = TRUE),
      ,
      drop = FALSE
    ]
  }
  if (nrow(fuzzy) == 0) {
    return(data.frame())
  }
  fuzzy$selection_rule <- "fuzzy_keyword_match"
  fuzzy[order(fuzzy$dataset_id), , drop = FALSE][1, , drop = FALSE]
}

selected <- do.call(rbind, lapply(seq_len(nrow(targets)), function(i) {
  x <- select_dataset(targets[i, ])
  if (nrow(x) == 0) {
    return(data.frame(
      tissue = targets$tissue[[i]],
      requested_dataset_id = targets$requested_dataset_id[[i]],
      dataset_id_hint = targets$dataset_id_hint[[i]],
      dataset_id = NA_character_,
      study_label = NA_character_,
      quant_method = NA_character_,
      sample_group = NA_character_,
      tissue_label = NA_character_,
      sample_size = NA_real_,
      selection_rule = "not_found",
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    tissue = targets$tissue[[i]],
    requested_dataset_id = targets$requested_dataset_id[[i]],
    dataset_id_hint = targets$dataset_id_hint[[i]],
    dataset_id = x$dataset_id[[1]],
    study_label = x$study_label[[1]],
    quant_method = x$quant_method[[1]],
    sample_group = x$sample_group[[1]],
    tissue_label = x$tissue_label[[1]],
    sample_size = safe_numeric(x$sample_size[[1]]),
    selection_rule = x$selection_rule[[1]],
    stringsAsFactors = FALSE
  )
}))

write_csv(selected, file.path(validation_dir, "coloc_eqtlcatalogue_dataset_selection.csv"))

add_status("## Dataset selection")
add_status("")
for (i in seq_len(nrow(selected))) {
  add_status(
    "- ", selected$tissue[[i]], ": hint `", selected$dataset_id_hint[[i]],
    "`, requested label `", selected$requested_dataset_id[[i]],
    "`, using `", selected$dataset_id[[i]], "` (", selected$selection_rule[[i]], "; sample_group=",
    selected$sample_group[[i]], "; sample_size=", selected$sample_size[[i]], ")"
  )
}
add_status("")

if (any(is.na(selected$dataset_id))) {
  write_status()
}

fetch_page <- function(dataset_id, region, start, size) {
  url <- api_url(
    api_base_v2,
    paste0("/datasets/", dataset_id, "/associations"),
    list(pos = region, start = start, size = size)
  )
  as_records(api_get(url, paste0(dataset_id, " ", region, " start=", start), timeout_sec = 180))
}

fetch_eqtl_v2 <- function(dataset_id) {
  full_region <- paste0(chrom, ":", locus_start, "-", locus_end)
  test_url <- api_url(
    api_base_v2,
    paste0("/datasets/", dataset_id, "/associations"),
    list(pos = full_region, gene_id = tyk2_gene_id, start = 0, size = 20)
  )
  test <- api_get(test_url, paste0(dataset_id, " v2 full pos+gene test"), timeout_sec = 60)
  if (!is_api_error(test)) {
    tx <- as_records(test)
    if (nrow(tx) > 0) {
      append_log(run_log_path, dataset_id, " | v2 full pos+gene returned rows")
    }
  } else {
    append_log(run_log_path, dataset_id, " | v2 pos+gene unavailable; falling back to chunked pos-only")
  }

  all_rows <- list()
  idx <- 1L
  starts <- seq(locus_start, locus_end, by = chunk_size)
  for (chunk_start in starts) {
    chunk_end <- min(locus_end, chunk_start + chunk_size - 1L)
    region <- paste0(chrom, ":", chunk_start, "-", chunk_end)
    offset <- 0L
    repeat {
      page <- fetch_page(dataset_id, region, offset, page_size)
      if (nrow(page) == 0) {
        break
      }
      all_rows[[idx]] <- page
      idx <- idx + 1L
      if (nrow(page) < page_size) {
        break
      }
      offset <- offset + page_size
      Sys.sleep(0.1)
    }
  }
  if (length(all_rows) == 0) {
    return(data.frame())
  }
  out <- data.table::rbindlist(all_rows, fill = TRUE)
  out <- as.data.frame(out)
  out <- out[strip_ensembl_version(out$gene_id) == tyk2_gene_id, , drop = FALSE]
  out <- unique(out)
  out
}

fetch_eqtl_v3_gene_fallback <- function(dataset_id) {
  out <- list()
  idx <- 1L
  offset <- 0L
  complete <- TRUE
  reason <- NA_character_
  repeat {
    url <- api_url(
      api_base_v3,
      paste0("/datasets/", dataset_id, "/associations"),
      list(gene_id = tyk2_gene_id, chromosome = chrom, start = offset, size = page_size)
    )
    raw <- api_get(url, paste0(dataset_id, " v3 gene fallback start=", offset), timeout_sec = 240)
    if (is_api_error(raw)) {
      complete <- FALSE
      reason <- paste0("v3 dataset gene pagination failed at start=", offset, ": ", raw$error)
      break
    }
    page <- as_records(raw)
    if (nrow(page) == 0) {
      break
    }
    out[[idx]] <- page
    idx <- idx + 1L
    if (nrow(page) < page_size) {
      break
    }
    offset <- offset + page_size
    Sys.sleep(0.1)
  }
  if (length(out) == 0) {
    empty <- data.frame()
    attr(empty, "api_complete") <- complete
    attr(empty, "api_reason") <- reason
    return(empty)
  }
  x <- as.data.frame(data.table::rbindlist(out, fill = TRUE))
  x <- x[x$position >= locus_start & x$position <= locus_end, , drop = FALSE]
  x <- unique(x)
  attr(x, "api_complete") <- complete
  attr(x, "api_reason") <- reason
  x
}

standardize_eqtl <- function(x, sample_size) {
  if (nrow(x) == 0) {
    return(data.frame())
  }
  need <- c("rsid", "variant", "pvalue", "beta", "se", "chromosome", "position", "ref", "alt", "maf")
  for (nm in setdiff(need, names(x))) {
    x[[nm]] <- NA
  }
  out <- data.frame(
    rsid = as.character(x$rsid),
    variant = as.character(x$variant),
    chr = normalize_chr(x$chromosome),
    pos = as.integer(x$position),
    ref = toupper(as.character(x$ref)),
    alt = toupper(as.character(x$alt)),
    pvalue = standardize_p(x$pvalue),
    beta = safe_numeric(x$beta),
    se = safe_numeric(x$se),
    maf = safe_numeric(x$maf),
    N = safe_numeric(sample_size),
    stringsAsFactors = FALSE
  )
  out$key <- make_variant_key(out$chr, out$pos, out$ref, out$alt)
  out <- out[complete.cases(out[, c("chr", "pos", "ref", "alt", "beta", "se")]), , drop = FALSE]
  out <- out[is.finite(out$se) & out$se > 0, , drop = FALSE]
  out <- out[!duplicated(out$key), , drop = FALSE]
  out
}

read_finngen <- function() {
  locus_path <- file.path(validation_dir, "finngen_tyk2_locus_psoriasis.tsv.gz")
  check_file(locus_path, "FinnGen TYK2 psoriasis locus")
  g <- data.table::fread(locus_path, data.table = FALSE)
  cols <- names(g)
  chr_col <- find_col(cols, c("chr", "chromosome", "#chrom"), "FinnGen chromosome")
  pos_col <- find_col(cols, c("pos", "position", "bp"), "FinnGen position")
  ref_col <- find_col(cols, c("ref", "nea", "non_effect_allele"), "FinnGen ref")
  alt_col <- find_col(cols, c("alt", "ea", "effect_allele"), "FinnGen alt/effect allele")
  rsid_col <- find_col(cols, c("rsid", "rs_id", "snp"), "FinnGen rsid")
  beta_col <- find_col(cols, c("beta", "b", "effect", "log_or"), "FinnGen beta")
  se_col <- find_col(cols, c("se", "stderr", "standard_error"), "FinnGen se")
  p_col <- find_col(cols, c("pval", "p", "pvalue"), "FinnGen p")
  eaf_col <- find_col(cols, c("eaf", "af_alt", "effect_allele_frequency"), "FinnGen eaf")

  endpoints <- read_csv_required(cfg$data$finngen$outcome_index_path, "FinnGen outcome index")
  endpoint_row <- endpoints[endpoints$endpoint == "L12_PSORIASIS" | endpoints$disease_query == "psoriasis", , drop = FALSE]
  if (nrow(endpoint_row) == 0) {
    fail_run("Could not find psoriasis cases/controls in the FinnGen outcome index.")
  }
  ncases <- safe_numeric(endpoint_row$cases[[1]])
  ncontrols <- safe_numeric(endpoint_row$controls[[1]])
  if (!is.finite(ncases) || !is.finite(ncontrols) || ncases <= 0 || ncontrols <= 0) {
    fail_run("FinnGen psoriasis cases/controls are missing or invalid.")
  }

  out <- data.frame(
    rsid = as.character(g[[rsid_col]]),
    chr = normalize_chr(g[[chr_col]]),
    pos = as.integer(g[[pos_col]]),
    ref = toupper(as.character(g[[ref_col]])),
    alt = toupper(as.character(g[[alt_col]])),
    beta = safe_numeric(g[[beta_col]]),
    se = safe_numeric(g[[se_col]]),
    pvalue = standardize_p(g[[p_col]]),
    eaf = safe_numeric(g[[eaf_col]]),
    Ncases = ncases,
    Ncontrols = ncontrols,
    stringsAsFactors = FALSE
  )
  out$key <- make_variant_key(out$chr, out$pos, out$ref, out$alt)
  out
}

merge_loci <- function(eqtl, gwas) {
  exact <- merge(eqtl, gwas, by = "key", suffixes = c(".eqtl", ".gwas"))
  if (nrow(exact) > 0) {
    exact$allele_alignment <- "matched_chr_pos_ref_alt"
    exact$beta.gwas.aligned <- exact$beta.gwas
    exact$eaf.aligned <- exact$eaf
  }

  eqtl$swap_key <- make_variant_key(eqtl$chr, eqtl$pos, eqtl$alt, eqtl$ref)
  swapped <- merge(eqtl, gwas, by.x = "swap_key", by.y = "key", suffixes = c(".eqtl", ".gwas"))
  if (nrow(swapped) > 0) {
    if ("key.eqtl" %in% names(swapped)) {
      swapped$key <- swapped$key.eqtl
    } else if ("key" %in% names(swapped)) {
      swapped$key <- swapped$key
    } else {
      swapped$key <- make_variant_key(swapped$chr.eqtl, swapped$pos.eqtl, swapped$ref.eqtl, swapped$alt.eqtl)
    }
    swapped$allele_alignment <- "flipped_gwas_to_eqtl_alt"
    swapped$beta.gwas.aligned <- -swapped$beta.gwas
    swapped$eaf.aligned <- 1 - swapped$eaf
  }

  out <- data.table::rbindlist(list(exact, swapped), fill = TRUE)
  out <- as.data.frame(out)
  if (nrow(out) == 0) {
    return(out)
  }
  out <- out[!duplicated(out$key), , drop = FALSE]
  keep <- complete.cases(out[, c("beta.eqtl", "se.eqtl", "beta.gwas.aligned", "se.gwas")])
  out <- out[keep, , drop = FALSE]
  out <- out[is.finite(out$se.eqtl) & out$se.eqtl > 0 & is.finite(out$se.gwas) & out$se.gwas > 0, , drop = FALSE]
  out
}

run_coloc_one <- function(row, gwas) {
  tissue <- row$tissue
  out_csv <- file.path(validation_dir, paste0("coloc_tyk2_psoriasis_", tissue, ".csv"))
  fig_path <- file.path(fig_validation_dir, paste0("coloc_", tissue, "_sensitivity.png"))
  if (file.exists(fig_path)) {
    file.remove(fig_path)
  }

  if (is.na(row$dataset_id)) {
    res <- data.frame(
      tissue = tissue,
      dataset_id = NA_character_,
      status = "unavailable",
      reason = "Dataset not found in eQTL Catalogue REST API.",
      n_eqtl = 0L,
      n_overlap = 0L,
      PP.H0 = NA_real_,
      PP.H1 = NA_real_,
      PP.H2 = NA_real_,
      PP.H3 = NA_real_,
      PP.H4 = NA_real_,
      coloc_support = NA_character_,
      stringsAsFactors = FALSE
    )
    write_csv(res, out_csv)
    return(res)
  }

  eqtl_raw <- fetch_eqtl_v2(row$dataset_id)
  fetch_method <- "v2_chunked_pos_local_gene_filter"
  if (nrow(eqtl_raw) == 0) {
    append_log(run_log_path, row$dataset_id, " | v2 returned no TYK2 rows; trying v3 gene fallback")
    eqtl_raw <- fetch_eqtl_v3_gene_fallback(row$dataset_id)
    fetch_method <- "v3_dataset_gene_chromosome_fallback"
    if (identical(attr(eqtl_raw, "api_complete"), FALSE)) {
      res <- data.frame(
        tissue = tissue,
        dataset_id = row$dataset_id,
        status = "unavailable",
        reason = attr(eqtl_raw, "api_reason"),
        n_eqtl = nrow(eqtl_raw),
        n_overlap = 0L,
        PP.H0 = NA_real_,
        PP.H1 = NA_real_,
        PP.H2 = NA_real_,
        PP.H3 = NA_real_,
        PP.H4 = NA_real_,
        coloc_support = NA_character_,
        stringsAsFactors = FALSE
      )
      write_csv(res, out_csv)
      return(res)
    }
  }

  eqtl <- standardize_eqtl(eqtl_raw, row$sample_size)
  eqtl_path <- file.path(validation_dir, paste0("eqtlcatalogue_tyk2_", tissue, ".csv"))
  write_csv(eqtl, eqtl_path)

  if (nrow(eqtl) == 0) {
    res <- data.frame(
      tissue = tissue,
      dataset_id = row$dataset_id,
      status = "unavailable",
      reason = paste0("No TYK2 eQTL rows returned by ", fetch_method, "."),
      n_eqtl = 0L,
      n_overlap = 0L,
      PP.H0 = NA_real_,
      PP.H1 = NA_real_,
      PP.H2 = NA_real_,
      PP.H3 = NA_real_,
      PP.H4 = NA_real_,
      coloc_support = NA_character_,
      stringsAsFactors = FALSE
    )
    write_csv(res, out_csv)
    return(res)
  }

  merged <- merge_loci(eqtl, gwas)
  merged_path <- file.path(validation_dir, paste0("coloc_tyk2_psoriasis_", tissue, "_merged_variants.csv"))
  write_csv(merged, merged_path)

  if (nrow(merged) < min_overlap) {
    res <- data.frame(
      tissue = tissue,
      dataset_id = row$dataset_id,
      status = "unavailable",
      reason = paste0("Fewer than ", min_overlap, " overlapping complete variants: ", nrow(merged), "."),
      n_eqtl = nrow(eqtl),
      n_overlap = nrow(merged),
      PP.H0 = NA_real_,
      PP.H1 = NA_real_,
      PP.H2 = NA_real_,
      PP.H3 = NA_real_,
      PP.H4 = NA_real_,
      coloc_support = NA_character_,
      stringsAsFactors = FALSE
    )
    write_csv(res, out_csv)
    return(res)
  }

  ncases <- median(merged$Ncases, na.rm = TRUE)
  ncontrols <- median(merged$Ncontrols, na.rm = TRUE)
  d1 <- list(
    beta = merged$beta.eqtl,
    varbeta = merged$se.eqtl^2,
    snp = ifelse(
      !is.na(merged$rsid.eqtl) & nzchar(merged$rsid.eqtl) & !duplicated(merged$rsid.eqtl),
      merged$rsid.eqtl,
      merged$key
    ),
    N = median(merged$N.eqtl, na.rm = TRUE),
    type = "quant",
    sdY = 1
  )
  maf_eqtl <- pmin(merged$maf, 1 - merged$maf)
  if (any(is.finite(maf_eqtl))) {
    d1$MAF <- maf_eqtl
  }
  d2 <- list(
    beta = merged$beta.gwas.aligned,
    varbeta = merged$se.gwas^2,
    snp = d1$snp,
    N = ncases + ncontrols,
    s = ncases / (ncases + ncontrols),
    type = "cc"
  )
  maf_gwas <- pmin(merged$eaf.aligned, 1 - merged$eaf.aligned)
  if (any(is.finite(maf_gwas))) {
    d2$MAF <- maf_gwas
  }

  coloc_res <- coloc::coloc.abf(d1, d2, p1 = cfg$coloc$p1, p2 = cfg$coloc$p2, p12 = cfg$coloc$p12)
  pp <- as.numeric(coloc_res$summary[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")])
  names(pp) <- c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")
  support <- if (pp[["PP.H4"]] >= cfg$coloc$pp_h4_strong) {
    "strong_support"
  } else if (pp[["PP.H4"]] >= cfg$coloc$pp_h4_support) {
    "support"
  } else {
    "not_supported"
  }

  sensitivity_status <- "skipped_pph4_not_supported"
  if (pp[["PP.H4"]] >= cfg$coloc$pp_h4_support) {
    tryCatch({
      grDevices::png(fig_path, width = 1800, height = 1400, res = 200)
      coloc::sensitivity(coloc_res, rule = paste0("H4 >= ", cfg$coloc$pp_h4_support))
      grDevices::dev.off()
      sensitivity_status <- "complete"
    }, error = function(e) {
      if (names(grDevices::dev.cur()) != "null device") {
        grDevices::dev.off()
      }
      sensitivity_status <<- paste0("unavailable: ", conditionMessage(e))
    })
  }

  res <- data.frame(
    tissue = tissue,
    dataset_id = row$dataset_id,
    status = "complete",
    reason = NA_character_,
    n_eqtl = nrow(eqtl),
    n_overlap = nrow(merged),
    PP.H0 = pp[["PP.H0"]],
    PP.H1 = pp[["PP.H1"]],
    PP.H2 = pp[["PP.H2"]],
    PP.H3 = pp[["PP.H3"]],
    PP.H4 = pp[["PP.H4"]],
    coloc_support = support,
    sensitivity_plot = sensitivity_status,
    fetch_method = fetch_method,
    stringsAsFactors = FALSE
  )
  write_csv(res, out_csv)
  res
}

gwas <- read_finngen()
results <- do.call(rbind, lapply(seq_len(nrow(selected)), function(i) run_coloc_one(selected[i, ], gwas)))
write_csv(results, summary_path)

add_status("## Coloc run")
add_status("")
for (i in seq_len(nrow(results))) {
  pp <- ifelse(is.na(results$PP.H4[[i]]), "NA", signif(results$PP.H4[[i]], 4))
  add_status(
    "- ", results$tissue[[i]], ": ", results$status[[i]],
    "; overlap=", results$n_overlap[[i]],
    "; PP.H4=", pp,
    "; support=", results$coloc_support[[i]]
  )
  if (!is.na(results$reason[[i]])) {
    add_status("  - Reason: ", results$reason[[i]])
  }
}
add_status("")

complete <- results[results$status == "complete", , drop = FALSE]
max_pph4 <- if (nrow(complete) > 0) max(complete$PP.H4, na.rm = TRUE) else NA_real_
new_recommendation <- if (is.finite(max_pph4) && max_pph4 >= cfg$coloc$pp_h4_strong) {
  "coloc_strong_support_pending_external_replication"
} else if (is.finite(max_pph4) && max_pph4 >= cfg$coloc$pp_h4_support) {
  "coloc_support_pending_external_replication"
} else {
  "needs_external_replication"
}

add_status("## Recommendation rule")
add_status("")
add_status("- PP.H4 >= ", cfg$coloc$pp_h4_support, ": support colocalization")
add_status("- PP.H4 >= ", cfg$coloc$pp_h4_strong, ": strong support")
add_status("- Recommendation after eQTL Catalogue coloc: **", new_recommendation, "**")
add_status("")
write_status()

append_validation_summary <- function() {
  top_summary <- file.path(validation_dir, "top_hit_validation_summary.md")
  if (!file.exists(top_summary)) {
    append_log(run_log_path, "top_hit_validation_summary.md not found; skipping append")
    return(invisible(FALSE))
  }
  backup_dir <- file.path(
    "backups",
    paste0("coloc-eqtlcatalogue-", format(Sys.time(), "%Y-%m-%d-%H%M")),
    "results",
    "validation"
  )
  ensure_dir(backup_dir)
  file.copy(top_summary, file.path(backup_dir, basename(top_summary)), overwrite = TRUE)

  existing <- readLines(top_summary, warn = FALSE)
  prior_section <- grep("^## 9\\. eQTL Catalogue REST TYK2 coloc", existing)
  if (length(prior_section) > 0) {
    existing <- existing[seq_len(prior_section[[1]] - 1L)]
    writeLines(existing, top_summary, useBytes = TRUE)
  }

  block <- c(
    "",
    paste0("## 9. eQTL Catalogue REST TYK2 coloc (", format(Sys.time(), "%Y-%m-%d %H:%M"), ")"),
    "",
    paste0(
      "Public eQTL Catalogue REST API was queried for GTEx gene-expression datasets at ",
      "chr19:", locus_start, "-", locus_end, " and TYK2 (", tyk2_gene_id, ")."
    ),
    "",
    "| tissue | dataset_id | status | overlap variants | PP.H4 | coloc interpretation |",
    "| --- | --- | --- | ---: | ---: | --- |"
  )
  for (i in seq_len(nrow(results))) {
    pp <- ifelse(is.na(results$PP.H4[[i]]), "NA", format(signif(results$PP.H4[[i]], 4), scientific = TRUE))
    block <- c(block, paste0(
      "| ", results$tissue[[i]], " | ", results$dataset_id[[i]], " | ", results$status[[i]],
      " | ", results$n_overlap[[i]], " | ", pp, " | ", results$coloc_support[[i]], " |"
    ))
  }
  block <- c(
    block,
    "",
    paste0("Recommendation after this coloc attempt: **", new_recommendation, "**."),
    "If any tissue is unavailable, the original external-replication caution remains for that tissue.",
    ""
  )
  cat(paste(block, collapse = "\n"), file = top_summary, append = TRUE)
  invisible(TRUE)
}

append_validation_summary()

cat("\n=== results/validation/coloc_tyk2_psoriasis_summary.csv ===\n")
print(data.table::fread(summary_path, data.table = FALSE))
cat("\n=== results/validation/coloc_eqtlcatalogue_status.md ===\n")
cat(paste(readLines(status_path, warn = FALSE), collapse = "\n"), "\n")
cat("\n=== top_hit_validation_summary.md appended section ===\n")
tail_lines <- tail(readLines(file.path(validation_dir, "top_hit_validation_summary.md"), warn = FALSE), 30)
cat(paste(tail_lines, collapse = "\n"), "\n")
quit(save = "no", status = 0)
