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

status_path <- file.path(validation_dir, "coloc_pde4a_status.md")
summary_path <- file.path(validation_dir, "coloc_pde4a_summary.csv")
run_log_path <- file.path(log_dir, "coloc_pde4a_api.log")

api_base_v2 <- "https://www.ebi.ac.uk/eqtl/api/v2"
api_base_v3 <- "https://www.ebi.ac.uk/eqtl/api/v3"
page_size <- 1000L
chunk_size <- 100000L
min_overlap <- 50L

result_columns <- c(
  "tissue", "dataset_id", "status", "reason", "n_eqtl", "n_overlap",
  "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "coloc_support",
  "sensitivity_plot", "fetch_method"
)

targets <- data.frame(
  gene = "PDE4A",
  gene_id = "ENSG00000065989",
  chrom = "19",
  locus_start = 9933000L,
  locus_end = 10954000L,
  instrument_chr = "19",
  instrument_pos = 10430708L,
  instrument_ref = "C",
  instrument_alt = "G",
  tissue = "skin_sun_exposed_lower_leg_skin",
  requested_dataset_id = "GTEx_V8.skin_sun_exposed_lower_leg_skin",
  dataset_id_hint = "QTD000316",
  sample_group_pattern = "skin_sun_exposed",
  tissue_pattern = "skin",
  stringsAsFactors = FALSE
)

status_lines <- c(
  "# eQTL Catalogue PDE4A coloc status",
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

empty_result <- function(row, reason, n_eqtl = 0L, n_overlap = 0L, fetch_method = NA_character_) {
  out <- data.frame(
    tissue = row$tissue,
    dataset_id = if ("dataset_id" %in% names(row)) row$dataset_id else NA_character_,
    status = "unavailable",
    reason = reason,
    n_eqtl = as.integer(n_eqtl),
    n_overlap = as.integer(n_overlap),
    PP.H0 = NA_real_,
    PP.H1 = NA_real_,
    PP.H2 = NA_real_,
    PP.H3 = NA_real_,
    PP.H4 = NA_real_,
    coloc_support = NA_character_,
    sensitivity_plot = NA_character_,
    fetch_method = fetch_method,
    stringsAsFactors = FALSE
  )
  out[, result_columns]
}

if (length(missing) > 0) {
  fail_run(paste0("Missing required R packages: ", paste(missing, collapse = ", ")))
}

append_log(run_log_path, "Starting PDE4A eQTL Catalogue REST coloc run")

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
    if (FALSE && is.null(res) && nzchar(Sys.which("curl.exe"))) {
      tmp <- tempfile(pattern = "eqtlcatalogue_", tmpdir = log_dir, fileext = ".json")
      args <- c(
        "--location", "--fail", "--silent", "--show-error", "--ssl-no-revoke",
        "--max-time", as.character(timeout_sec), "--output", tmp, url
      )
      code <- tryCatch(system2(Sys.which("curl.exe"), args = args), error = function(e) {
        last_error <<- conditionMessage(e)
        1L
      })
      if (identical(code, 0L) && file.exists(tmp) && file.info(tmp)$size > 0) {
        res <- tryCatch(
          jsonlite::fromJSON(tmp, simplifyDataFrame = TRUE),
          error = function(e) {
            last_error <<- conditionMessage(e)
            NULL
          }
        )
      } else if (file.exists(tmp)) {
        last_error <- paste0(last_error, " curl_exit=", code)
      }
      if (file.exists(tmp)) {
        unlink(tmp)
      }
    }
    if (is.null(res) && file.exists(file.path(script_dir, "eqtl_api_get.py"))) {
      tmp <- tempfile(pattern = "eqtlcatalogue_py_", tmpdir = log_dir, fileext = ".json")
      py <- Sys.getenv("PYTHON", unset = "")
      if (!nzchar(py)) {
        candidates <- Sys.which(c("python3", "python"))
        candidates <- unname(candidates[nzchar(candidates)])
        py <- if (length(candidates)) candidates[[1]] else ""
      }
      if (nzchar(py)) {
        code <- tryCatch(
          system2(py, args = c(file.path(script_dir, "eqtl_api_get.py"), url, tmp, as.character(timeout_sec))),
          error = function(e) {
            last_error <<- conditionMessage(e)
            1L
          }
        )
        if (identical(code, 0L) && file.exists(tmp) && file.info(tmp)$size > 0) {
          res <- tryCatch(
            jsonlite::fromJSON(tmp, simplifyDataFrame = TRUE),
            error = function(e) {
              last_error <<- conditionMessage(e)
              NULL
            }
          )
        } else {
          last_error <- paste0(last_error, " python_requests_exit=", code)
        }
      }
      if (file.exists(tmp)) {
        unlink(tmp)
      }
    }
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
    datasets$quant_method == "ge" & grepl("GTEx", datasets$study_label, ignore.case = TRUE),
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
      gene = targets$gene[[i]],
      gene_id = targets$gene_id[[i]],
      chrom = targets$chrom[[i]],
      locus_start = targets$locus_start[[i]],
      locus_end = targets$locus_end[[i]],
      instrument_chr = targets$instrument_chr[[i]],
      instrument_pos = targets$instrument_pos[[i]],
      instrument_ref = targets$instrument_ref[[i]],
      instrument_alt = targets$instrument_alt[[i]],
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
    gene = targets$gene[[i]],
    gene_id = targets$gene_id[[i]],
    chrom = targets$chrom[[i]],
    locus_start = targets$locus_start[[i]],
    locus_end = targets$locus_end[[i]],
    instrument_chr = targets$instrument_chr[[i]],
    instrument_pos = targets$instrument_pos[[i]],
    instrument_ref = targets$instrument_ref[[i]],
    instrument_alt = targets$instrument_alt[[i]],
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

write_csv(selected, file.path(validation_dir, "coloc_pde4a_dataset_selection.csv"))

add_status("## Dataset selection")
add_status("")
for (i in seq_len(nrow(selected))) {
  add_status(
    "- ", selected$gene[[i]], " / ", selected$tissue[[i]], ": requested `",
    selected$requested_dataset_id[[i]], "`, hint `", selected$dataset_id_hint[[i]],
    "`, using `", selected$dataset_id[[i]],
    "` (", selected$selection_rule[[i]], "; sample_group=", selected$sample_group[[i]],
    "; sample_size=", selected$sample_size[[i]], ")"
  )
}
add_status("")

fetch_page <- function(dataset_id, region, start, size) {
  url <- api_url(
    api_base_v2,
    paste0("/datasets/", dataset_id, "/associations"),
    list(pos = region, start = start, size = size)
  )
  as_records(api_get(url, paste0(dataset_id, " ", region, " start=", start), timeout_sec = 180))
}

fetch_eqtl_v2 <- function(dataset_id, target_row) {
  full_region <- paste0(target_row$chrom, ":", target_row$locus_start, "-", target_row$locus_end)
  test_url <- api_url(
    api_base_v2,
    paste0("/datasets/", dataset_id, "/associations"),
    list(pos = full_region, gene_id = target_row$gene_id, start = 0, size = 20)
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
  starts <- seq(target_row$locus_start, target_row$locus_end, by = chunk_size)
  for (chunk_start in starts) {
    chunk_end <- min(target_row$locus_end, chunk_start + chunk_size - 1L)
    region <- paste0(target_row$chrom, ":", chunk_start, "-", chunk_end)
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
  out <- out[strip_ensembl_version(out$gene_id) == target_row$gene_id, , drop = FALSE]
  out <- unique(out)
  out
}

fetch_eqtl_v3_gene_fallback <- function(dataset_id, target_row) {
  out <- list()
  idx <- 1L
  offset <- 0L
  complete <- TRUE
  reason <- NA_character_
  repeat {
    url <- api_url(
      api_base_v3,
      paste0("/datasets/", dataset_id, "/associations"),
      list(gene_id = target_row$gene_id, chromosome = target_row$chrom, start = offset, size = page_size)
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
  x <- x[x$position >= target_row$locus_start & x$position <= target_row$locus_end, , drop = FALSE]
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

get_psoriasis_counts <- function() {
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
  list(ncases = ncases, ncontrols = ncontrols)
}

ensure_finngen_locus <- function(target_row) {
  outcome_path <- file.path(cfg$data$finngen$outcome_summary_dir, "finngen_R13_L12_PSORIASIS.gz")
  outcome_url <- "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_L12_PSORIASIS.gz"
  download_if_missing(outcome_url, outcome_path, "FinnGen R13 psoriasis summary statistics")

  out_path <- file.path(validation_dir, paste0("finngen_", tolower(target_row$gene), "_locus_psoriasis.tsv.gz"))
  if (file.exists(out_path) && file.info(out_path)$size > 0) {
    return(out_path)
  }

  append_log(run_log_path, "Subsetting FinnGen locus for ", target_row$gene)
  con <- gzfile(outcome_path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1L, warn = FALSE)
  if (length(header) == 0) {
    fail_run("FinnGen psoriasis summary statistics file is empty.")
  }

  chunks <- list()
  idx <- 1L
  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (length(lines) == 0) {
      break
    }
    dt <- data.table::fread(text = paste(c(header, lines), collapse = "\n"), data.table = FALSE, showProgress = FALSE)
    chr_col <- find_col(names(dt), c("chr", "chromosome", "#chrom"), "FinnGen chromosome")
    pos_col <- find_col(names(dt), c("pos", "position", "bp"), "FinnGen position")
    keep <- normalize_chr(dt[[chr_col]]) == target_row$chrom &
      as.integer(dt[[pos_col]]) >= target_row$locus_start &
      as.integer(dt[[pos_col]]) <= target_row$locus_end
    if (any(keep, na.rm = TRUE)) {
      chunks[[idx]] <- dt[keep, , drop = FALSE]
      idx <- idx + 1L
    }
  }

  if (length(chunks) == 0) {
    fail_run(paste0("No FinnGen psoriasis variants found for ", target_row$gene, " locus."))
  }

  g <- as.data.frame(data.table::rbindlist(chunks, fill = TRUE))
  cols <- names(g)
  chr_col <- find_col(cols, c("chr", "chromosome", "#chrom"), "FinnGen chromosome")
  pos_col <- find_col(cols, c("pos", "position", "bp"), "FinnGen position")
  ref_col <- find_col(cols, c("ref", "nea", "non_effect_allele"), "FinnGen ref")
  alt_col <- find_col(cols, c("alt", "ea", "effect_allele"), "FinnGen alt/effect allele")
  rsid_col <- find_col(cols, c("rsids", "rsid", "rs_id", "snp"), "FinnGen rsid")
  beta_col <- find_col(cols, c("beta", "b", "effect", "log_or"), "FinnGen beta")
  se_col <- find_col(cols, c("sebeta", "se", "stderr", "standard_error"), "FinnGen se")
  p_col <- find_col(cols, c("pval", "p", "pvalue"), "FinnGen p")
  eaf_col <- find_col(cols, c("af_alt", "eaf", "effect_allele_frequency"), "FinnGen eaf")

  out <- data.frame(
    chr = normalize_chr(g[[chr_col]]),
    pos = as.integer(g[[pos_col]]),
    ref = toupper(as.character(g[[ref_col]])),
    alt = toupper(as.character(g[[alt_col]])),
    rsid = as.character(g[[rsid_col]]),
    beta = safe_numeric(g[[beta_col]]),
    se = safe_numeric(g[[se_col]]),
    p = standardize_p(g[[p_col]]),
    eaf = safe_numeric(g[[eaf_col]]),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$chr, out$pos), , drop = FALSE]
  gz <- gzfile(out_path, open = "wt")
  on.exit(close(gz), add = TRUE)
  utils::write.table(out, gz, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  out_path
}

read_finngen <- function(target_row) {
  locus_path <- ensure_finngen_locus(target_row)
  g <- data.table::fread(locus_path, data.table = FALSE)
  counts <- get_psoriasis_counts()
  cols <- names(g)
  chr_col <- find_col(cols, c("chr", "chromosome", "#chrom"), "FinnGen chromosome")
  pos_col <- find_col(cols, c("pos", "position", "bp"), "FinnGen position")
  ref_col <- find_col(cols, c("ref", "nea", "non_effect_allele"), "FinnGen ref")
  alt_col <- find_col(cols, c("alt", "ea", "effect_allele"), "FinnGen alt/effect allele")
  rsid_col <- find_col(cols, c("rsid", "rsids", "rs_id", "snp"), "FinnGen rsid")
  beta_col <- find_col(cols, c("beta", "b", "effect", "log_or"), "FinnGen beta")
  se_col <- find_col(cols, c("se", "sebeta", "stderr", "standard_error"), "FinnGen se")
  p_col <- find_col(cols, c("p", "pval", "pvalue"), "FinnGen p")
  eaf_col <- find_col(cols, c("eaf", "af_alt", "effect_allele_frequency"), "FinnGen eaf")

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
    Ncases = counts$ncases,
    Ncontrols = counts$ncontrols,
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

safe_file_gene <- function(gene) tolower(gene)

run_coloc_one <- function(row) {
  gene_slug <- safe_file_gene(row$gene)
  tissue <- row$tissue
  out_csv <- file.path(validation_dir, paste0("coloc_", gene_slug, "_psoriasis_", tissue, ".csv"))
  merged_path <- file.path(validation_dir, paste0("coloc_", gene_slug, "_psoriasis_", tissue, "_merged_variants.csv"))
  fig_path <- file.path(fig_validation_dir, paste0("coloc_", gene_slug, "_", tissue, "_sensitivity.png"))
  if (file.exists(out_csv) && file.exists(merged_path)) {
    prior <- data.table::fread(out_csv, data.table = FALSE)
    if (nrow(prior) == 1L && identical(prior$status[[1]], "complete")) {
      return(prior[, result_columns])
    }
  }
  if (file.exists(fig_path)) {
    file.remove(fig_path)
  }

  if (is.na(row$dataset_id)) {
    res <- empty_result(row, "Dataset not found in eQTL Catalogue REST API.")
    write_csv(res, out_csv)
    return(res)
  }

  gwas <- read_finngen(row)
  instrument_key <- make_variant_key(row$instrument_chr, row$instrument_pos, row$instrument_ref, row$instrument_alt)
  if (!instrument_key %in% gwas$key) {
    res <- empty_result(row, paste0("MR instrument not found in FinnGen locus: ", instrument_key, "."))
    write_csv(res, out_csv)
    return(res)
  }

  eqtl_raw <- fetch_eqtl_v2(row$dataset_id, row)
  fetch_method <- "v2_chunked_pos_local_gene_filter"
  if (nrow(eqtl_raw) == 0) {
    append_log(run_log_path, row$dataset_id, " | v2 returned no ", row$gene, " rows; trying v3 gene fallback")
    eqtl_raw <- fetch_eqtl_v3_gene_fallback(row$dataset_id, row)
    fetch_method <- "v3_dataset_gene_chromosome_fallback"
    if (identical(attr(eqtl_raw, "api_complete"), FALSE)) {
      res <- empty_result(row, attr(eqtl_raw, "api_reason"), n_eqtl = nrow(eqtl_raw), fetch_method = fetch_method)
      write_csv(res, out_csv)
      return(res)
    }
  }

  eqtl <- standardize_eqtl(eqtl_raw, row$sample_size)
  eqtl_path <- file.path(validation_dir, paste0("eqtlcatalogue_", gene_slug, "_", tissue, ".csv"))
  write_csv(eqtl, eqtl_path)

  if (nrow(eqtl) == 0) {
    res <- empty_result(row, paste0("No ", row$gene, " eQTL rows returned by ", fetch_method, "."), fetch_method = fetch_method)
    write_csv(res, out_csv)
    return(res)
  }

  merged <- merge_loci(eqtl, gwas)
  write_csv(merged, merged_path)

  if (nrow(merged) < min_overlap) {
    res <- empty_result(
      row,
      paste0("Fewer than ", min_overlap, " overlapping complete variants: ", nrow(merged), "."),
      n_eqtl = nrow(eqtl),
      n_overlap = nrow(merged),
      fetch_method = fetch_method
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
    N = median(merged$N, na.rm = TRUE),
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
  write_csv(res[, result_columns], out_csv)
  res[, result_columns]
}

results <- do.call(rbind, lapply(seq_len(nrow(selected)), function(i) run_coloc_one(selected[i, ])))
write_csv(results[, result_columns], summary_path)

add_status("## Coloc run")
add_status("")
for (i in seq_len(nrow(results))) {
  pp <- ifelse(is.na(results$PP.H4[[i]]), "NA", signif(results$PP.H4[[i]], 4))
  add_status(
    "- ", targets$gene[[i]], " / ", results$tissue[[i]], ": ", results$status[[i]],
    "; overlap=", results$n_overlap[[i]],
    "; PP.H4=", pp,
    "; support=", results$coloc_support[[i]]
  )
  if (!is.na(results$reason[[i]])) {
    add_status("  - Reason: ", results$reason[[i]])
  }
}
add_status("")
write_status()

cat("\n=== results/validation/coloc_pde4a_summary.csv ===\n")
print(data.table::fread(summary_path, data.table = FALSE))
cat("\n=== results/validation/coloc_pde4a_status.md ===\n")
cat(paste(readLines(status_path, warn = FALSE), collapse = "\n"), "\n")
quit(save = "no", status = 0)
