suppressPackageStartupMessages({
  options(stringsAsFactors = FALSE)
})

args_file <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])
if (is.na(args_file) || !nzchar(args_file)) {
  args_file <- "scripts/10_opengwas_psoriasis_replication.R"
}
root <- normalizePath(file.path(dirname(args_file), ".."), winslash = "/", mustWork = TRUE)
validation_dir <- file.path(root, "results", "validation")
logs_dir <- file.path(root, "logs", "validation")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

status_path <- file.path(validation_dir, "opengwas_replication_status.md")
candidates_path <- file.path(validation_dir, "opengwas_psoriasis_candidates.csv")
replication_path <- file.path(validation_dir, "tyk2_opengwas_replication.csv")

write_status <- function(status, detail, extra = character()) {
  writeLines(c("# OpenGWAS psoriasis association-lookup status", "", paste0("Status: **", status, "**"), "", detail, extra, ""), status_path, useBytes = TRUE)
}

get_col <- function(df, choices, default = NA_character_) {
  hit <- choices[choices %in% names(df)]
  if (length(hit) == 0) return(rep(default, nrow(df)))
  as.character(df[[hit[[1]]]])
}

get_num <- function(x) suppressWarnings(as.numeric(x))

empty_candidates <- function() {
  data.frame(
    id = character(), trait = character(), population = character(), sample_size = numeric(),
    ncase = numeric(), ncontrol = numeric(), source = character(), has_summary_stats = character(),
    excluded_as_finngen = logical(), priority_score = numeric(), selected_for_replication = logical(),
    stringsAsFactors = FALSE
  )
}

empty_replication <- function() {
  data.frame(
    outcome_id = character(), outcome_trait = character(), tissue = character(),
    requested_rsid = character(), matched_variant = character(), association_type = character(),
    proxy_r2 = character(), effect_allele.outcome = character(), other_allele.outcome = character(),
    beta.outcome = character(), se.outcome = character(), p.outcome = character(), eaf.outcome = character(),
    beta.exposure = character(), se.exposure = character(), Wald_beta = character(), Wald_se = character(),
    Wald_p = character(), OR = character(), OR_lci = character(), OR_uci = character(),
    status = character(), note = character(), stringsAsFactors = FALSE
  )
}

audit_path <- file.path(validation_dir, "top_hit_instrument_audit.csv")
if (!file.exists(audit_path)) {
  write.csv(empty_candidates(), candidates_path, row.names = FALSE)
  write.csv(empty_replication(), replication_path, row.names = FALSE)
  write_status("unavailable", "top_hit_instrument_audit.csv is unavailable. Run scripts/10_top_hit_validation.py first.")
  quit(status = 0)
}

if (!requireNamespace("ieugwasr", quietly = TRUE)) {
  write.csv(empty_candidates(), candidates_path, row.names = FALSE)
  write.csv(empty_replication(), replication_path, row.names = FALSE)
  write_status("unavailable", "ieugwasr is not installed in the local R library. No OpenGWAS query was attempted.")
  quit(status = 0)
}

if (!nzchar(Sys.getenv("OPENGWAS_JWT"))) {
  write.csv(empty_candidates(), candidates_path, row.names = FALSE)
  write.csv(empty_replication(), replication_path, row.names = FALSE)
  write_status("unavailable", "OPENGWAS_JWT is not set in the R environment. No OpenGWAS query was attempted.")
  quit(status = 0)
}

audit <- read.csv(audit_path, check.names = FALSE)
tyk2 <- audit[audit$target_gene == "TYK2" & audit$disease == "psoriasis", , drop = FALSE]
rsid_col <- "rsid/rsids"
if (nrow(tyk2) == 0 || !(rsid_col %in% names(tyk2))) {
  write.csv(empty_candidates(), candidates_path, row.names = FALSE)
  write.csv(empty_replication(), replication_path, row.names = FALSE)
  write_status("unavailable", "TYK2 psoriasis rows or rsid column are unavailable in the audit table.")
  quit(status = 0)
}

tyk2$requested_rsid <- vapply(strsplit(as.character(tyk2[[rsid_col]]), ","), `[`, character(1), 1)
tyk2 <- tyk2[!is.na(tyk2$requested_rsid) & nzchar(tyk2$requested_rsid) & tyk2$requested_rsid != "unavailable", , drop = FALSE]
requested_rsids <- unique(tyk2$requested_rsid)
if (length(requested_rsids) == 0) {
  write.csv(empty_candidates(), candidates_path, row.names = FALSE)
  write.csv(empty_replication(), replication_path, row.names = FALSE)
  write_status("unavailable", "No usable TYK2 rsIDs were available for OpenGWAS lookup.")
  quit(status = 0)
}

phewas <- tryCatch(
  ieugwasr::phewas(requested_rsids, pval = 0.01, timeout = 60),
  error = function(e) e
)
if (inherits(phewas, "error")) {
  write.csv(empty_candidates(), candidates_path, row.names = FALSE)
  write.csv(empty_replication(), replication_path, row.names = FALSE)
  write_status("unavailable", paste("OpenGWAS PheWAS candidate search failed.", paste0("Error: ", conditionMessage(phewas))))
  quit(status = 0)
}

if (is.null(phewas) || nrow(phewas) == 0 || !("trait" %in% names(phewas)) || !("id" %in% names(phewas))) {
  write.csv(empty_candidates(), candidates_path, row.names = FALSE)
  write.csv(empty_replication(), replication_path, row.names = FALSE)
  write_status("unavailable", "OpenGWAS PheWAS returned no usable candidate metadata for the TYK2 variants.")
  quit(status = 0)
}

psoriasis_phewas <- phewas[grepl("psoriasis|psoriatic", phewas[["trait"]], ignore.case = TRUE), , drop = FALSE]
if (nrow(psoriasis_phewas) == 0) {
  write.csv(empty_candidates(), candidates_path, row.names = FALSE)
  write.csv(empty_replication(), replication_path, row.names = FALSE)
  write_status("unavailable", "OpenGWAS PheWAS returned no psoriasis/psoriatic associations for the TYK2 variants at p <= 0.01.")
  quit(status = 0)
}

candidate_ids_all <- unique(psoriasis_phewas$id)
candidate_info <- tryCatch(
  ieugwasr::gwasinfo(candidate_ids_all, timeout = 60),
  error = function(e) e
)
if (inherits(candidate_info, "error") || is.null(candidate_info) || nrow(candidate_info) == 0) {
  candidate_info <- unique(psoriasis_phewas[c("id", "trait")])
}

ids <- get_col(candidate_info, c("id"))
traits <- get_col(candidate_info, c("trait"))
population <- get_col(candidate_info, c("population", "pop", "ancestry"))
sample_size <- get_num(get_col(candidate_info, c("sample_size", "n", "n_total", "total_n", "ncase")))
ncase <- get_num(get_col(candidate_info, c("ncase", "cases", "n_cases")))
ncontrol <- get_num(get_col(candidate_info, c("ncontrol", "controls", "n_controls")))
source <- get_col(candidate_info, c("source", "author", "pmid", "consortium", "mr"))
has_summary_stats <- get_col(candidate_info, c("mr", "has_summary_stats", "public", "access"))

text_for_exclusion <- apply(candidate_info, 1, paste, collapse = " ")
excluded_as_finngen <- grepl("FinnGen|FINNGEN|finn-b|R[0-9]+_", text_for_exclusion, ignore.case = TRUE) | grepl("^finn", ids, ignore.case = TRUE)
is_eur <- grepl("EUR|European|Europe", population, ignore.case = TRUE)
cases_ok <- !is.na(ncase) & ncase > 1000
summary_ok <- is.na(has_summary_stats) | has_summary_stats == "" | grepl("1|TRUE|yes|public|available", has_summary_stats, ignore.case = TRUE)
priority_score <- (1000 * is_eur) + (500 * cases_ok) + (100 * summary_ok) + ifelse(is.na(ncase), 0, pmin(ncase, 100000) / 1000) + ifelse(is.na(sample_size), 0, pmin(sample_size, 1000000) / 100000)

candidates <- data.frame(
  id = ids, trait = traits, population = population, sample_size = sample_size,
  ncase = ncase, ncontrol = ncontrol, source = source, has_summary_stats = has_summary_stats,
  excluded_as_finngen = excluded_as_finngen, priority_score = priority_score,
  selected_for_replication = FALSE, stringsAsFactors = FALSE
)
candidates <- candidates[!is.na(candidates$id) & nzchar(candidates$id), , drop = FALSE]
candidates <- candidates[order(candidates$excluded_as_finngen, -candidates$priority_score, candidates$id), , drop = FALSE]
selected_ids <- head(unique(candidates$id[!candidates$excluded_as_finngen]), 10)
candidates$selected_for_replication <- candidates$id %in% selected_ids
write.csv(candidates, candidates_path, row.names = FALSE)

if (length(selected_ids) == 0) {
  write.csv(empty_replication(), replication_path, row.names = FALSE)
  write_status("unavailable", "Only FinnGen-like psoriasis outcomes were identifiable; no eligible non-FinnGen OpenGWAS record was selected.")
  quit(status = 0)
}

extract_assoc <- function(snps, outcomes, proxies) {
  tryCatch(
    ieugwasr::associations(
      variants = snps, id = outcomes, proxies = if (proxies) 1 else 0,
      r2 = 0.8, align_alleles = 1, palindromes = 1,
      assocs_per_request = 16, max_ids_per_request = 5, timeout = 60
    ),
    error = function(e) e
  )
}

exact_assoc <- extract_assoc(requested_rsids, selected_ids, proxies = FALSE)
exact_error <- NA_character_
if (inherits(exact_assoc, "error")) {
  exact_error <- conditionMessage(exact_assoc)
  exact_assoc <- data.frame()
}

assoc_key <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(character())
  paste(get_col(df, c("id", "id.outcome")), get_col(df, c("rsid", "SNP", "variant")), sep = "||")
}

all_pairs <- as.vector(outer(selected_ids, requested_rsids, paste, sep = "||"))
missing_pairs <- setdiff(all_pairs, assoc_key(exact_assoc))
proxy_assoc <- data.frame()
proxy_error <- NA_character_
if (length(missing_pairs) > 0) {
  proxy_assoc <- extract_assoc(unique(sub("^.*\\|\\|", "", missing_pairs)), unique(sub("\\|\\|.*$", "", missing_pairs)), proxies = TRUE)
  if (inherits(proxy_assoc, "error")) {
    proxy_error <- conditionMessage(proxy_assoc)
    proxy_assoc <- data.frame()
  }
}

standardise_assoc <- function(df, association_type) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())
  data.frame(
    outcome_id = get_col(df, c("id", "id.outcome")),
    outcome_trait = get_col(df, c("trait", "outcome")),
    requested_rsid = get_col(df, c("target_snp", "target_snp.outcome", "rsid", "SNP", "variant")),
    matched_variant = get_col(df, c("rsid", "SNP", "variant", "proxy", "proxy.outcome")),
    association_type = association_type,
    proxy_r2 = get_col(df, c("proxy_r2", "proxy_r2.outcome", "r2")),
    effect_allele.outcome = get_col(df, c("ea", "effect_allele.outcome", "effect_allele")),
    other_allele.outcome = get_col(df, c("nea", "other_allele.outcome", "other_allele")),
    beta.outcome = get_num(get_col(df, c("beta", "beta.outcome"))),
    se.outcome = get_num(get_col(df, c("se", "se.outcome"))),
    p.outcome = get_num(get_col(df, c("p", "pval", "pval.outcome"))),
    eaf.outcome = get_col(df, c("eaf", "eaf.outcome")),
    stringsAsFactors = FALSE
  )
}

assoc_std <- rbind(
  standardise_assoc(exact_assoc, "exact"),
  standardise_assoc(proxy_assoc, "proxy_r2_gt_0.8_1000G_EUR")
)

replication_rows <- list()
row_i <- 1
for (outcome_id in selected_ids) {
  outcome_trait <- candidates$trait[match(outcome_id, candidates$id)]
  for (rsid in requested_rsids) {
    instrument_rows <- tyk2[tyk2$requested_rsid == rsid, , drop = FALSE]
    assoc_rows <- assoc_std[assoc_std$outcome_id == outcome_id & (assoc_std$requested_rsid == rsid | assoc_std$matched_variant == rsid), , drop = FALSE]
    if (nrow(assoc_rows) == 0) {
      for (j in seq_len(nrow(instrument_rows))) {
        replication_rows[[row_i]] <- data.frame(
          outcome_id = outcome_id, outcome_trait = outcome_trait, tissue = instrument_rows$tissue[j],
          requested_rsid = rsid, matched_variant = "unavailable", association_type = "unavailable",
          proxy_r2 = "unavailable", effect_allele.outcome = "unavailable", other_allele.outcome = "unavailable",
          beta.outcome = "unavailable", se.outcome = "unavailable", p.outcome = "unavailable", eaf.outcome = "unavailable",
          beta.exposure = instrument_rows$`beta.exposure`[j], se.exposure = instrument_rows$`se.exposure`[j],
          Wald_beta = "unavailable", Wald_se = "unavailable", Wald_p = "unavailable",
          OR = "unavailable", OR_lci = "unavailable", OR_uci = "unavailable",
          status = "unavailable", note = "Exact variant unavailable and LD proxy unavailable or not returned.",
          stringsAsFactors = FALSE
        )
        row_i <- row_i + 1
      }
      next
    }

    assoc_row <- assoc_rows[order(assoc_rows$association_type != "exact"), , drop = FALSE][1, , drop = FALSE]
    for (j in seq_len(nrow(instrument_rows))) {
      beta_exp <- get_num(instrument_rows$`beta.exposure`[j])
      se_exp <- get_num(instrument_rows$`se.exposure`[j])
      beta_out <- get_num(assoc_row$beta.outcome[1])
      se_out <- get_num(assoc_row$se.outcome[1])
      if (!is.na(beta_exp) && beta_exp != 0 && !is.na(beta_out) && !is.na(se_out)) {
        wald_beta <- beta_out / beta_exp
        wald_se <- abs(se_out / beta_exp)
        wald_p <- 2 * pnorm(abs(wald_beta / wald_se), lower.tail = FALSE)
        OR <- exp(wald_beta)
        OR_lci <- exp(wald_beta - 1.96 * wald_se)
        OR_uci <- exp(wald_beta + 1.96 * wald_se)
      } else {
        wald_beta <- wald_se <- wald_p <- OR <- OR_lci <- OR_uci <- NA_real_
      }
      replication_rows[[row_i]] <- data.frame(
        outcome_id = outcome_id, outcome_trait = outcome_trait, tissue = instrument_rows$tissue[j],
        requested_rsid = rsid, matched_variant = assoc_row$matched_variant[1],
        association_type = assoc_row$association_type[1],
        proxy_r2 = ifelse(is.na(assoc_row$proxy_r2[1]) || assoc_row$proxy_r2[1] == "", "unavailable", assoc_row$proxy_r2[1]),
        effect_allele.outcome = assoc_row$effect_allele.outcome[1],
        other_allele.outcome = assoc_row$other_allele.outcome[1],
        beta.outcome = beta_out, se.outcome = se_out, p.outcome = assoc_row$p.outcome[1],
        eaf.outcome = assoc_row$eaf.outcome[1], beta.exposure = beta_exp, se.exposure = se_exp,
        Wald_beta = wald_beta, Wald_se = wald_se, Wald_p = wald_p,
        OR = OR, OR_lci = OR_lci, OR_uci = OR_uci,
        status = "available",
        note = "OpenGWAS association-lookup Wald ratio uses external outcome association and the existing GTEx exposure estimate; single-SNP estimate only.",
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1
    }
  }
}

replication <- if (length(replication_rows) == 0) empty_replication() else do.call(rbind, replication_rows)
write.csv(replication, replication_path, row.names = FALSE)

available_n <- sum(replication$status == "available", na.rm = TRUE)
unavailable_n <- sum(replication$status == "unavailable", na.rm = TRUE)
status <- if (available_n > 0) "completed_with_available_results" else "unavailable"
extra <- c(
  "",
  paste0("- Requested rsIDs: ", paste(requested_rsids, collapse = ", ")),
  paste0("- PheWAS psoriasis rows at p <= 0.01: ", nrow(psoriasis_phewas)),
  paste0("- Selected non-FinnGen candidate outcomes: ", length(selected_ids)),
  paste0("- Available instrument-outcome rows: ", available_n),
  paste0("- Unavailable instrument-outcome rows: ", unavailable_n)
)
if (!is.na(exact_error)) extra <- c(extra, paste0("- Exact extraction error: ", exact_error))
if (!is.na(proxy_error)) extra <- c(extra, paste0("- Proxy extraction error: ", proxy_error))

write_status(
  status,
  "OpenGWAS association lookup used authenticated ieugwasr calls. Candidate discovery used PheWAS for the two TYK2 rsIDs, excluded identifiable FinnGen outcomes, prioritized European ancestry/case count/summary-stat availability where metadata were present, then attempted exact rsID extraction followed by LD proxy lookup at r2 > 0.8 in 1000G EUR. Available records provide directional concordance and are not assumed to be statistically independent replication.",
  extra
)
writeLines("OpenGWAS association-lookup script completed", file.path(logs_dir, "opengwas_replication.log"), useBytes = TRUE)
quit(status = 0)
