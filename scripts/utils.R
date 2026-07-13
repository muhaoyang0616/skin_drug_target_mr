bootstrap_project <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
    root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
    setwd(root)
    return(invisible(root))
  }
  invisible(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
}

bootstrap_project()

required_packages <- c("data.table", "yaml")
missing_required <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_required, collapse = ", "),
    ". Run scripts/00_setup.R first.",
    call. = FALSE
  )
}

check_file <- function(path, label = path) {
  if (is.null(path) || is.na(path) || identical(path, "")) {
    stop("Missing configured path for ", label, ".", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(
      "Required input not found for ", label, ": ", path, "\n",
      "Place the real file at this path or update analysis_config.yaml.",
      call. = FALSE
    )
  }
  invisible(path)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

read_config <- function(path = "analysis_config.yaml") {
  check_file(path, "analysis config")
  yaml::read_yaml(path)
}

read_csv_required <- function(path, label = path) {
  check_file(path, label)
  data.table::fread(path, data.table = FALSE)
}

write_csv <- function(x, path) {
  ensure_dir(dirname(path))
  data.table::fwrite(as.data.frame(x), path)
  invisible(path)
}

write_unavailable <- function(path, step, reason) {
  out <- data.frame(
    step = step,
    status = "unavailable",
    reason = reason,
    stringsAsFactors = FALSE
  )
  write_csv(out, path)
}

normalize_colname <- function(x) {
  tolower(gsub("[^[:alnum:]]", "", x))
}

find_col <- function(cols, candidates, label, required = TRUE) {
  idx <- match(normalize_colname(candidates), normalize_colname(cols))
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) {
    return(cols[idx[[1]]])
  }
  if (required) {
    stop(
      "Could not find required column for ", label, ". Tried: ",
      paste(candidates, collapse = ", "), ". Available columns: ",
      paste(cols, collapse = ", "),
      call. = FALSE
    )
  }
  NA_character_
}

coalesce_col <- function(df, candidates, default = NA) {
  col <- find_col(names(df), candidates, label = paste(candidates, collapse = "/"), required = FALSE)
  if (is.na(col)) {
    return(rep(default, nrow(df)))
  }
  df[[col]]
}

standardize_p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[x < 0 | x > 1] <- NA_real_
  x
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

empty_df <- function(cols) {
  as.data.frame(setNames(replicate(length(cols), logical(0), simplify = FALSE), cols))
}

status_message <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

append_log <- function(path, ...) {
  ensure_dir(dirname(path))
  line <- paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(line, "\n", file = path, append = TRUE)
  invisible(line)
}

download_if_missing <- function(url, path, label = path) {
  ensure_dir(dirname(path))
  if (file.exists(path) && file.info(path)$size > 0) {
    status_message(label, " already exists: ", path)
    return(invisible(path))
  }
  if (is.null(url) || is.na(url) || !nzchar(url)) {
    stop("Missing download URL for ", label, ".", call. = FALSE)
  }
  status_message("Downloading ", label, " from ", url)
  tmp <- paste0(path, ".download")
  if (file.exists(tmp)) unlink(tmp)
  curl <- Sys.which("curl.exe")
  if (nzchar(curl)) {
    args <- c(
      "--location",
      "--fail",
      "--retry", "20",
      "--retry-delay", "10",
      "--retry-all-errors",
      "--continue-at", "-",
      "--http1.1",
      "--speed-time", "300",
      "--speed-limit", "1024",
      "--ssl-no-revoke",
      "--output", tmp,
      url
    )
    code <- system2(curl, args = args)
    if (identical(code, 0L) && file.exists(tmp) && file.info(tmp)$size > 0) {
      if (file.exists(path)) unlink(path)
      file.rename(tmp, path)
      return(invisible(path))
    }
  }
  old_timeout <- getOption("timeout")
  options(timeout = max(3600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  utils::download.file(url, tmp, mode = "wb", quiet = FALSE, method = "libcurl")
  if (!file.exists(tmp) || file.info(tmp)$size == 0) {
    stop("Download failed or produced an empty file for ", label, ": ", url, call. = FALSE)
  }
  if (file.exists(path)) unlink(path)
  file.rename(tmp, path)
  invisible(path)
}

strip_ensembl_version <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

normalize_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  toupper(x)
}

make_variant_key <- function(chr, pos, ref, alt) {
  paste(normalize_chr(chr), as.integer(pos), toupper(as.character(ref)), toupper(as.character(alt)), sep = ":")
}
