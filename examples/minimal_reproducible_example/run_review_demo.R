#!/usr/bin/env Rscript

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]))
root <- dirname(script_file)
local_lib <- file.path(root, ".example_lib")
package_source <- normalizePath(file.path(root, "..", ".."))
data_dir <- file.path(root, "data")
results_dir <- file.path(root, "results")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("Rcpp", "CellChat")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop(
    "Install prerequisite package(s) before running: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

stopifnot(file.exists(file.path(package_source, "DESCRIPTION")))
message("[1/4] Installing repository CellChatESS source into .example_lib/")
r_bin <- file.path(R.home("bin"), "R")
install_status <- system2(
  r_bin,
  args = c(
    "CMD", "INSTALL", "--preclean",
    paste0("--library=", local_lib),
    package_source
  ),
  stdout = TRUE,
  stderr = TRUE
)
status_code <- attr(install_status, "status")
if (!is.null(status_code) && status_code != 0L) {
  cat(paste(install_status, collapse = "\n"), "\n")
  stop("Bundled CellChatESS installation failed.", call. = FALSE)
}

.libPaths(c(local_lib, .libPaths()))
Sys.setenv(
  OMP_NUM_THREADS = "4",
  OMP_DYNAMIC = "FALSE",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1"
)
suppressPackageStartupMessages({
  library(CellChat)
  library(CellChatESS)
})
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")

prepared_path <- file.path(data_dir, "pbmc3k_prepared_cellchat.rds")
reference_path <- file.path(data_dir, "pbmc3k_original_reference.rds")
stopifnot(file.exists(prepared_path), file.exists(reference_path))
reference <- readRDS(reference_path)

run_ess <- function(threads) {
  object <- readRDS(prepared_path)
  object@idents <- droplevels(object@idents)
  start <- proc.time()[["elapsed"]]
  object <- CellChatESS::computeCommunProbESS(
    object,
    nboot = 100,
    seed.use = 1L,
    threads = threads,
    algorithm = "ess"
  )
  object <- CellChat::filterCommunication(object, min.cells = 10)
  object <- CellChatESS::computeCommunProbPathwayESS(object)
  object <- CellChatESS::aggregateNetESS(object)
  list(
    object = object,
    elapsed_sec = proc.time()[["elapsed"]] - start
  )
}

message("[2/4] Running the true serial path (threads = 1L)")
one <- run_ess(1L)
message("[3/4] Running the conditional OpenMP path (threads = 4L)")
four <- run_ess(4L)

safe_max_diff <- function(x, y) {
  if (!identical(dim(x), dim(y))) return(Inf)
  max(abs(x - y), na.rm = TRUE)
}
summarize_run <- function(label, requested, x) {
  object <- x$object
  data.frame(
    run = label,
    threads_requested = requested,
    threads_used = object@options$CellChatESS$threads_used,
    elapsed_sec = x$elapsed_sec,
    prob_max_abs_diff_vs_original =
      safe_max_diff(object@net$prob, reference@net$prob),
    pval_max_abs_diff_vs_original =
      safe_max_diff(object@net$pval, reference@net$pval),
    pathway_max_abs_diff_vs_original =
      safe_max_diff(object@netP$prob, reference@netP$prob),
    weight_max_abs_diff_vs_original =
      safe_max_diff(object@net$weight, reference@net$weight),
    count_identical_vs_original =
      identical(object@net$count, reference@net$count),
    stringsAsFactors = FALSE
  )
}

metrics <- rbind(
  summarize_run("serial_1t", 1L, one),
  summarize_run("openmp_4t_requested", 4L, four)
)
checks <- data.frame(
  prob_exact_1t_vs_4t =
    identical(one$object@net$prob, four$object@net$prob),
  pval_exact_1t_vs_4t =
    identical(one$object@net$pval, four$object@net$pval),
  pathway_exact_1t_vs_4t =
    identical(one$object@netP$prob, four$object@netP$prob),
  weight_exact_1t_vs_4t =
    identical(one$object@net$weight, four$object@net$weight),
  count_exact_1t_vs_4t =
    identical(one$object@net$count, four$object@net$count),
  serial_path_reported_one_thread =
    identical(one$object@options$CellChatESS$threads_used, 1L),
  requested_four_threads_used =
    identical(four$object@options$CellChatESS$threads_used, 4L),
  speedup_1t_over_4t = one$elapsed_sec / four$elapsed_sec,
  stringsAsFactors = FALSE
)

write.table(
  metrics,
  file.path(results_dir, "reviewer_demo_metrics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  checks,
  file.path(results_dir, "reviewer_demo_checks.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "sessionInfo.txt")
)

message("[4/4] Validating exactness and original-CellChat agreement")
print(metrics, row.names = FALSE)
print(checks, row.names = FALSE)

stopifnot(
  all(metrics$prob_max_abs_diff_vs_original <= 1e-12),
  all(metrics$pval_max_abs_diff_vs_original == 0),
  all(metrics$pathway_max_abs_diff_vs_original <= 1e-12),
  all(metrics$weight_max_abs_diff_vs_original <= 1e-12),
  all(metrics$count_identical_vs_original),
  checks$prob_exact_1t_vs_4t,
  checks$pval_exact_1t_vs_4t,
  checks$pathway_exact_1t_vs_4t,
  checks$weight_exact_1t_vs_4t,
  checks$count_exact_1t_vs_4t,
  checks$serial_path_reported_one_thread
)

if (!checks$requested_four_threads_used) {
  message(
    "NOTE: this compiler has no active OpenMP support; the four-thread ",
    "request safely fell back to one thread."
  )
}
cat(
  "\nPASS: serial and OpenMP results are identical; ",
  "both match original CellChat.\n",
  sep = ""
)
