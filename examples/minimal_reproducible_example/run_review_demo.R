#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
root <- dirname(script_path)
package_source <- normalizePath(file.path(root, "..", ".."))
local_lib <- file.path(root, ".example_lib")
results_dir <- file.path(root, "results")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

required <- c("Rcpp", "Matrix", "CellChat")
missing <- required[
  !vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing)) {
  stop("Install prerequisite package(s) before running: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

Sys.setenv(
  OMP_NUM_THREADS = "4", OMP_DYNAMIC = "FALSE",
  OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1"
)
message("[1/5] Installing repository source into the example's private library")
install_output <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", "--preclean", shQuote(paste0("--library=", local_lib)),
    shQuote(package_source)),
  stdout = TRUE, stderr = TRUE
)
install_status <- attr(install_output, "status")
if (!is.null(install_status) && install_status != 0L) {
  cat(paste(install_output, collapse = "\n"), "\n")
  stop("CellChatESS installation failed.", call. = FALSE)
}
.libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages({
  library(CellChat)
  library(CellChatESS)
})
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")

prepared_path <- file.path(root, "data", "pbmc3k_prepared_cellchat.rds")
parameters <- list(
  type = "triMean", raw.use = TRUE, population.size = FALSE,
  nboot = 100L, seed.use = 1L, Kh = 0.5, n = 1
)
run_method <- function(method, threads = 1L) {
  object <- readRDS(prepared_path)
  object@idents <- droplevels(object@idents)
  start <- proc.time()[["elapsed"]]
  if (method == "original") {
    object <- do.call(CellChat::computeCommunProb, c(list(object = object), parameters))
  } else {
    object <- do.call(CellChatESS::computeCommunProbESS,
                      c(list(object = object), parameters,
                        list(threads = threads, algorithm = "ess")))
  }
  object <- CellChat::filterCommunication(object, min.cells = 10)
  if (method == "original") {
    object <- CellChat::computeCommunProbPathway(object)
    object <- CellChat::aggregateNet(object)
  } else {
    object <- CellChatESS::computeCommunProbPathwayESS(object)
    object <- CellChatESS::aggregateNetESS(object)
  }
  list(object = object, elapsed_sec = proc.time()[["elapsed"]] - start)
}

message("[2/5] Computing the original CellChat reference in the installed environment")
reference <- run_method("original")
message("[3/5] Computing CellChatESS with one thread")
one <- run_method("ess", 1L)
message("[4/5] Computing CellChatESS with four threads requested")
four <- run_method("ess", 4L)

max_difference <- function(x, y) {
  if (!identical(dim(x), dim(y)) || !identical(dimnames(x), dimnames(y))) return(Inf)
  if (!length(x)) return(0)
  max(abs(x - y))
}
summarize_run <- function(label, requested, result) {
  object <- result$object
  original <- reference$object
  data.frame(
    run = label,
    CellChat_version = as.character(packageVersion("CellChat")),
    CellChatESS_version = as.character(packageVersion("CellChatESS")),
    threads_requested = requested,
    threads_used = object@options$CellChatESS$threads_used,
    elapsed_sec = result$elapsed_sec,
    prob_max_abs_diff_vs_original = max_difference(object@net$prob, original@net$prob),
    pval_identical_vs_original = identical(object@net$pval, original@net$pval),
    pathway_max_abs_diff_vs_original = max_difference(object@netP$prob, original@netP$prob),
    weight_max_abs_diff_vs_original = max_difference(object@net$weight, original@net$weight),
    count_identical_vs_original = identical(object@net$count, original@net$count)
  )
}
metrics <- rbind(summarize_run("serial_1t", 1L, one),
                 summarize_run("four_threads_requested", 4L, four))
checks <- data.frame(
  prob_exact_1t_vs_4t = identical(one$object@net$prob, four$object@net$prob),
  pval_exact_1t_vs_4t = identical(one$object@net$pval, four$object@net$pval),
  pathway_exact_1t_vs_4t = identical(one$object@netP$prob, four$object@netP$prob),
  weight_exact_1t_vs_4t = identical(one$object@net$weight, four$object@net$weight),
  count_exact_1t_vs_4t = identical(one$object@net$count, four$object@net$count),
  serial_path_used_one_thread = identical(one$object@options$CellChatESS$threads_used, 1L),
  four_threads_were_used = identical(four$object@options$CellChatESS$threads_used, 4L)
)
write.table(metrics, file.path(results_dir, "demo_metrics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(checks, file.path(results_dir, "demo_checks.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo.txt"))
writeLines(capture.output(dput(parameters)), file.path(results_dir, "parameters.R"))
message("[5/5] Checking original-CellChat agreement and thread agreement")
print(metrics, row.names = FALSE)
print(checks, row.names = FALSE)
stopifnot(
  all(metrics$prob_max_abs_diff_vs_original <= 1e-12),
  all(metrics$pval_identical_vs_original),
  all(metrics$pathway_max_abs_diff_vs_original <= 1e-12),
  all(metrics$weight_max_abs_diff_vs_original <= 1e-12),
  all(metrics$count_identical_vs_original),
  checks$prob_exact_1t_vs_4t, checks$pval_exact_1t_vs_4t,
  checks$pathway_exact_1t_vs_4t, checks$weight_exact_1t_vs_4t,
  checks$count_exact_1t_vs_4t, checks$serial_path_used_one_thread
)
if (!checks$four_threads_were_used) {
  message("Four threads were requested. This build used the serial fallback.")
}
cat("\nPASS: CellChatESS agrees with the original CellChat reference and its one- and four-thread requests give identical outputs.\n")
