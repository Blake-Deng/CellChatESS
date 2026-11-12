#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- sub("^--", "", x[[i]])
    if (!startsWith(x[[i]], "--") || i == length(x)) stop("Arguments must be --key value")
    out[[key]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}
a <- parse_args(args)
need <- function(k) {
  if (is.null(a[[k]]) || !nzchar(a[[k]])) stop("Missing --", k, call. = FALSE)
  a[[k]]
}

method <- need("method")
if (!method %in% c("original", "ess")) stop("--method must be original or ess", call. = FALSE)
experiment_id <- need("id")
input_path <- need("input")
reference_path <- need("reference")
out_path <- need("out")
result_path <- need("result")
session_path <- sub("\\.tsv$", ".sessionInfo.txt", out_path)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

archived_lib <- Sys.getenv("CELLCHAT_MAIN_R_LIB", "environment/main/library")
.libPaths(c(archived_lib, .libPaths()))
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1"
)
suppressPackageStartupMessages({
  library(CellChat)
  library(CellChatESS)
})
stopifnot(as.character(packageVersion("CellChat")) == "1.6.1")
stopifnot(as.character(packageVersion("CellChatESS")) == "0.1.1")
writeLines(capture.output(sessionInfo()), session_path)

obj <- readRDS(input_path)
ref <- readRDS(reference_path)
pars <- ref@options$parameter
if (is.null(pars)) stop("Reference object has no parameter record", call. = FALSE)

type <- if (is.null(pars$type.mean)) "triMean" else as.character(pars$type.mean)
raw.use <- if (is.null(pars$raw.use)) TRUE else isTRUE(pars$raw.use)
population.size <- if (is.null(pars$population.size)) FALSE else isTRUE(pars$population.size)
nboot <- if (is.null(pars$nboot)) 100L else as.integer(pars$nboot)
seed.use <- if (is.null(pars$seed.use)) 1L else as.integer(pars$seed.use)
Kh <- if (is.null(pars$Kh)) 0.5 else as.numeric(pars$Kh)
n <- if (is.null(pars$n)) 1 else as.numeric(pars$n)

started <- proc.time()[["elapsed"]]
if (identical(method, "original")) {
  out <- CellChat::computeCommunProb(
    obj, type = type, raw.use = raw.use, population.size = population.size,
    nboot = nboot, seed.use = seed.use, Kh = Kh, n = n
  )
} else {
  out <- CellChatESS::computeCommunProbESS(
    obj, type = type, raw.use = raw.use, population.size = population.size,
    nboot = nboot, seed.use = seed.use, Kh = Kh, n = n, algorithm = "ess"
  )
}
elapsed <- unname(proc.time()[["elapsed"]] - started)

prob <- out@net$prob
pval <- out@net$pval
ref_prob <- ref@net$prob
ref_pval <- ref@net$pval
if (!identical(dim(prob), dim(ref_prob)) || !identical(dim(pval), dim(ref_pval))) {
  stop("Output dimensions differ from reference", call. = FALSE)
}
prob_diff <- abs(prob - ref_prob)
pval_diff <- abs(pval - ref_pval)
max_prob_diff <- max(prob_diff, na.rm = TRUE)
max_pval_diff <- max(pval_diff, na.rm = TRUE)
prob_exact <- sum(prob != ref_prob, na.rm = TRUE) == 0L
pval_exact <- sum(pval != ref_pval, na.rm = TRUE) == 0L

saveRDS(list(prob = prob, pval = pval), result_path, compress = FALSE)
internal <- tryCatch(as.numeric(out@options$run.time), error = function(e) NA_real_)
row <- data.frame(
  experiment_id = experiment_id,
  method = method,
  status = "ok",
  inference_call_sec = elapsed,
  package_internal_run_time_sec = internal,
  cells = ncol(obj@data.signaling),
  groups = length(levels(obj@idents)),
  n_LR = nrow(obj@LR$LRsig),
  nboot = nboot,
  seed_use = seed.use,
  max_abs_prob_diff_vs_reference = max_prob_diff,
  max_abs_pvalue_diff_vs_reference = max_pval_diff,
  probability_elementwise_exact_vs_reference = prob_exact,
  pvalue_elementwise_exact_vs_reference = pval_exact,
  probability_payload_bytes = as.numeric(object.size(prob)),
  pvalue_payload_bytes = as.numeric(object.size(pval)),
  CellChat_version = as.character(packageVersion("CellChat")),
  CellChatESS_version = as.character(packageVersion("CellChatESS")),
  threads_requested = 1L,
  input_rds = normalizePath(input_path, winslash = "/", mustWork = FALSE),
  reference_rds = normalizePath(reference_path, winslash = "/", mustWork = FALSE),
  stringsAsFactors = FALSE
)
write.table(row, out_path, sep = "\t", row.names = FALSE, quote = TRUE)
cat(sprintf("%s\t%s\t%.6f\n", experiment_id, method, elapsed))
