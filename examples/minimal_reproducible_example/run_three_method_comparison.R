#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(CellChatESS)
})

all_args <- commandArgs(trailingOnly = FALSE)
script_file <- normalizePath(sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[[1]]))
example_dir <- dirname(script_file)
output_dir <- file.path(example_dir, "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
prepared <- readRDS(file.path(example_dir, "data", "pbmc3k_prepared_cellchat.rds"))

run_original <- function(object) {
  object <- CellChat::computeCommunProb(object, nboot = 100, seed.use = 1L)
  object <- CellChat::filterCommunication(object, min.cells = 10)
  object <- CellChat::computeCommunProbPathway(object)
  CellChat::aggregateNet(object)
}

run_rcpp <- function(object, algorithm) {
  object <- computeCommunProbESS(object, nboot = 100, seed.use = 1L, algorithm = algorithm)
  object <- CellChat::filterCommunication(object, min.cells = 10)
  object <- computeCommunProbPathwayESS(object)
  aggregateNetESS(object)
}

original_time <- system.time(original <- run_original(prepared))[["elapsed"]]
dense_time <- system.time(dense <- run_rcpp(prepared, "dense"))[["elapsed"]]
ess_time <- system.time(ess <- run_rcpp(prepared, "ess"))[["elapsed"]]

summarize <- function(name, elapsed, object) {
  data.frame(
    method = name,
    elapsed_sec = elapsed,
    prob_max_abs_diff = max(abs(object@net$prob - original@net$prob), na.rm = TRUE),
    pval_max_abs_diff = max(abs(object@net$pval - original@net$pval), na.rm = TRUE),
    pathway_max_abs_diff = max(abs(object@netP$prob - original@netP$prob), na.rm = TRUE),
    weight_max_abs_diff = max(abs(object@net$weight - original@net$weight), na.rm = TRUE),
    count_identical = identical(object@net$count, original@net$count)
  )
}

metrics <- rbind(
  summarize("Original CellChat R", original_time, original),
  summarize("Dense Rcpp", dense_time, dense),
  summarize("ESS", ess_time, ess)
)
write.table(metrics, file.path(output_dir, "three_method_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(metrics)
stopifnot(
  max(metrics$prob_max_abs_diff) <= 1e-12,
  max(metrics$pval_max_abs_diff) == 0,
  max(metrics$pathway_max_abs_diff) <= 1e-12,
  max(metrics$weight_max_abs_diff) <= 1e-12,
  all(metrics$count_identical)
)
