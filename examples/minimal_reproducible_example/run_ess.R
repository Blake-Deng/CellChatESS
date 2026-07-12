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
reference <- readRDS(file.path(example_dir, "data", "pbmc3k_original_reference.rds"))

elapsed <- system.time({
  result <- computeCommunProbESS(prepared, nboot = 100, seed.use = 1L, algorithm = "ess")
  result <- computeCommunProbPathwayESS(result)
  result <- aggregateNetESS(result)
})[["elapsed"]]

metrics <- data.frame(
  method = "ESS",
  elapsed_sec = elapsed,
  prob_max_abs_diff_vs_original = max(abs(result@net$prob - reference@net$prob), na.rm = TRUE),
  pval_max_abs_diff_vs_original = max(abs(result@net$pval - reference@net$pval), na.rm = TRUE),
  pathway_max_abs_diff_vs_original = max(abs(result@netP$prob - reference@netP$prob), na.rm = TRUE),
  weight_max_abs_diff_vs_original = max(abs(result@net$weight - reference@net$weight), na.rm = TRUE),
  count_identical = identical(result@net$count, reference@net$count),
  n_genes_full = result@options$CellChatESS$n_genes_full,
  n_genes_kernel = result@options$CellChatESS$n_genes_kernel
)
write.table(metrics, file.path(output_dir, "ess_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(metrics)
stopifnot(
  metrics$prob_max_abs_diff_vs_original <= 1e-12,
  metrics$pval_max_abs_diff_vs_original == 0,
  metrics$pathway_max_abs_diff_vs_original <= 1e-12,
  metrics$weight_max_abs_diff_vs_original <= 1e-12,
  metrics$count_identical
)
