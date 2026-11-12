#!/usr/bin/env Rscript

parse_args <- function(argv) {
  args <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- sub("^--", "", key)
    args[[key]] <- argv[[i + 1]]
    i <- i + 2
  }
  args
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
grid_path <- args[["grid"]]
task_id <- as.integer(args[["task-id"]])
root <- args[["root"]]
if (is.null(root)) root <- "/home/dzf/cellchat_acceleration"
if (is.na(task_id) || task_id < 1) stop("--task-id must be a positive integer")

grid <- utils::read.csv(grid_path, stringsAsFactors = FALSE)
if (task_id > nrow(grid)) stop("task_id exceeds grid size")
row <- grid[task_id, ]

rscript <- file.path(R.home("bin"), "Rscript")
bench <- file.path(root, "code", "cellchat_benchmark.R")
out_dir <- file.path(root, "results", "runs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cmd_args <- c(
  bench,
  "--input", row$input_path,
  "--out-dir", out_dir,
  "--experiment-id", row$experiment_id,
  "--dataset-id", row$dataset_id,
  "--engine", row$engine,
  "--ablation", row$ablation,
  "--n-cells", as.character(row$n_cells),
  "--seed", as.character(row$seed),
  "--label-col", row$label_col,
  "--resume", "true"
)

status <- system2(rscript, cmd_args)
quit(save = "no", status = status)

