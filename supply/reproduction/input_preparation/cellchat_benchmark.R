#!/usr/bin/env Rscript

parse_args <- function(argv) {
  args <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- sub("^--", "", key)
    if (i == length(argv) || startsWith(argv[[i + 1]], "--")) {
      args[[key]] <- TRUE
      i <- i + 1
    } else {
      args[[key]] <- argv[[i + 1]]
      i <- i + 2
    }
  }
  args
}

need_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing required R package: ", pkg, call. = FALSE)
  }
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x) || x == "") y else x

safe_num <- function(x, default = NA_real_) {
  if (is.null(x) || is.na(x) || x == "") return(default)
  as.numeric(x)
}

load_seurat <- function(path) {
  obj <- readRDS(path)
  if (!inherits(obj, "Seurat")) {
    stop("Input is not a Seurat object: ", path)
  }
  obj
}

choose_label_col <- function(seu, requested = "auto") {
  meta <- seu@meta.data
  if (!identical(requested, "auto") && requested %in% colnames(meta)) return(requested)
  candidates <- c(
    "cell_type", "celltype", "cell.type", "cell_type_refined",
    "annotation", "annot", "seurat_clusters", "cluster", "CellType"
  )
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) > 0) return(hit[[1]])
  NA_character_
}

prepare_groups <- function(seu, label_col, dims = 1:20, resolution = 0.5, seed = 1) {
  set.seed(seed)
  if (!is.na(label_col) && label_col %in% colnames(seu@meta.data)) {
    group <- as.character(seu@meta.data[[label_col]])
    group[is.na(group) | group == ""] <- "Unknown"
    seu$cellchat_group <- paste0("Group_", make.names(group))
    return(seu)
  }

  Seurat::DefaultAssay(seu) <- "RNA"
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  seu <- Seurat::RunPCA(seu, verbose = FALSE)
  seu <- Seurat::FindNeighbors(seu, dims = dims, verbose = FALSE)
  seu <- Seurat::FindClusters(seu, resolution = resolution, verbose = FALSE)
  seu$cellchat_group <- paste0("Cluster_", as.character(seu$seurat_clusters))
  seu
}

downsample_cells <- function(seu, n_cells, seed) {
  if (identical(n_cells, "all") || is.na(n_cells)) return(seu)
  n_cells <- as.integer(n_cells)
  cells <- colnames(seu)
  if (length(cells) <= n_cells) return(seu)
  set.seed(seed)
  keep <- sample(cells, n_cells)
  subset(seu, cells = keep)
}

get_rna_data <- function(seu) {
  Seurat::DefaultAssay(seu) <- "RNA"
  data <- tryCatch(
    Seurat::GetAssayData(seu, assay = "RNA", layer = "data"),
    error = function(e) NULL
  )
  if (is.null(data)) {
    data <- tryCatch(
      Seurat::GetAssayData(seu, assay = "RNA", slot = "data"),
      error = function(e) NULL
    )
  }
  if (is.null(data) || nrow(data) == 0 || ncol(data) == 0) {
    stop("Cannot read RNA normalized data from Seurat object.", call. = FALSE)
  }
  data
}

make_cellchat <- function(seu) {
  data.input <- get_rna_data(seu)
  meta <- seu@meta.data
  meta$cellchat_group <- as.character(seu$cellchat_group)
  cellchat <- CellChat::createCellChat(
    object = data.input,
    meta = meta,
    group.by = "cellchat_group"
  )
  cellchat@DB <- CellChat::CellChatDB.human
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat
}

find_accel_fun <- function(pkg, explicit_fun = "") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing accelerated package: ", pkg)
  }
  if (nzchar(explicit_fun)) {
    if (grepl("::", explicit_fun, fixed = TRUE)) {
      parts <- strsplit(explicit_fun, "::", fixed = TRUE)[[1]]
      return(getExportedValue(parts[[1]], parts[[2]]))
    }
    return(getExportedValue(pkg, explicit_fun))
  }
  candidates <- c(
    "computeCommunProbAccelRcpp",
    "computeCommunProbAccel",
    "computeCommunProb_accel",
    "computeCommunProb_Rcpp",
    "computeCommunProb_cpp",
    "accelerated_computeCommunProb"
  )
  ns <- getNamespace(pkg)
  for (nm in candidates) {
    if (exists(nm, envir = ns, inherits = FALSE)) return(get(nm, envir = ns))
  }
  stop(
    "Cannot find accelerated function in ", pkg,
    ". Set CELLCHAT_ACCEL_FUN or pass --accel-fun."
  )
}

get_optional_export <- function(pkg, fun) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NULL)
  ns <- getNamespace(pkg)
  if (!exists(fun, envir = ns, inherits = FALSE)) return(NULL)
  get(fun, envir = ns)
}

run_engine <- function(cellchat, engine, min_cells, accel_pkg, accel_fun, ablation) {
  if (identical(engine, "baseline")) {
    cellchat <- CellChat::computeCommunProb(cellchat)
    cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells)
    cellchat <- CellChat::computeCommunProbPathway(cellchat)
    cellchat <- CellChat::aggregateNet(cellchat)
  } else if (identical(engine, "accelerated")) {
    Sys.setenv(CELLCHAT_ACCEL_ABLATION = ablation)
    options(CellChatAccelRcpp.ablation = ablation)
    use_accel_kernel <- !(ablation %in% c("no_accel_kernel", "baseline_kernel"))
    use_accel_pathway <- !(ablation %in% c("no_accel_pathway", "baseline_pathway"))
    use_accel_aggregate <- !(ablation %in% c("no_accel_aggregate", "baseline_aggregate"))

    if (use_accel_kernel) {
      fun <- find_accel_fun(accel_pkg, accel_fun)
      cellchat <- fun(cellchat)
    } else {
      cellchat <- CellChat::computeCommunProb(cellchat)
    }
    cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells)

    accel_pathway <- get_optional_export(accel_pkg, "computeCommunProbPathwayAccelRcpp")
    if (use_accel_pathway && !is.null(accel_pathway)) {
      cellchat <- accel_pathway(cellchat)
    } else {
      cellchat <- CellChat::computeCommunProbPathway(cellchat)
    }

    accel_aggregate <- get_optional_export(accel_pkg, "aggregateNetAccelRcpp")
    if (use_accel_aggregate && !is.null(accel_aggregate)) {
      cellchat <- accel_aggregate(cellchat)
    } else {
      cellchat <- CellChat::aggregateNet(cellchat)
    }
  } else {
    stop("Unsupported engine: ", engine)
  }
  cellchat
}

flatten_prob <- function(cellchat) {
  prob <- cellchat@net$prob
  if (is.null(prob)) return(numeric())
  as.numeric(prob)
}

compare_objects <- function(a, b) {
  x <- flatten_prob(a)
  y <- flatten_prob(b)
  n <- min(length(x), length(y))
  if (n == 0) {
    return(list(max_abs_prob_diff = NA_real_, pearson_prob = NA_real_))
  }
  x <- x[seq_len(n)]
  y <- y[seq_len(n)]
  list(
    max_abs_prob_diff = max(abs(x - y), na.rm = TRUE),
    pearson_prob = suppressWarnings(stats::cor(x, y, use = "pairwise.complete.obs"))
  )
}

metric_row <- function(args, engine, status, elapsed, cellchat = NULL, err = "", extra = list()) {
  n_lr <- NA_integer_
  mean_weight <- NA_real_
  max_weight <- NA_real_
  n_groups <- NA_integer_
  if (!is.null(cellchat)) {
    lr <- try(CellChat::subsetCommunication(cellchat), silent = TRUE)
    if (!inherits(lr, "try-error")) n_lr <- nrow(lr)
    w <- try(as.numeric(cellchat@net$weight), silent = TRUE)
    if (!inherits(w, "try-error")) {
      mean_weight <- mean(w, na.rm = TRUE)
      max_weight <- max(w, na.rm = TRUE)
    }
    n_groups <- tryCatch(length(levels(cellchat@idents)), error = function(e) NA_integer_)
  }
  data.frame(
    experiment_id = args[["experiment-id"]] %||% "",
    dataset_id = args[["dataset-id"]] %||% "",
    input_path = args[["input"]] %||% "",
    engine = engine,
    ablation = args[["ablation"]] %||% "none",
    n_cells_target = args[["n-cells"]] %||% "all",
    seed = args[["seed"]] %||% "",
    status = status,
    elapsed_sec = round(elapsed, 4),
    n_groups = n_groups,
    total_lr_pairs = n_lr,
    mean_weight = mean_weight,
    max_weight = max_weight,
    max_abs_prob_diff = extra$max_abs_prob_diff %||% NA_real_,
    pearson_prob = extra$pearson_prob %||% NA_real_,
    error = err,
    stringsAsFactors = FALSE
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  need_pkg("Seurat")
  need_pkg("CellChat")
  suppressPackageStartupMessages({
    library(Seurat)
    library(CellChat)
  })

  input <- args[["input"]] %||% stop("--input is required")
  out_dir <- args[["out-dir"]] %||% "/home/dzf/cellchat_acceleration/results/runs"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  exp_id <- args[["experiment-id"]] %||% paste0("manual_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  metrics_path <- file.path(out_dir, paste0(exp_id, ".metrics.tsv"))
  if (identical(args[["resume"]], "true") && file.exists(metrics_path)) {
    message("skip_existing=", metrics_path)
    return(invisible(NULL))
  }
  checkpoint_dir <- file.path(out_dir, "checkpoints", exp_id)
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

  seed <- as.integer(args[["seed"]] %||% 1)
  n_cells <- args[["n-cells"]] %||% "all"
  min_cells <- as.integer(args[["min-cells"]] %||% 10)
  engine <- args[["engine"]] %||% "both"
  accel_pkg <- Sys.getenv("CELLCHAT_ACCEL_PACKAGE", args[["accel-package"]] %||% "CellChatAccelRcpp")
  accel_fun <- Sys.getenv("CELLCHAT_ACCEL_FUN", args[["accel-fun"]] %||% "")
  ablation <- args[["ablation"]] %||% "full"
  label_col <- args[["label-col"]] %||% "auto"

  rows <- list()
  prepared_cp <- file.path(checkpoint_dir, "prepared_cellchat.rds")
  if (identical(args[["resume"]], "true") && file.exists(prepared_cp)) {
    message("resume_checkpoint=", prepared_cp)
    base_cellchat <- readRDS(prepared_cp)
  } else {
    set.seed(seed)
    seu <- load_seurat(input)
    seu <- downsample_cells(seu, n_cells, seed)
    chosen_label <- choose_label_col(seu, label_col)
    seu <- prepare_groups(seu, chosen_label, seed = seed)
    base_cellchat <- make_cellchat(seu)
    saveRDS(base_cellchat, prepared_cp)
  }

  engines <- if (engine == "both") c("baseline", "accelerated") else engine
  objects <- list()
  for (eng in engines) {
    engine_cp <- file.path(checkpoint_dir, paste0(eng, ".computed.rds"))
    t0 <- proc.time()[["elapsed"]]
    obj <- NULL
    err <- ""
    status <- "ok"
    if (identical(args[["resume"]], "true") && file.exists(engine_cp)) {
      message("resume_checkpoint=", engine_cp)
      obj <- readRDS(engine_cp)
      status <- "ok_resumed"
    } else {
      obj <- tryCatch(
        run_engine(base_cellchat, eng, min_cells, accel_pkg, accel_fun, ablation),
        error = function(e) {
          status <<- "error"
          err <<- conditionMessage(e)
          NULL
        }
      )
      if (!is.null(obj)) saveRDS(obj, engine_cp)
    }
    elapsed <- proc.time()[["elapsed"]] - t0
    objects[[eng]] <- obj
    rows[[length(rows) + 1]] <- metric_row(args, eng, status, elapsed, obj, err)
    if (!is.null(obj) && identical(args[["save-objects"]], "true")) {
      saveRDS(obj, file.path(out_dir, paste0(exp_id, ".", eng, ".cellchat.rds")))
    }
  }

  if (!is.null(objects$baseline) && !is.null(objects$accelerated)) {
    cmp <- compare_objects(objects$baseline, objects$accelerated)
    rows[[length(rows) + 1]] <- metric_row(
      args, "comparison", "ok", 0, NULL, "", cmp
    )
  }

  metrics <- do.call(rbind, rows)
  utils::write.table(metrics, metrics_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("metrics=", metrics_path)
}

main()
