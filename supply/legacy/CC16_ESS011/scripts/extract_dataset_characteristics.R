args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop('usage: extract_dataset_characteristics.R MANIFEST OUT_TSV SESSION_TXT')
manifest_file <- args[[1L]]; out_file <- args[[2L]]; session_file <- args[[3L]]
archived_lib <- Sys.getenv("CELLCHAT_MAIN_R_LIB", "environment/main/library")
.libPaths(c(archived_lib, .libPaths()))
suppressPackageStartupMessages({library(CellChat); library(CellChatESS); library(Matrix)})
manifest <- read.delim(manifest_file, check.names=FALSE, stringsAsFactors=FALSE)
manifest <- manifest[grepl('engine-both__ablation-full$', manifest$experiment_id), , drop=FALSE]
if (nrow(manifest) != 216L) stop('Expected 216 full-ablation prepared objects; found ', nrow(manifest))
if (any(!file.exists(manifest$prepared))) stop('One or more prepared object paths are missing')
rows <- vector('list', nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  x <- manifest[i, ]
  obj <- readRDS(x$prepared)
  expr <- obj@data.signaling
  groups <- obj@idents
  group_sizes <- as.numeric(table(groups))
  k <- nlevels(groups)
  nlr <- nrow(obj@LR$LRsig)
  total <- as.double(nlr) * as.double(k) * as.double(k)
  nnz <- if (inherits(expr, 'Matrix')) Matrix::nnzero(expr) else sum(expr != 0, na.rm=TRUE)
  active_fraction <- NA_real_
  checked <- tryCatch(
    CellChatESS::computeCommunProbESS(obj, type='triMean', raw.use=TRUE, population.size=FALSE,
      nboot=1L, seed.use=1L, Kh=0.5, n=1, algorithm='ess'),
    error=function(e) e)
  if (inherits(checked, 'error')) stop('ESS support count failed for ', x$experiment_id, ': ', conditionMessage(checked))
  metrics <- checked@options$CellChatESS
  active_fraction <- as.numeric(metrics$active_fraction)
  if (length(active_fraction) != 1L || !is.finite(active_fraction)) stop('Missing active_fraction for ', x$experiment_id)
  rows[[i]] <- data.frame(
    experiment_id=x$experiment_id, dataset_id=x$dataset_id, target_cells=x$n_cells_target,
    downsample_seed=x$seed, actual_cells=ncol(expr), cell_groups=k,
    min_group_cells=min(group_sizes), median_group_cells=median(group_sizes), max_group_cells=max(group_sizes),
    evaluable_lr_pairs=nlr, matrix_genes=nrow(expr), matrix_nonzero=nnz,
    matrix_sparsity=1 - nnz/(as.double(nrow(expr))*as.double(ncol(expr))),
    lr_group_pairs_total=total, lr_group_pairs_after_pruning=round(total*active_fraction),
    group_pairs_retained_fraction=active_fraction,
    data_signaling_size_mib=as.numeric(object.size(expr))/1024^2,
    stringsAsFactors=FALSE
  )
  if (i %% 12L == 0L) cat('processed ', i, '/', nrow(manifest), ': ', x$experiment_id, '\n', sep='')
  rm(obj, expr, checked); gc(verbose=FALSE)
}
out <- do.call(rbind, rows)
write.table(out, out_file, sep='\t', row.names=FALSE, quote=TRUE)
# Replicate-aware summaries at dataset by target scale for supplementary reporting.
key <- interaction(out$dataset_id, out$target_cells, drop=TRUE, lex.order=TRUE)
summ <- do.call(rbind, lapply(split(out, key), function(d) data.frame(
  dataset_id=d$dataset_id[[1]], target_cells=d$target_cells[[1]], n_replicates=nrow(d),
  actual_cells_median=median(d$actual_cells), cell_groups_median=median(d$cell_groups),
  evaluable_lr_pairs_median=median(d$evaluable_lr_pairs), matrix_sparsity_median=median(d$matrix_sparsity),
  group_pairs_retained_fraction_median=median(d$group_pairs_retained_fraction),
  min_group_cells_median=median(d$min_group_cells), median_group_cells_median=median(d$median_group_cells),
  max_group_cells_median=median(d$max_group_cells), data_signaling_size_mib_median=median(d$data_signaling_size_mib)
)))
write.table(summ, sub('\\.tsv$', '_dataset_scale_summary.tsv', out_file), sep='\t', row.names=FALSE, quote=TRUE)
writeLines(capture.output(sessionInfo()), session_file)
cat('R=', as.character(getRversion()), ' CellChat=', as.character(packageVersion('CellChat')),
    ' CellChatESS=', as.character(packageVersion('CellChatESS')), '\n')
cat('Wrote ', nrow(out), ' per-experiment records and ', nrow(summ), ' dataset-scale summaries.\n', sep='')
