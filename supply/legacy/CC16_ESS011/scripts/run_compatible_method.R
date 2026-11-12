args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5L) {
  stop('usage: run_cellchat22_multidataset_method.R METHOD INPUT_RDS REFERENCE_RDS_OR_DASH OUT_TSV OUT_RDS')
}
method <- args[[1L]]
input_path <- args[[2L]]
reference_path <- args[[3L]]
out_path <- args[[4L]]
result_path <- args[[5L]]
if (!method %in% c('original', 'autozyme', 'ess')) stop('METHOD must be original, autozyme or ess')

cellchat22_lib <- Sys.getenv("CELLCHAT_COMPAT_R_LIB", "environment/compatible/library")
ess_lib <- Sys.getenv("CELLCHAT_ESS_R_LIB", "environment/compatible/library")
archived_lib <- Sys.getenv("CELLCHAT_MAIN_R_LIB", "environment/main/library")
.libPaths(c(ess_lib, cellchat22_lib, archived_lib, .libPaths()))
Sys.setenv(AUTOZYME_THREADS='1', ZYME_THREADS='1', OMP_NUM_THREADS='1',
           OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
           VECLIB_MAXIMUM_THREADS='1', RCPP_PARALLEL_NUM_THREADS='1')
suppressPackageStartupMessages(library(CellChat))
if (method == 'ess') suppressPackageStartupMessages(library(CellChatESS))
if (method == 'autozyme') suppressPackageStartupMessages(library(autozyme))
stopifnot(as.character(packageVersion('CellChat')) == '2.2.0.9001')
if (method == 'ess') stopifnot(as.character(packageVersion('CellChatESS')) == '0.2.2')
if (method == 'autozyme') stopifnot(as.character(packageVersion('autozyme')) == '0.3.1')

obj <- readRDS(input_path)
obj@idents <- droplevels(obj@idents)
if (!identical(obj@options$datatype, 'RNA')) stop('input CellChat object must use RNA datatype')
if (method == 'autozyme') {
  autozyme::set_threads(1L)
  autozyme::activate('cellchat')
  if (!identical(unname(autozyme::status()[['cellchat']]), 'active')) stop('AutoZyme patch did not activate')
}
future::plan('sequential')
gc()
started <- proc.time()[['elapsed']]
result <- if (method == 'ess') {
  CellChatESS::computeCommunProbESS(obj, type='triMean', raw.use=TRUE,
    population.size=FALSE, nboot=100L, seed.use=1L, Kh=0.5, n=1,
    algorithm='ess')
} else {
  CellChat::computeCommunProb(obj, type='triMean', raw.use=TRUE,
    population.size=FALSE, nboot=100L, seed.use=1L, Kh=0.5, n=1)
}
elapsed <- proc.time()[['elapsed']] - started
saveRDS(list(prob=result@net$prob, pval=result@net$pval), result_path, compress=FALSE)

prob_max_abs_diff <- NA_real_
pval_max_abs_diff <- NA_real_
prob_diff_count <- NA_integer_
pval_diff_count <- NA_integer_
prob_exact <- NA
pval_exact <- NA
if (reference_path != '-' && file.exists(reference_path)) {
  reference <- readRDS(reference_path)
  if (is.null(reference$prob) || is.null(reference$pval)) stop('reference RDS must contain prob and pval arrays')
  if (!identical(dim(result@net$prob), dim(reference$prob)) || !identical(dim(result@net$pval), dim(reference$pval))) {
    stop('result and reference tensor dimensions differ')
  }
  prob_diff <- result@net$prob - reference$prob
  pval_diff <- result@net$pval - reference$pval
  prob_max_abs_diff <- max(abs(prob_diff), na.rm=TRUE)
  pval_max_abs_diff <- max(abs(pval_diff), na.rm=TRUE)
  prob_diff_count <- sum(prob_diff != 0, na.rm=TRUE)
  pval_diff_count <- sum(pval_diff != 0, na.rm=TRUE)
  prob_exact <- all(prob_diff == 0, na.rm=TRUE)
  pval_exact <- all(pval_diff == 0, na.rm=TRUE)
}

summary <- data.frame(
  method=method, status='ok', elapsed_sec=elapsed,
  cells=ncol(obj@data.signaling), signaling_genes=nrow(obj@data.signaling),
  groups=nlevels(obj@idents), n_lr=nrow(obj@LR$LRsig), nboot=100L,
  seed_use=1L, threads_requested=1L,
  prob_max_abs_diff=prob_max_abs_diff, pval_max_abs_diff=pval_max_abs_diff,
  prob_diff_count=prob_diff_count, pval_diff_count=pval_diff_count,
  prob_exact=prob_exact, pval_exact=pval_exact,
  probability_object_bytes=as.numeric(object.size(result@net$prob)),
  pvalue_object_bytes=as.numeric(object.size(result@net$pval)),
  CellChat_version=as.character(packageVersion('CellChat')),
  CellChatESS_version=if (requireNamespace('CellChatESS', quietly=TRUE)) as.character(packageVersion('CellChatESS')) else NA_character_,
  AutoZyme_version=if (requireNamespace('autozyme', quietly=TRUE)) as.character(packageVersion('autozyme')) else NA_character_,
  autozyme_patch_status=if (method == 'autozyme') unname(autozyme::status()[['cellchat']]) else NA_character_,
  input_path=normalizePath(input_path), reference_path=if (reference_path == '-') NA_character_ else normalizePath(reference_path),
  stringsAsFactors=FALSE
)
write.table(summary, out_path, sep='\t', row.names=FALSE, quote=TRUE)
writeLines(capture.output(sessionInfo()), sub('\\.tsv$', '_sessionInfo.txt', out_path))
cat('method=', method, ' elapsed_sec=', elapsed, ' cells=', ncol(obj@data.signaling),
    ' groups=', nlevels(obj@idents), ' n_lr=', nrow(obj@LR$LRsig), '\n', sep='')
