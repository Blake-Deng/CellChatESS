args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=5L) stop('usage: run_ablation_method.R METHOD LIBRARY INPUT_RDS ORIGINAL_RDS OUTPUT_TSV')
method <- args[[1L]]; pkg_lib <- args[[2L]]; input_path <- args[[3L]]; ref_path <- args[[4L]]; output_path <- args[[5L]]
archived_lib <- Sys.getenv("CELLCHAT_MAIN_R_LIB", "environment/main/library")
.libPaths(c(pkg_lib, archived_lib, .libPaths()))
suppressPackageStartupMessages({library(CellChat); library(CellChatESS)})
obj <- readRDS(input_path); ref <- readRDS(ref_path)
if (!identical(as.character(packageVersion('CellChat')), '1.6.1')) stop('Unexpected CellChat version')
if (!as.character(packageVersion('CellChatESS')) %in% c('0.1.1.9000','0.1.1')) stop('Unexpected CellChatESS version')
params <- ref@options$parameter
param_or <- function(key, default) { value <- params[[key]]; if (is.null(value)) default else value }
nboot <- as.integer(param_or('nboot',100L)); seed.use <- as.integer(param_or('seed.use',1L))
Kh <- as.numeric(param_or('Kh',0.5)); n_power <- as.numeric(param_or('n',1))
start <- proc.time()[['elapsed']]
if (method == 'dense_projected') {
  result <- CellChatESS::computeCommunProbESS(obj, type='triMean', raw.use=TRUE, population.size=FALSE,
    nboot=nboot, seed.use=seed.use, Kh=Kh, n=n_power, algorithm='dense')
} else if (method %in% c('ess_pruned','ess_no_pruning','ess_no_cache')) {
  result <- CellChatESS::computeCommunProbESS(obj, type='triMean', raw.use=TRUE, population.size=FALSE,
    nboot=nboot, seed.use=seed.use, Kh=Kh, n=n_power, algorithm='ess')
} else if (method == 'dense_full_gene') {
  data.raw <- obj@data.signaling
  max.data <- max(data.raw)
  if (!is.finite(max.data) || max.data <= 0) stop('invalid maximum expression')
  data.full <- as.matrix(data.raw / max.data)
  lr <- obj@LR$LRsig; data.genes <- rownames(data.full)
  cofactor <- obj@DB$cofactor; complex <- obj@DB$complex
  cofactor.cols <- grepl('cofactor', colnames(cofactor))
  make.index <- getFromNamespace('make_index_matrix','CellChatESS')
  ligand.idx <- make.index(as.character(lr$ligand), data.genes, complex_input=complex, keep_duplicates=TRUE)
  receptor.idx <- make.index(as.character(lr$receptor), data.genes, complex_input=complex, keep_duplicates=TRUE)
  coA.idx <- make.index(as.character(lr$co_A_receptor), data.genes, cofactor_input=cofactor, cofactor_cols=cofactor.cols)
  coI.idx <- make.index(as.character(lr$co_I_receptor), data.genes, cofactor_input=cofactor, cofactor_cols=cofactor.cols)
  agonist.names <- as.character(lr$agonist); antagonist.names <- as.character(lr$antagonist)
  has.agonist <- !is.na(agonist.names) & nzchar(agonist.names)
  has.antagonist <- !is.na(antagonist.names) & nzchar(antagonist.names)
  agonist.idx <- make.index(agonist.names, data.genes, cofactor_input=cofactor, cofactor_cols=cofactor.cols)
  antagonist.idx <- make.index(antagonist.names, data.genes, cofactor_input=cofactor, cofactor_cols=cofactor.cols)
  agonist.idx[!has.agonist,] <- 0L; antagonist.idx[!has.antagonist,] <- 0L
  group <- obj@idents; group.int <- as.integer(group); K <- nlevels(group); nC <- ncol(data.full)
  set.seed(seed.use)
  group.boot <- replicate(nboot, group.int[sample.int(nC, size=nC)])
  group.boot <- matrix(group.boot, nrow=nC, ncol=nboot)
  avg <- CellChatESS:::group_tri_mean_cpp(data.full, group.int, K)
  avg.boot <- CellChatESS:::group_tri_mean_boot_cpp(data.full, group.boot, K)$avg_boot
  res <- CellChatESS:::cellchat_prob_from_avg_cpp(avg, avg.boot, ligand.idx, receptor.idx, coA.idx, coI.idx,
    agonist.idx, antagonist.idx, has.agonist, has.antagonist, Kh, n_power)
  result <- obj; dimnames(res$prob) <- list(levels(group),levels(group),rownames(lr)); dimnames(res$pval) <- dimnames(res$prob)
  result@net <- list(prob=res$prob,pval=res$pval)
} else stop('unknown method: ', method)
elapsed <- proc.time()[['elapsed']] - start
prob.ref <- ref@net$prob; pval.ref <- ref@net$pval
prob <- result@net$prob; pval <- result@net$pval
if (!identical(dim(prob), dim(prob.ref)) || !identical(dim(pval), dim(pval.ref))) stop('output dimensions differ from reference')
pd <- max(abs(prob-prob.ref),na.rm=TRUE); qd <- max(abs(pval-pval.ref),na.rm=TRUE)
row <- data.frame(method=method, experiment_id=basename(dirname(input_path)),
  n_cells=ncol(obj@data.signaling), n_groups=nlevels(obj@idents), n_lr=nrow(obj@LR$LRsig), nboot=nboot,
  seed_use=seed.use, elapsed_sec=elapsed, prob_max_abs_diff_vs_original=pd,
  pval_max_abs_diff_vs_original=qd, prob_exact=sum(prob!=prob.ref,na.rm=TRUE)==0L,
  pval_exact=sum(pval!=pval.ref,na.rm=TRUE)==0L,
  active_fraction=if (is.null(result@options$CellChatESS)) NA_real_ else result@options$CellChatESS$active_fraction,
  stringsAsFactors=FALSE)
write.table(row,output_path,sep='\t',row.names=FALSE,quote=TRUE)
cat('method=',method,' elapsed_sec=',elapsed,' prob_max_diff=',pd,' pval_max_diff=',qd,'\n',sep='')
