#!/usr/bin/env Rscript
# Memory-bounded, chunked full-array comparison and independent aggregateNet count check.
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 5L) stop('usage: ORIGINAL_RDS ESS_RDS TARGETS_TSV OUTPUT_DIR EXPECTED_LR_COUNT')
orig_path <- normalizePath(a[[1L]], mustWork = TRUE)
ess_path <- normalizePath(a[[2L]], mustWork = TRUE)
targets_path <- normalizePath(a[[3L]], mustWork = TRUE)
out <- normalizePath(a[[4L]], mustWork = TRUE)
expected_lr <- as.integer(a[[5L]])
if (!dir.exists(out)) stop('Create a fresh comparison output directory first')
if (file.exists(file.path(out, 'full627_comparison_summary.tsv'))) stop('Refusing to overwrite an existing comparison')
if (!requireNamespace('jsonlite', quietly = TRUE)) stop('jsonlite required')
orig <- readRDS(orig_path)
ess <- readRDS(ess_path)
if (!identical(names(orig), names(ess)) || !identical(names(orig), c('prob','pval','count','groups','lr_names','parameters')))
  stop('Saved output fields differ')
shared_parameters <- c('type','raw.use','population.size','nboot','seed.use','Kh','n','threads','future_workers')
if (!identical(orig$parameters[shared_parameters], ess$parameters[shared_parameters]))
  stop('Methods do not have identical shared inference parameters')
if (!identical(ess$parameters$algorithm, 'ess')) stop('ESS algorithm is not ess')
shape_ok <- identical(dim(orig$prob), dim(ess$prob)) && identical(dim(orig$pval), dim(ess$pval))
names_ok <- identical(dimnames(orig$prob), dimnames(ess$prob)) && identical(dimnames(orig$pval), dimnames(ess$pval)) &&
  identical(orig$groups, ess$groups) && identical(orig$lr_names, ess$lr_names)
if (!shape_ok || !names_ok) stop('Output dimensions or names/order differ; inspect run identities and arrays')
K <- dim(orig$prob)[1L]; L <- dim(orig$prob)[3L]
if (L != expected_lr) stop('Expected LR count ', expected_lr, '; observed ', L)
source_groups <- orig$groups; lr_names <- orig$lr_names
src_chunk <- 64L
max_abs <- -Inf; max_loc <- c(NA_integer_, NA_integer_, NA_integer_)
prob_diff_n <- 0; pval_diff_n <- 0
nf <- function(x) c(NA_count = sum(is.na(x) & !is.nan(x)), NaN_count = sum(is.nan(x)),
  positive_Inf_count = sum(is.infinite(x) & x > 0), negative_Inf_count = sum(is.infinite(x) & x < 0))
prob_nf_o <- prob_nf_e <- pval_nf_o <- pval_nf_e <- c(NA_count=0,NaN_count=0,positive_Inf_count=0,negative_Inf_count=0)
prob_nf_diff_rows <- list()
prob_diff_file <- file.path(out, 'probability_nonfinite_differences.tsv')
pval_diff_file <- file.path(out, 'pvalue_differences.tsv')
count_orig_calc <- matrix(0, nrow=K, ncol=K, dimnames=dimnames(orig$prob)[1:2])
count_ess_calc <- matrix(0, nrow=K, ncol=K, dimnames=dimnames(ess$prob)[1:2])
writeLines('lr_index\tlr_name\tsource_group\ttarget_group\toriginal_class\tESS_class', prob_diff_file)
writeLines(paste(c('lr_index','lr_name','source_group','target_group','original_value','ESS_value',
  'original_hex','ESS_hex','difference_kind','original_observed_probability','ESS_observed_probability',
  'original_observed_probability_hex','ESS_observed_probability_hex'),collapse='\t'), pval_diff_file)
value_class <- function(x) if (is.nan(x)) 'NaN' else if (is.na(x)) 'NA' else if (is.infinite(x)) if (x>0) '+Inf' else '-Inf' else 'finite'
to_text <- function(x) if (is.finite(x)) sprintf('%.17g',x) else as.character(x)
to_hex <- function(x) if (is.finite(x)) sprintf('%a',as.double(x)) else as.character(x)
nonfinite_same <- function(x,y) (is.nan(x)&is.nan(y)) | (is.na(x)&!is.nan(x)&is.na(y)&!is.nan(y)) |
  (is.infinite(x)&is.infinite(y)&sign(x)==sign(y))
write_rows <- function(path, x) {
  if (!nrow(x)) return(invisible(NULL))
  con <- file(path, open='at'); on.exit(close(con),add=TRUE)
  write.table(x, con, sep='\t', row.names=FALSE, col.names=FALSE, quote=TRUE, na='NA')
}
for (ell in seq_len(L)) {
  for (s0 in seq.int(1L,K,by=src_chunk)) {
    sr <- s0:min(K,s0+src_chunk-1L)
    po <- orig$prob[sr,,ell,drop=FALSE][,,1L]
    pe <- ess$prob[sr,,ell,drop=FALSE][,,1L]
    pvo <- orig$pval[sr,,ell,drop=FALSE][,,1L]
    pve <- ess$pval[sr,,ell,drop=FALSE][,,1L]
    prob_nf_o <- prob_nf_o + nf(po); prob_nf_e <- prob_nf_e + nf(pe)
    pval_nf_o <- pval_nf_o + nf(pvo); pval_nf_e <- pval_nf_e + nf(pve)
    prob_finite <- is.finite(po) & is.finite(pe)
    prob_unequal <- (prob_finite & (po != pe)) | (!prob_finite & !nonfinite_same(po,pe))
    prob_diff_n <- prob_diff_n + sum(prob_unequal)
    d <- abs(po-pe); d[!prob_finite] <- NA_real_
    if (any(is.finite(d))) {
      md <- max(d[is.finite(d)])
      if (md > max_abs) {
        q <- which(d == md, arr.ind=TRUE)[1L,]
        max_abs <- md; max_loc <- c(sr[q[1L]],q[2L],ell)
      }
    }
    nf_mismatch <- !prob_finite & !nonfinite_same(po,pe)
    if (any(nf_mismatch)) {
      ix <- which(nf_mismatch,arr.ind=TRUE)
      rows <- data.frame(lr_index=ell,lr_name=lr_names[ell],source_group=source_groups[sr[ix[,1L]]],
        target_group=source_groups[ix[,2L]],original_class=vapply(po[nf_mismatch],value_class,''),
        ESS_class=vapply(pe[nf_mismatch],value_class,''),stringsAsFactors=FALSE)
      write_rows(prob_diff_file,rows)
    }
    p_finite <- is.finite(pvo) & is.finite(pve)
    p_unequal <- (p_finite & (pvo != pve)) | (!p_finite & !nonfinite_same(pvo,pve))
    pval_diff_n <- pval_diff_n + sum(p_unequal)
    if (any(p_unequal)) {
      ix <- which(p_unequal,arr.ind=TRUE)
      ov <- pvo[p_unequal]; ev <- pve[p_unequal]
      kinds <- ifelse(p_finite[p_unequal],'finite_unequal',paste0('nonfinite:',vapply(seq_along(ov),function(i) paste(value_class(ov[i]),value_class(ev[i]),sep='_vs_'),'')))
      rows <- data.frame(lr_index=ell,lr_name=lr_names[ell],source_group=source_groups[sr[ix[,1L]]],
        target_group=source_groups[ix[,2L]],original_value=vapply(ov,to_text,''),ESS_value=vapply(ev,to_text,''),
        original_hex=vapply(ov,to_hex,''),ESS_hex=vapply(ev,to_hex,''),difference_kind=kinds,
        original_observed_probability=vapply(po[p_unequal],to_text,''),
        ESS_observed_probability=vapply(pe[p_unequal],to_text,''),
        original_observed_probability_hex=vapply(po[p_unequal],to_hex,''),
        ESS_observed_probability_hex=vapply(pe[p_unequal],to_hex,''),
        stringsAsFactors=FALSE)
      write_rows(pval_diff_file,rows)
    }
    # Reproduce CellChat 2.2.0.9001 aggregateNet's default rule per chunk:
    # pval[prob == 0] <- 1; prob[pval >= 0.05] <- 0; apply(prob > 0, c(1,2), sum).
    pvo[po == 0] <- 1; pve[pe == 0] <- 1
    po[pvo >= 0.05] <- 0; pe[pve >= 0.05] <- 0
    # The full CellChat apply() reduces the third (LR) dimension. At a single-LR
    # page this is exactly an elementwise 0/1 increment, including NA propagation.
    count_orig_calc[sr,] <- count_orig_calc[sr,] + (po > 0)
    count_ess_calc[sr,] <- count_ess_calc[sr,] + (pe > 0)
  }
}
count_orig_calc[is.na(count_orig_calc)] <- 0
count_ess_calc[is.na(count_ess_calc)] <- 0
storage.mode(count_orig_calc) <- 'integer'
storage.mode(count_ess_calc) <- 'integer'
count_orig_saved <- orig$count; count_ess_saved <- ess$count
count_matches_saved <- function(recalc,saved) {
  if (!identical(dim(recalc),dim(saved)) || !identical(dimnames(recalc),dimnames(saved))) return(FALSE)
  eq <- as.vector(recalc) == as.vector(saved)
  !anyNA(eq) && all(eq)
}
orig_count_recalc_ok <- count_matches_saved(count_orig_calc,count_orig_saved)
ess_count_recalc_ok <- count_matches_saved(count_ess_calc,count_ess_saved)
recalc_diff <- list()
for (spec in list(list('Original',count_orig_calc,count_orig_saved),list('ESS',count_ess_calc,count_ess_saved))) {
  method_name <- spec[[1L]]; calc <- spec[[2L]]; saved <- spec[[3L]]
  if (identical(dim(calc),dim(saved))) {
    bad <- is.na(calc) != is.na(saved) | (!is.na(calc) & !is.na(saved) & calc != saved)
    bad[is.na(bad)] <- TRUE
    if (any(bad)) {
      ix <- which(bad,arr.ind=TRUE)
      recalc_diff[[length(recalc_diff)+1L]] <- data.frame(method=method_name,
        source_group=rownames(calc)[ix[,1L]],target_group=colnames(calc)[ix[,2L]],
        recalculated_count=calc[ix],saved_aggregateNet_count=saved[ix],stringsAsFactors=FALSE)
    }
  }
}
if (length(recalc_diff)) write.table(do.call(rbind,recalc_diff),file.path(out,'aggregate_count_recomputation_differences.tsv'),
  sep='\t',row.names=FALSE,quote=TRUE) else writeLines('method\tsource_group\ttarget_group\trecalculated_count\tsaved_aggregateNet_count',
  file.path(out,'aggregate_count_recomputation_differences.tsv'))
count_same <- identical(dim(count_orig_calc),dim(count_ess_calc)) &&
  identical(dimnames(count_orig_calc),dimnames(count_ess_calc)) && identical(count_orig_calc,count_ess_calc)
count_diff <- if (identical(dim(count_orig_calc),dim(count_ess_calc))) count_orig_calc != count_ess_calc else matrix(TRUE,K,K)
count_diff[is.na(count_diff)] <- TRUE
if (any(count_diff)) {
  ix <- which(count_diff,arr.ind=TRUE)
  write.table(data.frame(source_group=source_groups[ix[,1]],target_group=source_groups[ix[,2]],
    original_count=count_orig_calc[ix],ESS_count=count_ess_calc[ix]),
    file.path(out,'aggregate_count_differences.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
} else writeLines('source_group\ttarget_group\toriginal_count\tESS_count',file.path(out,'aggregate_count_differences.tsv'))

# Independent paper-workflow check, after inference and all performance timers.
# CellChat modeling.R:916-920 excludes groups with <= min.cells, not < min.cells.
# For default min.samples=NULL (effective 1), filterCommunication changes only
# these full sender/receiver pages. Masking the final count axes is equivalent.
sizes_orig_path <- file.path(dirname(orig_path),'group_sizes.tsv')
sizes_ess_path <- file.path(dirname(ess_path),'group_sizes.tsv')
if (file.exists(sizes_orig_path) && file.exists(sizes_ess_path)) {
  sizes_o <- read.delim(sizes_orig_path,stringsAsFactors=FALSE)
  sizes_e <- read.delim(sizes_ess_path,stringsAsFactors=FALSE)
  stopifnot(identical(sizes_o,sizes_e),identical(sizes_o$group,source_groups))
  excluded <- which(sizes_o$n_cells <= 10L)
  filtered_o <- count_orig_calc; filtered_e <- count_ess_calc
  if (length(excluded)) {
    filtered_o[excluded,] <- 0L; filtered_o[,excluded] <- 0L
    filtered_e[excluded,] <- 0L; filtered_e[,excluded] <- 0L
  }
  filtered_bad <- filtered_o != filtered_e
  filtered_bad[is.na(filtered_bad)] <- TRUE
  write.table(data.frame(policy=c('aggregateNet_default','filterCommunication_min_cells_10_then_aggregateNet'),
    pvalue_threshold=0.05,min_cells=c(NA_integer_,10L),
    excluded_groups=c(0L,length(excluded)),
    original_total_significant_LR_edges=c(sum(count_orig_calc),sum(filtered_o)),
    ESS_total_significant_LR_edges=c(sum(count_ess_calc),sum(filtered_e)),
    different_group_pairs=c(sum(count_diff),sum(filtered_bad)),
    counts_identical=c(count_same,identical(filtered_o,filtered_e))),
    file.path(out,'aggregate_count_policy_summary.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
  write.table(sizes_o[excluded,,drop=FALSE],file.path(out,'groups_excluded_by_min_cells_10.tsv'),
    sep='\t',row.names=FALSE,quote=TRUE)
  if (any(filtered_bad)) {
    ix <- which(filtered_bad,arr.ind=TRUE)
    write.table(data.frame(source_group=source_groups[ix[,1]],target_group=source_groups[ix[,2]],
      original_count=filtered_o[ix],ESS_count=filtered_e[ix]),
      file.path(out,'aggregate_count_min_cells_10_differences.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
  } else writeLines('source_group\ttarget_group\toriginal_count\tESS_count',
    file.path(out,'aggregate_count_min_cells_10_differences.tsv'))
  saveRDS(list(count=filtered_o,thresh=0.05,min.cells=10L,excluded_groups=sizes_o$group[excluded]),
    file.path(out,'original_aggregate_count_min_cells_10.rds'),compress='gzip')
  saveRDS(list(count=filtered_e,thresh=0.05,min.cells=10L,excluded_groups=sizes_o$group[excluded]),
    file.path(out,'ESS_aggregate_count_min_cells_10.rds'),compress='gzip')
}

if (!is.finite(max_abs)) max_abs <- NA_real_
max_info <- if (all(is.finite(max_loc))) list(lr_index=max_loc[3],lr_name=lr_names[max_loc[3]],
  source_group=source_groups[max_loc[1]],target_group=source_groups[max_loc[2]],
  original_value=sprintf('%.17g',orig$prob[max_loc[1],max_loc[2],max_loc[3]]),
  ESS_value=sprintf('%.17g',ess$prob[max_loc[1],max_loc[2],max_loc[3]]),
  original_hex=sprintf('%a',orig$prob[max_loc[1],max_loc[2],max_loc[3]]),
  ESS_hex=sprintf('%a',ess$prob[max_loc[1],max_loc[2],max_loc[3]])) else list(lr_index=NA,lr_name=NA,source_group=NA,target_group=NA,original_value=NA,ESS_value=NA,original_hex=NA,ESS_hex=NA)
summary <- data.frame(check=c('dimensions_match','names_and_order_match','probability_entries','probability_different_elements',
  'maximum_absolute_probability_difference','pvalue_entries','pvalue_different_elements','aggregate_count_entries',
  'original_aggregate_count_recomputed_matches_saved','ESS_aggregate_count_recomputed_matches_saved',
  'aggregate_counts_identical','original_probability_NA','original_probability_NaN','original_probability_posInf','original_probability_negInf',
  'ESS_probability_NA','ESS_probability_NaN','ESS_probability_posInf','ESS_probability_negInf',
  'original_pvalue_NA','original_pvalue_NaN','original_pvalue_posInf','original_pvalue_negInf',
  'ESS_pvalue_NA','ESS_pvalue_NaN','ESS_pvalue_posInf','ESS_pvalue_negInf'),
  value=as.character(c(TRUE,TRUE,length(orig$prob),prob_diff_n,
    if (is.finite(max_abs)) sprintf('%.17g',max_abs) else 'NA',length(orig$pval),pval_diff_n,length(count_orig_calc),
    orig_count_recalc_ok,ess_count_recalc_ok,count_same,
    prob_nf_o,prob_nf_e,pval_nf_o,pval_nf_e)),stringsAsFactors=FALSE)
# Keep typed values as explicit key/value rows to avoid coercion ambiguity.
write.table(summary,file.path(out,'full627_comparison_summary.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
write.table(data.frame(metric=c('max_probability_lr_index','max_probability_lr_name','max_probability_source_group',
  'max_probability_target_group','max_probability_original','max_probability_ESS','max_probability_original_hex','max_probability_ESS_hex'),
  value=unlist(max_info,use.names=FALSE),stringsAsFactors=FALSE),
  file.path(out,'maximum_probability_difference.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
write.table(data.frame(method=rep(c('Original','ESS'),each=8),
  array=rep(c(rep('probability',4),rep('pvalue',4)),times=2),
  class=rep(c('NA','NaN','positive_Inf','negative_Inf'),times=4),
  n=c(prob_nf_o,pval_nf_o,prob_nf_e,pval_nf_e)),file.path(out,'nonfinite_summary.tsv'),
  sep='\t',row.names=FALSE,quote=TRUE)

if (expected_lr == 627L && file.exists(targets_path) && file.info(targets_path)$size > 0L) {
  targets <- read.delim(targets_path,check.names=FALSE,stringsAsFactors=FALSE)
  if (all(c('lr_name','source_group','target_group') %in% names(targets))) {
    rr <- lapply(seq_len(nrow(targets)),function(i) {
      t <- targets[i,]; ell <- match(t$lr_name,lr_names); si <- match(t$source_group,source_groups); ri <- match(t$target_group,source_groups)
      if (anyNA(c(ell,si,ri))) stop('Could not locate target by LR+sender+receiver: ',t$lr_name,'/',t$source_group,'/',t$target_group)
      po <- orig$prob[si,ri,ell]; pe <- ess$prob[si,ri,ell]; pvo <- orig$pval[si,ri,ell]; pve <- ess$pval[si,ri,ell]
      data.frame(lr_name=t$lr_name,source_group=t$source_group,target_group=t$target_group,lr_index_full_627=ell,
        original_probability=sprintf('%.17g',po),ESS_probability=sprintf('%.17g',pe),
        probability_absolute_difference=sprintf('%.17g',abs(po-pe)),original_probability_hex=sprintf('%a',po),ESS_probability_hex=sprintf('%a',pe),
        original_pvalue=sprintf('%.17g',pvo),ESS_pvalue=sprintf('%.17g',pve),
        original_pvalue_hex=sprintf('%a',pvo),ESS_pvalue_hex=sprintf('%a',pve),
        pvalue_identical=identical(as.double(pvo),as.double(pve)),
        original_nReject_from_p_if_probability_positive=if (is.finite(po) && po>0) round(100*pvo) else NA_real_,
        ESS_nReject_from_p_if_probability_positive=if (is.finite(pe) && pe>0) round(100*pve) else NA_real_,
        original_aggregate_count_at_pair=count_orig_calc[si,ri],ESS_aggregate_count_at_pair=count_ess_calc[si,ri],
        stringsAsFactors=FALSE)
    })
    write.table(do.call(rbind,rr),file.path(out,'five_prior_targets.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
  }
}
cat('Comparison complete: probability unequal=',prob_diff_n,'; max absolute difference=',sprintf('%.17g',max_abs),
  '; p-value unequal=',pval_diff_n,'; aggregate counts identical=',count_same,'\n',sep='')
if (!orig_count_recalc_ok || !ess_count_recalc_ok)
  stop('Independent aggregate count recomputation disagrees with saved aggregateNet count; do not accept this comparison')
