args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop('usage: trace_xenium_ess_bootstraps.R INPUT_RDS SAVED_FIVE_LR_RDS DIAGNOSTICS_TSV OUTDIR')
input_path <- args[[1L]]
saved_path <- args[[2L]]
diagnostics_path <- args[[3L]]
out_dir <- args[[4L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

profile_lib <- Sys.getenv("CELLCHAT_TRACE_R_LIB", "environment/historical_trace/library")
cellchat16_lib <- Sys.getenv("CELLCHAT_TRACE_R_LIB", "environment/historical_trace/library")
.libPaths(c(profile_lib, cellchat16_lib, .libPaths()))
suppressPackageStartupMessages({ library(CellChat); library(CellChatESS) })
if (as.character(packageVersion('CellChat')) != '1.6.1') stop('Expected CellChat 1.6.1')
if (as.character(packageVersion('CellChatESS')) != '0.2.2') stop('Expected profiled CellChatESS 0.2.2')

obj <- readRDS(input_path)
saved <- readRDS(saved_path)
diagnostics <- read.delim(diagnostics_path, check.names = FALSE, stringsAsFactors = FALSE)
all_pairs <- obj@LR$LRsig
lr_order <- saved$lr_names
pair_rows <- match(lr_order, rownames(all_pairs))
if (anyNA(pair_rows)) stop('A saved LR name is missing from the prepared Xenium object')
pairs <- all_pairs[pair_rows, , drop = FALSE]
obj@LR$LRsig <- pairs

complex <- obj@DB$complex
cofactor <- obj@DB$cofactor
cofactor_cols <- grep('cofactor', colnames(cofactor), value = TRUE)
resolve_items <- function(items) {
  needed <- character()
  for (x in items) {
    if (is.na(x) || !nzchar(x)) next
    if (x %in% rownames(complex)) {
      genes <- unlist(complex[x, grep('^subunit', colnames(complex)), drop = FALSE], use.names = FALSE)
      needed <- c(needed, genes[!is.na(genes) & nzchar(genes)])
    } else if (x %in% rownames(cofactor)) {
      genes <- unlist(cofactor[x, cofactor_cols, drop = FALSE], use.names = FALSE)
      needed <- c(needed, genes[!is.na(genes) & nzchar(genes)])
    } else needed <- c(needed, x)
  }
  unique(needed)
}
items <- unique(unlist(pairs[, intersect(c('ligand', 'receptor', 'co_A_receptor', 'co_I_receptor',
  'agonist', 'antagonist'), colnames(pairs)), drop = FALSE], use.names = FALSE))
needed <- intersect(resolve_items(items), rownames(obj@data.signaling))
global_max <- max(obj@data.signaling)
max_position <- which(obj@data.signaling == global_max, arr.ind = TRUE)[1L, 1L]
max_gene <- rownames(obj@data.signaling)[max_position]
keep_genes <- unique(c(needed, max_gene))
obj@data.signaling <- obj@data.signaling[keep_genes, , drop = FALSE]
if (max(obj@data.signaling) != global_max) stop('Gene projection changed the global normalization')

K <- nlevels(obj@idents)
targets <- vector('list', nrow(diagnostics))
target_map <- file.path(out_dir, 'ess_trace_targets.tsv')
for (j in seq_len(nrow(diagnostics))) {
  lr_idx <- match(diagnostics$lr_name[[j]], lr_order)
  sender <- match(diagnostics$source_group[[j]], levels(obj@idents))
  receiver <- match(diagnostics$target_group[[j]], levels(obj@idents))
  if (anyNA(c(lr_idx, sender, receiver))) stop('Failed to resolve a target LR or group')
  pair_idx0 <- (sender - 1L) + K * (receiver - 1L)
  targets[[j]] <- data.frame(discrepancy_id = j, lr_name = diagnostics$lr_name[[j]],
    lr_index_1based = lr_idx, pair_index_0based = pair_idx0,
    source_group = diagnostics$source_group[[j]], target_group = diagnostics$target_group[[j]])
}
target_table <- unique(do.call(rbind, targets))
write.table(target_table, file.path(out_dir, 'ess_trace_targets_resolved.tsv'), sep = '\t',
  row.names = FALSE, quote = TRUE)
write.table(data.frame(lr_index_0based = target_table$lr_index_1based - 1L,
  pair_index_0based = target_table$pair_index_0based), target_map, sep = '\t',
  row.names = FALSE, col.names = FALSE, quote = FALSE)
trace_path <- file.path(out_dir, 'xenium_ess_bootstrap_trace_raw.tsv')
write.table(data.frame(lr_index_1based = 'lr_index_1based', pair_index_0based = 'pair_index_0based',
  bootstrap = 'bootstrap', pnull_ess = 'pnull_ess', pboot_ess = 'pboot_ess',
  delta_ess = 'delta_ess', reject_ess = 'reject_ess'), trace_path, sep = '\t',
  row.names = FALSE, col.names = FALSE, quote = FALSE)
Sys.setenv(CELLCHATESS_TRACE_TARGETS = target_map, CELLCHATESS_TRACE_OUTPUT = trace_path,
  OMP_NUM_THREADS = '1', OPENBLAS_NUM_THREADS = '1', MKL_NUM_THREADS = '1')

cat('Tracing ', nrow(pairs), ' LR rows; ', nrow(obj@data.signaling), ' genes; ',
  ncol(obj@data.signaling), ' cells; ', K, ' groups; ', nrow(target_table),
  ' unique LR/group targets.\n', sep = '')
start <- proc.time()[['elapsed']]
ess <- CellChatESS::computeCommunProbESS(obj, type = 'triMean', raw.use = TRUE,
  population.size = FALSE, nboot = 100L, seed.use = 1L, Kh = 0.5, n = 1,
  algorithm = 'ess', threads = 1L)
elapsed <- proc.time()[['elapsed']] - start

output_checks <- lapply(seq_len(nrow(target_table)), function(i) {
  z <- target_table[i, ]
  j <- match(z$lr_name, lr_order)
  s <- match(z$source_group, dimnames(ess@net$prob)[[1L]])
  r <- match(z$target_group, dimnames(ess@net$prob)[[2L]])
  ess_p <- ess@net$pval[s, r, j]
  saved_p <- saved$ess_pval[s, r, j]
  data.frame(discrepancy_id = z$discrepancy_id, lr_name = z$lr_name,
    source_group = z$source_group, target_group = z$target_group,
    ess_profiled_probability = ess@net$prob[s, r, j], saved_ess_probability = saved$ess_prob[s, r, j],
    ess_profiled_pvalue = ess_p, saved_ess_pvalue = saved_p,
    ess_profiled_reject_count = as.integer(round(ess_p * 100)),
    saved_ess_reject_count = as.integer(round(saved_p * 100)),
    pvalue_matches_saved = identical(as.numeric(ess_p), as.numeric(saved_p)),
    stringsAsFactors = FALSE)
})
write.table(do.call(rbind, output_checks), file.path(out_dir, 'ess_trace_output_validation.tsv'),
  sep = '\t', row.names = FALSE, quote = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, 'sessionInfo.txt'))
write.table(data.frame(elapsed_sec = elapsed, cells = ncol(obj@data.signaling),
  genes = nrow(obj@data.signaling), groups = K, n_lr = nrow(pairs), nboot = 100L,
  CellChat_version = as.character(packageVersion('CellChat')),
  CellChatESS_version = as.character(packageVersion('CellChatESS'))),
  file.path(out_dir, 'trace_run_summary.tsv'), sep = '\t', row.names = FALSE, quote = TRUE)
cat('ESS bootstrap trace complete in ', elapsed, ' sec.\n', sep = '')
