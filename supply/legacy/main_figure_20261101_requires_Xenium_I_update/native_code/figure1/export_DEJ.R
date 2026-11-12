a <- commandArgs(TRUE)
stopifnot(length(a) == 2L)
src <- a[1]; out <- a[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
original <- readRDS(file.path(src, 'original/result.rds'))
ess <- readRDS(file.path(src, 'ess/result.rds'))
stopifnot(identical(dim(original$prob), c(9L, 9L, 1128L)),
          identical(dimnames(original$prob), dimnames(ess$prob)),
          identical(original$pval, ess$pval), identical(original$count, ess$count),
          all(is.finite(original$prob)), all(is.finite(ess$prob)))
write_tsv <- function(x, name) write.table(x, file.path(out, name), sep = '\t',
                                         row.names = FALSE, quote = TRUE)
all_long <- list(); weights <- list(); counts <- list()
for (engine in c('Original CellChat', 'CellChatESS')) {
  z <- if (engine == 'Original CellChat') original else ess
  stopifnot(identical(dimnames(z$prob), dimnames(z$pval)))
  dn <- dimnames(z$prob)
  ind <- which(z$prob > 0, arr.ind = TRUE)
  all_long[[engine]] <- data.frame(engine = engine,
    source_id = dn[[1]][ind[,1]], target_id = dn[[2]][ind[,2]],
    interaction_name = dn[[3]][ind[,3]],
    probability = sprintf('%.17g', z$prob[ind]),
    pvalue = sprintf('%.17g', z$pval[ind]), stringsAsFactors = FALSE)
  filtered <- z$prob
  filtered[z$pval >= 0.05] <- 0
  weight <- apply(filtered, c(1, 2), sum)
  count <- apply(filtered > 0, c(1, 2), sum)
  stopifnot(identical(count, z$count))
  ij <- which(matrix(TRUE, nrow = 9, ncol = 9), arr.ind = TRUE)
  weights[[engine]] <- data.frame(engine = engine, source_id = dn[[1]][ij[,1]],
    target_id = dn[[2]][ij[,2]], weight = sprintf('%.17g', weight[ij]))
  counts[[engine]] <- data.frame(engine = engine, source_id = dn[[1]][ij[,1]],
    target_id = dn[[2]][ij[,2]], count = count[ij])
}
write_tsv(do.call(rbind, all_long), 'DEJ_nonzero_probabilities.tsv')
write_tsv(do.call(rbind, weights), 'DE_network_weights.tsv')
write_tsv(do.call(rbind, counts), 'DE_network_counts.tsv')
metrics <- data.frame(metric = c('probability_elements', 'probability_max_abs_difference',
  'probability_different_elements', 'pvalue_different_elements', 'count_different_elements'),
  value = c(as.character(length(original$prob)),
    sprintf('%.17g', max(abs(original$prob - ess$prob))),
    as.character(sum(original$prob != ess$prob)),
    as.character(sum(original$pval != ess$pval)), as.character(sum(original$count != ess$count))))
write_tsv(metrics, 'DEJ_numerical_comparison.tsv')
cat('Export complete from saved results; no inference rerun.\n')
