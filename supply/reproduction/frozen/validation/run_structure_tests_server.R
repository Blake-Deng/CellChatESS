args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop('usage: run_support_pruning_matrix.R OUTDIR')
out_dir <- args[[1L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(out_dir, 'support_pruning_matrix_results.tsv')))
  stop('Refusing to overwrite prior support-pruning matrix results')

ess_lib <- Sys.getenv("CELLCHAT_ESS_R_LIB", "environment/compatible/library")
cellchat_lib <- Sys.getenv("CELLCHAT_COMPAT_R_LIB", "environment/compatible/library")
.libPaths(c(ess_lib, cellchat_lib, .libPaths()))
suppressPackageStartupMessages({ library(CellChat); library(CellChatESS); library(Matrix) })
if (as.character(packageVersion('CellChat')) != Sys.getenv('CELLCHAT_REFERENCE_VERSION','2.2.0.9001')) stop('Unexpected CellChat reference version')
if (as.character(packageVersion('CellChatESS')) != '0.2.2') stop('Expected CellChatESS 0.2.2')

input <- Sys.getenv("CELLCHAT_TEST_INPUT", "data/inputs/3CA_data__10X__cells-1000__rep-1__engine-both__ablation-full/prepared_cellchat.rds")
base <- readRDS(input)
pairs <- base@LR$LRsig
complex <- base@DB$complex
cofactor <- base@DB$cofactor
cofactor_cols <- grep('cofactor', colnames(cofactor), value = TRUE)
has_value <- function(x) !is.na(x) && nzchar(as.character(x))
component_genes <- function(x) {
  if (!has_value(x)) return(character())
  x <- as.character(x)
  if (x %in% rownames(complex)) {
    z <- unlist(complex[x, grep('^subunit', colnames(complex)), drop = FALSE], use.names = FALSE)
  } else if (x %in% rownames(cofactor)) {
    z <- unlist(cofactor[x, cofactor_cols, drop = FALSE], use.names = FALSE)
  } else z <- x
  unique(intersect(z[!is.na(z) & nzchar(z)], rownames(base@data.signaling)))
}
positive_components <- function(x) {
  genes <- component_genes(x)
  if (!length(genes)) return(character())
  sums <- Matrix::rowSums(base@data.signaling[genes, , drop = FALSE])
  genes[as.numeric(sums) > 0]
}
reduce_to_lr_features <- function(obj, selected) {
  data <- obj@data.signaling
  normalizer <- max(data)
  if (!is.finite(normalizer) || normalizer <= 0) stop('Invalid expression normalizer')
  # Keep one row containing the full-input maximum so CellChat's own max(data)
  # normalization remains unchanged after feature projection.
  anchor_pos <- which(as.matrix(data) == normalizer, arr.ind = TRUE)[1L, 1L]
  anchor_gene <- rownames(data)[[anchor_pos]]
  fields <- intersect(c('ligand', 'receptor', 'co_A_receptor', 'co_I_receptor',
                        'agonist', 'antagonist'), colnames(selected))
  needed <- unique(unlist(lapply(fields, function(field)
    component_genes(as.character(selected[[field]][[1L]]))), use.names = FALSE))
  keep <- unique(c(intersect(needed, rownames(data)), anchor_gene))
  obj@data.signaling <- (data / normalizer)[keep, , drop = FALSE]
  list(object = obj, n_genes = length(keep), normalizer = normalizer,
       anchor_gene = anchor_gene)
}
probe <- CellChat::computeCommunProb(base, type = 'triMean', raw.use = TRUE,
  population.size = FALSE, nboot = 1L, seed.use = 1L, Kh = 0.5, n = 1)
active_ids <- which(vapply(seq_len(dim(probe@net$prob)[[3L]]), function(i)
  any(probe@net$prob[, , i] > 0), logical(1)))
ligand_complex <- vapply(pairs$ligand, function(x) length(component_genes(x)) >= 2L, logical(1))
receptor_complex <- vapply(pairs$receptor, function(x) length(component_genes(x)) >= 2L, logical(1))
modifier_has <- function(field) vapply(pairs[[field]], function(x) length(positive_components(x)) > 0L, logical(1))
simple <- !vapply(pairs$ligand, function(x) as.character(x) %in% rownames(complex), logical(1)) &
          !vapply(pairs$receptor, function(x) as.character(x) %in% rownames(complex), logical(1))

criteria <- list(
  single_ligand_receptor = list(ids = active_ids[simple[active_ids]], zero_field = 'ligand', zero_mode = 'all_components'),
  multimeric_ligand = list(ids = if (length(active_ids[ligand_complex[active_ids]]))
    active_ids[ligand_complex[active_ids]] else which(ligand_complex),
    zero_field = 'ligand', zero_mode = 'one_component'),
  multimeric_receptor = list(ids = active_ids[receptor_complex[active_ids]], zero_field = 'receptor', zero_mode = 'one_component'),
  co_activator = list(ids = active_ids[vapply(seq_along(active_ids), function(j)
    has_value(pairs$co_A_receptor[[active_ids[[j]]]]) && modifier_has('co_A_receptor')[[active_ids[[j]]]], logical(1))],
    zero_field = 'co_A_receptor', zero_mode = 'one_component'),
  co_inhibitor = list(ids = active_ids[vapply(seq_along(active_ids), function(j)
    has_value(pairs$co_I_receptor[[active_ids[[j]]]]) && modifier_has('co_I_receptor')[[active_ids[[j]]]], logical(1))],
    zero_field = 'co_I_receptor', zero_mode = 'one_component'),
  agonist = list(ids = active_ids[vapply(seq_along(active_ids), function(j)
    has_value(pairs$agonist[[active_ids[[j]]]]) && modifier_has('agonist')[[active_ids[[j]]]], logical(1))],
    zero_field = 'agonist', zero_mode = 'one_component'),
  antagonist = list(ids = active_ids[vapply(seq_along(active_ids), function(j)
    has_value(pairs$antagonist[[active_ids[[j]]]]) && modifier_has('antagonist')[[active_ids[[j]]]], logical(1))],
    zero_field = 'antagonist', zero_mode = 'one_component')
)
for (nm in names(criteria)) if (!length(criteria[[nm]]$ids)) stop('No candidate found for ', nm)
selected_ids <- vapply(criteria, function(x) x$ids[[1L]], integer(1))
synthetic_ids <- setNames(!selected_ids %in% active_ids, names(selected_ids))
write.table(data.frame(structure = names(selected_ids), lr_index = selected_ids,
  lr_name = rownames(pairs)[selected_ids], ligand = pairs$ligand[selected_ids],
  receptor = pairs$receptor[selected_ids],
  co_A_receptor = pairs$co_A_receptor[selected_ids],
  co_I_receptor = pairs$co_I_receptor[selected_ids],
  agonist = pairs$agonist[selected_ids], antagonist = pairs$antagonist[selected_ids],
  synthetic_expression_priming = synthetic_ids),
  file.path(out_dir, 'support_pruning_selected_lr.tsv'),
  sep = '\t', row.names = FALSE, quote = TRUE)

# Check that reducing to the selected interaction's required genes preserves
# the full-input CellChat result for a multimeric receptor with modifier fields.
validation_id <- selected_ids[['multimeric_receptor']]
validation_selected <- pairs[validation_id, , drop = FALSE]
validation_full <- base
validation_full@LR$LRsig <- validation_selected
validation_reduced <- reduce_to_lr_features(validation_full, validation_selected)
validation_original_full <- CellChat::computeCommunProb(validation_full,
  LR.use = validation_selected, type = 'triMean', raw.use = TRUE,
  population.size = FALSE, nboot = 100L, seed.use = 1L, Kh = 0.5, n = 1)
validation_original_reduced <- CellChat::computeCommunProb(validation_reduced$object,
  LR.use = validation_selected, type = 'triMean', raw.use = TRUE,
  population.size = FALSE, nboot = 100L, seed.use = 1L, Kh = 0.5, n = 1)
validation <- data.frame(structure = 'multimeric_receptor_with_modifiers',
  lr_index = validation_id, lr_name = rownames(validation_selected),
  n_genes_full = nrow(base@data.signaling), n_genes_reduced = validation_reduced$n_genes,
  normalization_max = validation_reduced$normalizer,
  normalization_anchor_gene = validation_reduced$anchor_gene,
  probability_max_abs_diff = max(abs(validation_original_full@net$prob -
    validation_original_reduced@net$prob), na.rm = TRUE),
  pvalue_max_abs_diff = max(abs(validation_original_full@net$pval -
    validation_original_reduced@net$pval), na.rm = TRUE),
  probability_exact = sum(validation_original_full@net$prob !=
    validation_original_reduced@net$prob, na.rm = TRUE) == 0L,
  pvalue_exact = sum(validation_original_full@net$pval !=
    validation_original_reduced@net$pval, na.rm = TRUE) == 0L)
write.table(validation, file.path(out_dir, 'feature_reduction_validation.tsv'),
  sep = '\t', row.names = FALSE, quote = TRUE)

run_case <- function(structure, state, selected_id, zero_field = '', zero_mode = '',
                     storage = 'sparse', prime_expression = FALSE) {
  selected <- pairs[selected_id, , drop = FALSE]
  obj <- base
  obj@LR$LRsig <- selected
  if (prime_expression) {
    sender_cells <- which(obj@idents == levels(obj@idents)[[1L]])
    receiver_cells <- which(obj@idents == levels(obj@idents)[[2L]])
    ligand_genes <- component_genes(selected$ligand[[1L]])
    receptor_genes <- component_genes(selected$receptor[[1L]])
    if (!length(sender_cells) || !length(receiver_cells) || !length(ligand_genes) || !length(receptor_genes))
      stop('Cannot synthesize a nonzero sender/receiver expression state for ', structure)
    obj@data.signaling[ligand_genes, sender_cells] <- 1
    obj@data.signaling[receptor_genes, receiver_cells] <- 1
  }
  zeroed <- character()
  if (state %in% c('zero_required_component', 'zero_ligand_with_modifier')) {
    genes <- component_genes(selected[[zero_field]][[1L]])
    if (zero_mode == 'one_component') genes <- head(genes, 1L)
    if (!length(genes)) stop('No measured component available for ', structure)
    obj@data.signaling[genes, ] <- 0
    zeroed <- genes
  }
  reduced <- reduce_to_lr_features(obj, selected)
  obj <- reduced$object
  if (storage == 'dense') obj@data.signaling <- as.matrix(obj@data.signaling)
  original_start <- proc.time()[['elapsed']]
  original <- CellChat::computeCommunProb(obj, LR.use = selected, type = 'triMean', raw.use = TRUE,
    population.size = FALSE, nboot = 100L, seed.use = 1L, Kh = 0.5, n = 1)
  original_sec <- proc.time()[['elapsed']] - original_start
  ess_start <- proc.time()[['elapsed']]
  ess <- CellChatESS::computeCommunProbESS(obj, type = 'triMean', raw.use = TRUE,
    population.size = FALSE, nboot = 100L, seed.use = 1L, Kh = 0.5, n = 1,
    algorithm = 'ess', threads = 1L)
  ess_sec <- proc.time()[['elapsed']] - ess_start
  original_count <- CellChat::aggregateNet(original)@net$count
  ess_count <- CellChatESS::aggregateNetESS(ess)@net$count
  prob_diff <- max(abs(ess@net$prob - original@net$prob), na.rm = TRUE)
  pval_diff <- max(abs(ess@net$pval - original@net$pval), na.rm = TRUE)
  data.frame(structure = structure, state = state, storage = storage,
    lr_index = selected_id, lr_name = rownames(selected), ligand = as.character(selected$ligand),
    receptor = as.character(selected$receptor),
    co_A_receptor = as.character(selected$co_A_receptor),
    co_I_receptor = as.character(selected$co_I_receptor),
    agonist = as.character(selected$agonist), antagonist = as.character(selected$antagonist),
    zeroed_field = if (state %in% c('zero_required_component', 'zero_ligand_with_modifier')) zero_field else '',
    zero_mode = if (state %in% c('zero_required_component', 'zero_ligand_with_modifier')) zero_mode else '',
    zeroed_components = paste(zeroed, collapse = ';'),
    cells = ncol(obj@data.signaling), groups_K = nlevels(obj@idents), n_LR = nrow(selected),
    full_signaling_genes = nrow(base@data.signaling),
    projected_signaling_genes = reduced$n_genes,
    normalization_max = reduced$normalizer,
    normalization_anchor_gene = reduced$anchor_gene,
    nboot = 100L, original_sec = original_sec, ess_sec = ess_sec,
    probability_max_abs_diff = prob_diff, pvalue_max_abs_diff = pval_diff,
    probability_different_elements = sum(ess@net$prob != original@net$prob),
    pvalue_different_elements = sum(ess@net$pval != original@net$pval),
    count_different_elements = sum(ess_count != original_count),
    count_exact = all(ess_count == original_count),
    dimensions_and_names_match = identical(dim(ess@net$prob),dim(original@net$prob)) && identical(dimnames(ess@net$prob),dimnames(original@net$prob)),
    CellChat_version = as.character(packageVersion('CellChat')),
    CellChatESS_version = as.character(packageVersion('CellChatESS')),
    probability_exact = sum(ess@net$prob != original@net$prob, na.rm = TRUE) == 0L,
    pvalue_exact = sum(ess@net$pval != original@net$pval, na.rm = TRUE) == 0L,
    probability_finite = all(is.finite(ess@net$prob)), pvalue_finite = all(is.finite(ess@net$pval)),
    synthetic_expression_priming = prime_expression,
    probability_max_original = max(original@net$prob, na.rm = TRUE),
    pvalue_min_original = min(original@net$pval, na.rm = TRUE),
    stringsAsFactors = FALSE)
}

results <- list()
append <- function(x) results[[length(results) + 1L]] <<- x
for (structure in names(criteria)) {
  spec <- criteria[[structure]]
  id <- selected_ids[[structure]]
  prime <- synthetic_ids[[structure]]
  append(run_case(structure, 'nonzero_components', id, spec$zero_field, spec$zero_mode, 'sparse', prime))
  append(run_case(structure, 'zero_required_component', id, spec$zero_field, spec$zero_mode, 'sparse', prime))
  if (structure %in% c('co_activator', 'co_inhibitor', 'agonist', 'antagonist')) {
    # Keep the annotated modifier in the selected LR row and remove the base
    # ligand support, testing that modifiers cannot create a signal from zero.
    append(run_case(structure, 'zero_ligand_with_modifier', id,
      'ligand', 'all_components', 'sparse', prime))
  }
}
# Compare dense and sparse storage for representative single and multimeric
# receptor interactions, with both required components present and zeroed.
for (structure in c('single_ligand_receptor', 'multimeric_receptor')) {
  spec <- criteria[[structure]]
  id <- selected_ids[[structure]]
  for (storage in c('sparse', 'dense')) {
    append(run_case(paste0(structure, '_storage_check'), 'nonzero_components', id,
      spec$zero_field, spec$zero_mode, storage, synthetic_ids[[structure]]))
    append(run_case(paste0(structure, '_storage_check'), 'zero_required_component', id,
      spec$zero_field, spec$zero_mode, storage, synthetic_ids[[structure]]))
  }
}
out <- do.call(rbind, results)
write.table(out, file.path(out_dir, 'support_pruning_matrix_results.tsv'), sep = '\t',
  row.names = FALSE, quote = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, 'sessionInfo.txt'))
cat('completed ', nrow(out), ' cases; all probability differences < 1e-12: ',
  all(out$probability_max_abs_diff < 1e-12), '; all p-values exact: ', all(out$pvalue_exact), '\n', sep = '')

stopifnot(nrow(out)>=26L,all(out$probability_max_abs_diff<1e-12),all(out$pvalue_exact),all(out$count_exact),all(out$dimensions_and_names_match),all(out$probability_finite),all(out$pvalue_finite))
