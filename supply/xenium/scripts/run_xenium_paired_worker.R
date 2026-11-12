#!/usr/bin/env Rscript
# Same-input, same-reader/writer full-Xenium worker for Original CellChat or ESS.
# Inference-call timing excludes process/package/input setup, aggregation and saveRDS.
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 9L) stop(paste(
  'usage: METHOD INPUT INPUT_SHA OUTPUT_DIR CELLCHAT_LIB ESS_LIB DEP_LIB',
  'CELLCHAT_MODELING_R ESS_ACCEL_R'))
method <- match.arg(a[[1L]], c('original', 'ess'))
input <- normalizePath(a[[2L]], mustWork = TRUE)
expected_input_sha <- a[[3L]]
out <- normalizePath(a[[4L]], mustWork = TRUE)
cellchat_lib <- normalizePath(a[[5L]], mustWork = TRUE)
ess_lib <- normalizePath(a[[6L]], mustWork = TRUE)
dep_lib <- normalizePath(a[[7L]], mustWork = TRUE)
cellchat_source <- normalizePath(a[[8L]], mustWork = TRUE)
ess_source <- normalizePath(a[[9L]], mustWork = TRUE)
expected_ess_dll_sha <- Sys.getenv('EXPECTED_ESS_DLL_SHA256', '')
if (!nzchar(expected_ess_dll_sha)) stop('EXPECTED_ESS_DLL_SHA256 must be set')
if (file.exists(file.path(out, 'result.rds'))) stop('Refusing to overwrite result.rds')

Sys.setenv(OMP_NUM_THREADS = '1', OPENBLAS_NUM_THREADS = '1', MKL_NUM_THREADS = '1',
  BLIS_NUM_THREADS = '1', VECLIB_MAXIMUM_THREADS = '1', RCPP_PARALLEL_NUM_THREADS = '1',
  AUTOZYME_DISABLED = 'true', AUTOZYME_THREADS = '1', ZYME_THREADS = '1')
base_libs <- unique(c(.Library.site, .Library))
.libPaths(unique(c(if (method == 'ess') ess_lib, cellchat_lib, dep_lib, base_libs)))
suppressPackageStartupMessages(library(CellChat))
if (method == 'ess') suppressPackageStartupMessages(library(CellChatESS))

sha_file <- function(path) {
  z <- system2('sha256sum', shQuote(path), stdout = TRUE, stderr = TRUE)
  if (!length(z) || !grepl('^[[:xdigit:]]{64}[[:space:]]', z[[1L]])) stop('sha256sum failed: ', path)
  sub('[[:space:]].*$', '', z[[1L]])
}
hash_text <- function(x) digest::digest(paste(x, collapse = '\n'), algo = 'sha256', serialize = FALSE)
function_source_expression <- function(path, symbol) {
  ex <- parse(path, keep.source = FALSE)
  found <- Filter(function(e) is.call(e) && identical(e[[1L]], as.name('<-')) &&
    identical(e[[2L]], as.name(symbol)) && is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name('function')), ex)
  if (length(found) != 1L) stop('Could not uniquely find ', symbol, ' in ', path)
  found[[1L]][[3L]]
}
fun_text <- function(f) paste(deparse(f, width.cutoff = 500L), collapse = '\n')
actual_pkg <- function(pkg) normalizePath(find.package(pkg), mustWork = TRUE)
if (!identical(sha_file(input), expected_input_sha)) stop('Input SHA256 mismatch')
if (!requireNamespace('digest', quietly = TRUE)) stop('digest package is required for provenance hashes')
stopifnot(as.character(packageVersion('CellChat')) == '2.2.0.9001',
  identical(actual_pkg('CellChat'), normalizePath(file.path(cellchat_lib, 'CellChat'))),
  !isNamespaceLoaded('autozyme'))
if (method == 'ess') stopifnot(as.character(packageVersion('CellChatESS')) == '0.2.2',
  identical(actual_pkg('CellChatESS'), normalizePath(file.path(ess_lib, 'CellChatESS'))))

cellchat_ns <- asNamespace('CellChat')
cellchat_compute <- get('computeCommunProb', envir = cellchat_ns, inherits = FALSE)
stopifnot(identical(cellchat_compute, getExportedValue('CellChat', 'computeCommunProb')),
  identical(environmentName(environment(cellchat_compute)), 'CellChat'))
cellchat_expr <- function_source_expression(cellchat_source, 'computeCommunProb')
stopifnot(identical(body(cellchat_compute), body(eval(cellchat_expr, envir = new.env(parent = baseenv())))))
upstream_env <- new.env(parent = cellchat_ns)
sys.source(cellchat_source, envir = upstream_env)
upstream_symbols <- c('computeCommunProb', 'triMean', 'computeExpr_LR',
  'computeExpr_complex', 'computeExpr_coreceptor', 'computeExpr_agonist',
  'computeExpr_antagonist', 'aggregateNet')
upstream_checks <- do.call(rbind, lapply(upstream_symbols, function(symbol) {
  pinned <- get(symbol, envir = upstream_env, inherits = FALSE)
  loaded <- get(symbol, envir = cellchat_ns, inherits = FALSE)
  data.frame(symbol = symbol,
    formals_match = identical(formals(pinned), formals(loaded)),
    body_matches = identical(deparse(body(pinned)), deparse(body(loaded))),
    loaded_body_sha256 = hash_text(fun_text(loaded)), stringsAsFactors = FALSE)
}))
write.table(upstream_checks, file.path(out, 'upstream_function_identity.tsv'),
  sep = '\t', row.names = FALSE, quote = TRUE)
stopifnot(all(upstream_checks$formals_match), all(upstream_checks$body_matches))
active_fun <- if (method == 'original') cellchat_compute else {
  ess_fun <- get('computeCommunProbESS', envir = asNamespace('CellChatESS'), inherits = FALSE)
  stopifnot(identical(ess_fun, getExportedValue('CellChatESS', 'computeCommunProbESS')),
    identical(environmentName(environment(ess_fun)), 'CellChatESS'))
  ess_expr <- function_source_expression(ess_source, 'computeCommunProbESS')
  stopifnot(identical(body(ess_fun), body(eval(ess_expr, envir = new.env(parent = baseenv())))))
  ess_fun
}

stage_file <- file.path(out, 'stage_events.tsv')
if (file.exists(stage_file)) stop('Refusing to overwrite stage events')
writeLines('timestamp_epoch_ms\tstage\tdetail', stage_file)
mark_stage <- function(stage, detail = '') {
  row <- paste0(sprintf('%.3f', as.numeric(Sys.time()) * 1000), '\t', stage, '\t', gsub('[\t\r\n]', ' ', detail))
  cat(row, '\n', file = stage_file, append = TRUE, sep = '')
  flush.console()
  cat('[STAGE] ', stage, ' ', detail, '\n', sep = '')
}
future::plan('sequential')
mark_stage('process_initialized', paste0('method=', method))
obj <- readRDS(input)
stopifnot(identical(obj@options$datatype, 'RNA'))
stopifnot(identical(obj@idents, droplevels(obj@idents)))
stopifnot(is.factor(obj@idents), nlevels(obj@idents) == length(unique(obj@idents)))
stopifnot(length(obj@idents) == ncol(obj@data.signaling),
  !is.null(colnames(obj@data.signaling)), !anyDuplicated(colnames(obj@data.signaling)),
  is.null(names(obj@idents)) || identical(names(obj@idents), colnames(obj@data.signaling)))
if (identical(Sys.getenv('REQUIRE_XENIUM_FULL', '1'), '1'))
  stopifnot(ncol(obj@data.signaling) == 1156091L, nrow(obj@data.signaling) == 648L,
    nlevels(obj@idents) == 2357L, nrow(obj@LR$LRsig) == 627L)
expected_K <- nlevels(obj@idents)
expected_L <- nrow(obj@LR$LRsig)

params <- list(type = 'triMean', raw.use = TRUE, population.size = FALSE,
  nboot = 100L, seed.use = 1L, Kh = 0.5, n = 1, threads = 1L,
  algorithm = if (method == 'ess') 'ess' else NA_character_, future_workers = 1L)
input_info <- list(input_path = input, input_sha256 = sha_file(input),
  cells = ncol(obj@data.signaling), signaling_genes = nrow(obj@data.signaling),
  groups = nlevels(obj@idents), n_LR = nrow(obj@LR$LRsig),
  expression_matrix_file_bytes = unname(file.info(input)$size),
  expression_matrix_object_bytes = as.numeric(object.size(obj@data.signaling)),
  prepared_CellChat_object_bytes = as.numeric(object.size(obj)),
  cell_identifiers_and_group_labels_aligned = TRUE,
  database_serialized_sha256 = digest::digest(obj@DB, algo = 'sha256'),
  LRsig_serialized_sha256 = digest::digest(obj@LR$LRsig, algo = 'sha256'),
  cell_names_sha256 = digest::digest(colnames(obj@data.signaling), algo = 'sha256'),
  cell_group_labels_sha256 = digest::digest(as.character(obj@idents), algo = 'sha256'),
  group_order_sha256 = digest::digest(levels(obj@idents), algo = 'sha256'),
  LR_order_sha256 = digest::digest(rownames(obj@LR$LRsig), algo = 'sha256'),
  group_order = levels(obj@idents), LR_order = rownames(obj@LR$LRsig))

capture_fun <- if (method == 'original') cellchat_compute else active_fun
fun_text_value <- fun_text(capture_fun)
function_info <- list(function_name = if (method == 'original') 'CellChat::computeCommunProb' else
    'CellChatESS::computeCommunProbESS', function_environment = environmentName(environment(capture_fun)),
  function_body_sha256 = hash_text(fun_text_value),
  function_body_matches_pinned_source = TRUE,
  CellChat_namespace_function_body_sha256 = hash_text(fun_text(cellchat_compute)),
  CellChat_source_file = cellchat_source, CellChat_source_sha256 = sha_file(cellchat_source),
  ESS_source_file = if (method == 'ess') ess_source else NA_character_,
  ESS_source_sha256 = if (method == 'ess') sha_file(ess_source) else NA_character_,
  AutoZyme_namespace_loaded = isNamespaceLoaded('autozyme'),
  CellChat_function_export_matches_namespace_binding = identical(cellchat_compute,
    getExportedValue('CellChat', 'computeCommunProb')))
if (method == 'ess') {
  dll <- getLoadedDLLs()[['CellChatESS']]
  if (is.null(dll)) stop('CellChatESS shared library is not loaded')
  dll_path <- normalizePath(dll[['path']], mustWork = TRUE)
  dll_sha <- sha_file(dll_path)
  stopifnot(identical(dll_sha, expected_ess_dll_sha))
  function_info$ESS_DLL_path <- dll_path
  function_info$ESS_DLL_sha256 <- dll_sha
  function_info$ESS_function_body_matches_pinned_source <- TRUE
}
writeLines(fun_text_value, file.path(out, 'active_compute_function.R'))

write.table(data.frame(group = levels(obj@idents), n_cells = as.integer(table(obj@idents)),
  stringsAsFactors = FALSE), file.path(out, 'group_sizes.tsv'), sep = '\t',
  row.names = FALSE, quote = TRUE)
write.table(obj@LR$LRsig, file.path(out, 'LRsig.tsv'), sep = '\t',
  row.names = TRUE, col.names = NA, quote = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo_before.txt'))
writeLines(capture.output(extSoftVersion()), file.path(out, 'external_software_versions.txt'))
writeLines(capture.output(Sys.getenv(c('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
  'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'RCPP_PARALLEL_NUM_THREADS'))),
  file.path(out, 'thread_environment.txt'))
identity <- list(method = method, package_versions = list(CellChat = as.character(packageVersion('CellChat')),
  CellChatESS = if (method == 'ess') as.character(packageVersion('CellChatESS')) else NA_character_),
  installation_paths = list(CellChat = actual_pkg('CellChat'),
    CellChatESS = if (method == 'ess') actual_pkg('CellChatESS') else NA_character_),
  source_function_identity = function_info, input = input_info, parameters = params,
  rng_kind_before_call = RNGkind(), blas_lapack = extSoftVersion()[c('BLAS', 'LAPACK')],
  r_version = R.version.string, platform = R.version$platform,
  libraries = .libPaths())
writeLines(jsonlite::toJSON(identity, auto_unbox = TRUE, pretty = TRUE, null = 'null'),
  file.path(out, 'run_identity_before.json'))
gc(verbose = FALSE)
mark_stage('before_inference', paste0('input_expression_object_bytes=', input_info$expression_matrix_object_bytes))
start_epoch <- as.numeric(Sys.time())
started <- proc.time()[['elapsed']]
if (method == 'original') {
  result <- CellChat::computeCommunProb(obj, type = 'triMean', raw.use = TRUE,
    population.size = FALSE, nboot = 100L, seed.use = 1L, Kh = 0.5, n = 1)
} else {
  result <- CellChatESS::computeCommunProbESS(obj, type = 'triMean', raw.use = TRUE,
    population.size = FALSE, nboot = 100L, seed.use = 1L, Kh = 0.5, n = 1,
    algorithm = 'ess', threads = 1L)
}
inference_sec <- proc.time()[['elapsed']] - started
end_epoch <- as.numeric(Sys.time())
random_seed_after_inference_sha256 <- digest::digest(.Random.seed, algo = 'sha256')
mark_stage('after_inference', paste0('function_elapsed_sec=', sprintf('%.9f', inference_sec)))
if (method == 'ess') stopifnot(identical(as.integer(result@options$CellChatESS$threads_used), 1L))
prob <- result@net$prob
pval <- result@net$pval
stopifnot(identical(dim(prob), dim(pval)), identical(dimnames(prob), dimnames(pval)),
  identical(as.integer(dim(prob)), c(expected_K, expected_K, expected_L)))
prob_object_bytes <- as.numeric(object.size(prob))
pval_object_bytes <- as.numeric(object.size(pval))
prob_payload_bytes <- 8 * length(prob)
pval_payload_bytes <- 8 * length(pval)
mark_stage('before_aggregateNet', paste0('prob_object_bytes=', prob_object_bytes,
  ';pval_object_bytes=', pval_object_bytes))
agg_start <- proc.time()[['elapsed']]
aggregated <- CellChat::aggregateNet(result, return.object = FALSE)
aggregate_sec <- proc.time()[['elapsed']] - agg_start
count <- aggregated$count
mark_stage('after_aggregateNet', paste0('aggregate_elapsed_sec=', sprintf('%.9f', aggregate_sec),
  ';count_bytes=', as.numeric(object.size(count))))
if (!identical(as.integer(dim(count)), c(expected_K, expected_K))) stop('Unexpected aggregate count shape')

result_payload <- list(prob = prob, pval = pval, count = count,
  groups = dimnames(prob)[[1L]], lr_names = dimnames(prob)[[3L]], parameters = params)
mark_stage('before_saveRDS', paste0('compression=gzip;output_fields=', paste(names(result_payload), collapse = ',')))
saveRDS(result_payload, file.path(out, 'result.rds'), compress = 'gzip', version = 3L)
mark_stage('after_saveRDS', paste0('output_bytes=', file.info(file.path(out, 'result.rds'))$size))
run_record <- list(method = method, status = 'COMPLETED', function_call_start_epoch = start_epoch,
  function_call_end_epoch = end_epoch, inference_function_elapsed_sec = inference_sec,
  aggregateNet_elapsed_sec = aggregate_sec, parameters = params,
  random_seed_after_inference_sha256 = random_seed_after_inference_sha256,
  actual_ESS_threads_used = if (method == 'ess') result@options$CellChatESS$threads_used else NA_integer_,
  input = input_info, outputs = list(dimensions = dim(prob), dimnames = dimnames(prob),
    probability_object_bytes = prob_object_bytes, pvalue_object_bytes = pval_object_bytes,
    probability_payload_bytes = prob_payload_bytes, pvalue_payload_bytes = pval_payload_bytes,
    count_object_bytes = as.numeric(object.size(count)), output_rds_bytes = file.info(file.path(out, 'result.rds'))$size,
    output_rds_sha256 = sha_file(file.path(out, 'result.rds'))),
  identity = identity)
writeLines(jsonlite::toJSON(run_record, auto_unbox = TRUE, pretty = TRUE, null = 'null'),
  file.path(out, 'run_metadata.json'))
write.table(data.frame(method = method, inference_function_elapsed_sec = inference_sec,
  aggregateNet_elapsed_sec = aggregate_sec, input_file_bytes = input_info$expression_matrix_file_bytes,
  input_expression_object_bytes = input_info$expression_matrix_object_bytes,
  probability_object_bytes = prob_object_bytes, pvalue_object_bytes = pval_object_bytes,
  probability_payload_bytes = prob_payload_bytes, pvalue_payload_bytes = pval_payload_bytes,
  count_object_bytes = as.numeric(object.size(count)), output_rds_bytes = file.info(file.path(out, 'result.rds'))$size,
  cells = ncol(obj@data.signaling), signaling_genes = nrow(obj@data.signaling),
  K = nlevels(obj@idents), n_LR = nrow(obj@LR$LRsig),
  nboot = 100L, seed_use = 1L, threads = 1L, stringsAsFactors = FALSE),
  file.path(out, 'measurement.tsv'), sep = '\t', row.names = FALSE, quote = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo_after.txt'))
mark_stage('worker_complete', 'all result files written')
cat('Completed ', method, '; inference_sec=', sprintf('%.9f', inference_sec),
  '; aggregateNet_sec=', sprintf('%.9f', aggregate_sec), '; output=', file.path(out, 'result.rds'), '\n', sep = '')
