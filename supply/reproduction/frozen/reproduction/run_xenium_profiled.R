args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop('usage: run_xenium_profiled.R INPUT_RDS OUTPUT_RDS SUMMARY_TSV SESSIONINFO')
input_path <- args[[1L]]
output_path <- args[[2L]]
summary_path <- args[[3L]]
session_path <- args[[4L]]
profile_file <- Sys.getenv('CELLCHATESS_STAGE_FILE', '')
if (!nzchar(profile_file)) stop('CELLCHATESS_STAGE_FILE must be set')

if(file.exists(profile_file) || file.exists(output_path) || file.exists(summary_path)) stop('Choose fresh profile output paths')
header <- data.frame(epoch_ms = 'epoch_ms', stage = 'stage', rss_kb = 'rss_kb',
  input_bytes = 'input_bytes', probability_bytes = 'probability_bytes',
  pvalue_bytes = 'pvalue_bytes', object_bytes = 'object_bytes', detail = 'detail')
write.table(header, profile_file, sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
profile_stage <- function(stage, input_bytes = NA_real_, prob_bytes = NA_real_,
                          pval_bytes = NA_real_, object_bytes = NA_real_, detail = '') {
  epoch_ms <- as.numeric(round(as.numeric(Sys.time()) * 1000))
  status <- tryCatch(readLines('/proc/self/status', warn = FALSE), error = function(e) character())
  rss <- grep('^VmRSS:', status, value = TRUE)
  rss_kb <- if (length(rss)) as.numeric(strsplit(trimws(sub('^VmRSS:', '', rss[[1L]])),
    '[[:space:]]+')[[1L]][[1L]]) else NA_real_
  event <- data.frame(epoch_ms = epoch_ms, stage = stage, rss_kb = rss_kb,
    input_bytes = input_bytes, probability_bytes = prob_bytes, pvalue_bytes = pval_bytes,
    object_bytes = object_bytes, detail = detail)
  write.table(event, profile_file, sep = '\t', row.names = FALSE, col.names = FALSE,
    quote = FALSE, append = TRUE, na = '')
  cat('[ESS_STAGE] ', stage, ' rss_kb=', rss_kb, ' detail=', detail, '\n', sep = '')
  flush.console()
}

cellchat22_lib <- Sys.getenv("CELLCHAT_COMPAT_R_LIB", "environment/compatible/library")
profile_lib <- Sys.getenv("CELLCHAT_ESS_R_LIB", "environment/compatible/library")
.libPaths(c(profile_lib, cellchat22_lib, .libPaths()))
Sys.setenv(OMP_NUM_THREADS = '1', OPENBLAS_NUM_THREADS = '1', MKL_NUM_THREADS = '1',
  RCPP_PARALLEL_NUM_THREADS = '1', BLIS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
suppressPackageStartupMessages({ library(CellChat); library(CellChatESS) })
stopifnot(as.character(packageVersion('CellChat')) == '2.2.0.9001',
          as.character(packageVersion('CellChatESS')) == '0.2.2')

stopifnot(!isNamespaceLoaded('autozyme'),as.character(utils::packageDescription('CellChatESS')[['Config/CellChatESS/Build']])=='PROFILE')
future::plan('sequential')
obj <- readRDS(input_path)
obj@idents <- droplevels(obj@idents)
if (!identical(obj@options$datatype, 'RNA')) stop("Prepared Xenium object must use datatype='RNA'")
gc(verbose = FALSE)
input_bytes <- as.numeric(object.size(obj@data.signaling))
profile_stage('R_input_loaded', input_bytes = input_bytes, object_bytes = as.numeric(object.size(obj)),
  detail = 'input expression matrix and prepared CellChat object loaded')
profile_stage('R_compute_enter', input_bytes = input_bytes, object_bytes = as.numeric(object.size(obj)),
  detail = 'immediately before computeCommunProbESS')
start <- proc.time()[['elapsed']]
result <- CellChatESS::computeCommunProbESS(obj, type = 'triMean', raw.use = TRUE,
  population.size = FALSE, nboot = 100L, seed.use = 1L, Kh = 0.5, n = 1,
  algorithm = 'ess', threads = 1L)
elapsed <- proc.time()[['elapsed']] - start
prob_bytes <- as.numeric(object.size(result@net$prob))
pval_bytes <- as.numeric(object.size(result@net$pval))
profile_stage('R_compute_return', input_bytes = input_bytes, prob_bytes = prob_bytes,
  pval_bytes = pval_bytes, object_bytes = as.numeric(object.size(result)),
  detail = paste0('inference returned; function_sec=', elapsed))
profile_stage('R_serialization_enter', input_bytes = input_bytes, prob_bytes = prob_bytes,
  pval_bytes = pval_bytes, object_bytes = as.numeric(object.size(result)),
  detail = 'before compressed saveRDS')
saveRDS(list(prob = result@net$prob, pval = result@net$pval,
  groups = dimnames(result@net$prob)[[1L]], lr_names = dimnames(result@net$prob)[[3L]],
  parameters = list(nboot = 100L, seed.use = 1L, Kh = 0.5, n = 1),
  engine = 'ess_profiled', versions = c(CellChat = as.character(packageVersion('CellChat')),
    CellChatESS = as.character(packageVersion('CellChatESS')))), output_path, compress = TRUE)
profile_stage('R_serialization_return', input_bytes = input_bytes, prob_bytes = prob_bytes,
  pval_bytes = pval_bytes, object_bytes = as.numeric(object.size(result)),
  detail = paste0('compressed result written; file_bytes=', file.info(output_path)$size))

summary <- data.frame(engine = 'CellChatESS_0.2.2_profiled', status = 'ok', elapsed_sec = elapsed,
  cells = ncol(obj@data.signaling), signaling_genes = nrow(obj@data.signaling),
  groups = nlevels(obj@idents), n_lr = nrow(obj@LR$LRsig), nboot = 100L, seed_use = 1L,
  input_bytes = input_bytes, probability_bytes = prob_bytes, pvalue_bytes = pval_bytes,
  probability_payload_bytes = 8*length(result@net$prob), pvalue_payload_bytes = 8*length(result@net$pval),
  threads_used = result@options$CellChatESS$threads_used,
  output_file_bytes = file.info(output_path)$size,
  CellChat_version = as.character(packageVersion('CellChat')),
  CellChatESS_version = as.character(packageVersion('CellChatESS')))
write.table(summary, summary_path, sep = '\t', row.names = FALSE, quote = TRUE)
writeLines(capture.output(sessionInfo()), session_path)
profile_stage('R_process_end', input_bytes = input_bytes, prob_bytes = prob_bytes,
  pval_bytes = pval_bytes, object_bytes = as.numeric(object.size(result)),
  detail = 'R runner completed all writes')
cat('profiled Xenium ESS completed; function_sec=', elapsed, ' output=', output_path, '\n', sep = '')
