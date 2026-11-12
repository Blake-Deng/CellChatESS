args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: verify_v022_input_guard.R PREPARED_CELLCHAT_RDS OUTPUT_TSV")
}

.libPaths(c(Sys.getenv("CELLCHAT_ESS_R_LIB", "environment/compatible/library"), Sys.getenv("CELLCHAT_COMPAT_R_LIB", "environment/compatible/library"), Sys.getenv("CELLCHAT_MAIN_R_LIB", "environment/main/library"), .libPaths()))
stopifnot(as.character(packageVersion("CellChat")) == "2.2.0.9001")
suppressPackageStartupMessages({
  library(CellChat)
  library(CellChatESS)
  library(Matrix)
})
stopifnot(as.character(packageVersion("CellChatESS")) == "0.2.2")

object <- readRDS(args[[1L]])
pairs <- object@LR$LRsig
ligand <- as.character(pairs$ligand[[which(!is.na(pairs$ligand) & nzchar(pairs$ligand))[[1L]]]])
stopifnot(ligand %in% rownames(object@data.signaling))
non_lr <- setdiff(rownames(object@data.signaling),
                  unique(c(as.character(pairs$ligand), as.character(pairs$receptor))))
stopifnot(length(non_lr) > 0L)

cases <- list(
  negative_ligand = list(gene = ligand, value = -0.1, expected = "must be non-negative"),
  negative_non_lr = list(gene = non_lr[[1L]], value = -0.1, expected = "must be non-negative"),
  na_ligand = list(gene = ligand, value = NA_real_, expected = "only finite values"),
  nan_ligand = list(gene = ligand, value = NaN, expected = "only finite values"),
  inf_ligand = list(gene = ligand, value = Inf, expected = "only finite values")
)

records <- lapply(names(cases), function(case) {
  spec <- cases[[case]]
  test_object <- object
  test_object@data.signaling[spec$gene, 1L] <- spec$value
  error <- tryCatch({
    CellChatESS::computeCommunProbESS(test_object, nboot = 1L, seed.use = 1L)
    NULL
  }, error = identity)
  message <- if (inherits(error, "error")) conditionMessage(error) else ""
  passed <- inherits(error, "error") && grepl(spec$expected, message, fixed = TRUE)
  if (!passed) stop("Input guard test failed for ", case, ": ", message)
  data.frame(
    case = case,
    gene = spec$gene,
    expected_guard = spec$expected,
    status = "rejected_as_expected",
    error = message,
    CellChat_version = as.character(packageVersion("CellChat")),
    CellChatESS_version = as.character(packageVersion("CellChatESS")),
    stringsAsFactors = FALSE
  )
})

result <- do.call(rbind, records)
write.table(result, args[[2L]], sep = "\t", row.names = FALSE, quote = TRUE)
writeLines(capture.output(sessionInfo()), sub("\\.tsv$", "_sessionInfo.txt", args[[2L]]))
cat("v0.2.2 public API input-guard checks passed for", nrow(result), "cases.\n")
