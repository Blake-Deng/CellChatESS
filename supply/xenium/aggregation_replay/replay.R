args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6L)
prepared_file <- args[[1L]]
result_file <- args[[2L]]
method <- match.arg(args[[3L]], c("CellChat_2.2.0.9001", "CellChatESS_0.2.2"))
pair <- as.integer(args[[4L]])
ord <- as.integer(args[[5L]])
outfile <- args[[6L]]
suppressPackageStartupMessages(library(CellChat))
stopifnot(as.character(packageVersion("CellChat")) == "2.2.0.9001")
future::plan("sequential")
f <- CellChat::aggregateNet
fh <- digest::digest(paste(deparse(f, width.cutoff=500L), collapse="\n"),
                     algo="sha256", serialize=FALSE)
stopifnot(fh == "ee8e72f56c39c0d8ed6df89e2057477401bc0e126a97ab93492a012ccb815d3f")
message("Loading saved results for ", method, " in PID ", Sys.getpid(), " at ", Sys.time())
obj <- readRDS(prepared_file)
saved <- readRDS(result_file)
stopifnot(identical(dim(saved$prob), c(2357L, 2357L, 627L)),
          identical(dim(saved$prob), dim(saved$pval)),
          identical(dimnames(saved$prob), dimnames(saved$pval)),
          identical(as.character(dimnames(saved$prob)[[1L]]), as.character(levels(obj@idents))),
          identical(as.character(dimnames(saved$prob)[[3L]]), as.character(rownames(obj@LR$LRsig))))
expected_count <- saved$count
obj@net <- list(prob=saved$prob, pval=saved$pval)
rm(saved)
invisible(gc.time(TRUE))
invisible(gc(full=TRUE, reset=TRUE))
message("AGGREGATE_START ", Sys.time())
gc_before <- gc.time()
time_before <- proc.time()
net <- CellChat::aggregateNet(obj, return.object=FALSE)
cpu <- proc.time() - time_before
gc_delta <- gc.time() - gc_before
message("AGGREGATE_END ", Sys.time(), " elapsed=", unname(cpu[["elapsed"]]))
stopifnot(identical(net$count, expected_count), sum(net$count) == 70162984)
row <- data.frame(pair=pair, position=ord, method=method,
                  elapsed_s=unname(cpu[["elapsed"]]), user_s=unname(cpu[["user.self"]]),
                  system_s=unname(cpu[["sys.self"]]), gc_elapsed_s=gc_delta[[3L]],
                  count_matches_saved=TRUE, count_sha256=digest::digest(net$count, algo="sha256"),
                  total_count=sum(net$count), pid=Sys.getpid(), fn_sha256=fh)
write.table(row, file=paste0(outfile, ".tmp"), quote=FALSE, sep="\t", row.names=FALSE)
stopifnot(file.rename(paste0(outfile, ".tmp"), outfile))
