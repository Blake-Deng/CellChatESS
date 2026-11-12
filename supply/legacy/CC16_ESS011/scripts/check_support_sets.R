#!/usr/bin/env Rscript
args <- commandArgs(TRUE);out <- normalizePath(args[[1]])
.libPaths(c(Sys.getenv("CELLCHAT_MAIN_R_LIB", "environment/main/library"),.libPaths()))
suppressPackageStartupMessages({library(CellChat);library(CellChatESS);library(Rcpp)})
stopifnot(as.character(packageVersion('CellChat'))=='1.6.1',as.character(packageVersion('CellChatESS'))=='0.1.1')
cppenv <- new.env()
sourceCpp('scripts/support_observed_only.cpp',env=cppenv)
env <- new.env(parent=asNamespace('CellChatESS'))
audit_result <- NULL
env$cellchat_prob_ess_cpp <- function(...) {audit_result <<- cppenv$cellchat_prob_ess_cpp(...); audit_result}
fun <- CellChatESS::computeCommunProbESS;environment(fun)<-env
meta <- read.delim('results/benchmark_216/common_external_timer_216.tsv')
prefix <- 'data/inputs'
rows <- vector('list',nrow(meta));bad <- list()
for(i in seq_len(nrow(meta))) {
  id <- meta$experiment_id[i]
  prepared <- file.path(prefix,id,'prepared_cellchat.rds')
  x <- readRDS(prepared);x@idents <- droplevels(x@idents)
  # nboot=1 only constructs a dummy label matrix; instrumented C++ skips bootstrap.
  invisible(fun(x,algorithm='ess',nboot=1,seed.use=1))
  a <- audit_result$positive_support;b <- audit_result$prob!=0
  saved <- file.path('data/timing_outputs',id,'ess.rds');s <- readRDS(saved)
  stopifnot(identical(dim(a),dim(b)),identical(unname(b),unname(s$prob!=0)),max(abs(audit_result$prob-s$prob))==0)
  rows[[i]] <- data.frame(experiment_id=id,dataset_id=meta$dataset_id[i],target_cells=meta$target_cells[i],positive_L_Radj_triples=sum(a),nonzero_Pobs_triples=sum(b),positive_but_zero_Pobs=sum(a & !b),nonzero_without_positive_support=sum(b & !a),sets_equal=identical(a,b),observed_prob_exact_saved=identical(audit_result$prob,unname(s$prob)),prepared_path=prepared,saved_output=saved)
  if(any(a!=b)) {
    ix<-which(a!=b,arr.ind=TRUE)
    bad[[length(bad)+1L]]<-data.frame(experiment_id=id,sender=ix[,1],receiver=ix[,2],LR=ix[,3],positive_support=a[ix],Pobs=audit_result$prob[ix])
  }
  write.table(do.call(rbind,rows[seq_len(i)]),file.path(out,'F5_support_set_counts_216.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  if(i%%12==0) {cat(sprintf('%d / %d inputs checked; mismatching inputs=%d\n',i,nrow(meta),length(bad)));flush.console()}
  rm(x,s);gc(FALSE)
}
if(length(bad))write.table(do.call(rbind,bad),file.path(out,'F5_mismatching_entries.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'F5_sessionInfo.txt'))
cat('Completed: observed support only, no bootstrap inference.\n')
