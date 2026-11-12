args <- commandArgs(TRUE)
outdir <- args[[1]]
root <- args[[2]]
project <- args[[3]]
.libPaths(c(Sys.getenv('CELLCHAT_MAIN_R_LIB','environment/main/library'), .libPaths()))
suppressPackageStartupMessages({library(CellChat); library(CellChatESS)})
stopifnot(as.character(packageVersion('CellChat'))=='1.6.1',
          as.character(packageVersion('CellChatESS'))=='0.1.1')
future::plan('sequential')
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
tab <- read.delim(file.path(root,'common_external_timer_216.tsv'),check.names=FALSE)
cmp <- function(a,b,tag) data.frame(comparison=tag,
  max_abs_prob_diff=max(abs(a$prob-b$prob)),
  probability_unequal_elements=sum(a$prob!=b$prob),
  pvalue_unequal_elements=sum(a$pval!=b$pval),
  max_abs_pvalue_diff=max(abs(a$pval-b$pval)))
# All input groups have >10 cells (checked against the audited characteristics).
# Thus filterCommunication(min.cells=10) does not change these arrays.
# This exactly follows aggregateNet's default count rule, thresh=0.05.
count_from_net <- function(net) {
  p <- net$prob; q <- net$pval; q[p==0] <- 1; p[q>=0.05] <- 0
  z <- apply(p>0,c(1,2),sum); z[is.na(z)] <- 0; z
}
if (!file.exists(file.path(outdir,'P1_accuracy_216.tsv'))) {
rows <- vector('list',nrow(tab))
for(i in seq_len(nrow(tab))) {
  id <- tab$experiment_id[i]
  a <- readRDS(file.path(root,'pairs',id,'original.rds'))
  b <- readRDS(file.path(root,'pairs',id,'ess.rds'))
  r <- cmp(a,b,id); ca <- count_from_net(a); cb <- count_from_net(b)
  r$aggregate_count_unequal_elements <- sum(ca!=cb)
  r$aggregate_count_exact <- identical(ca,cb)
  r$count_definition <- 'aggregateNet default: sum(prob>0 after pval>=0.05 removal); min group cells>10'
  r$original_source <- file.path('pairs',id,'original.rds')
  r$ess_source <- file.path('pairs',id,'ess.rds')
  rows[[i]] <- r
  rm(a,b,ca,cb); if(i%%36==0) {cat('P1',i,'/216\n');invisible(gc())}
}
res <- do.call(rbind,rows)
write.table(res,file.path(outdir,'P1_accuracy_216.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
summary <- data.frame(n_pairs=nrow(res),max_abs_prob_diff=max(res$max_abs_prob_diff),
  pvalue_unequal_elements=sum(res$pvalue_unequal_elements),
  aggregate_count_exact_pairs=sum(res$aggregate_count_exact),
  CellChat='1.6.1',CellChatESS='0.1.1')
write.table(summary,file.path(outdir,'P1_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
} else {
summary <- read.delim(file.path(outdir,'P1_summary.tsv'))
}

# P2: reuse exactly the data/LR/idents held by the saved K=9 representative.
b <- readRDS(file.path(project,'data/representative/clean_baseline.rds'))
s <- readRDS(file.path(project,'data/representative/clean_sparse_stream.rds'))
stopifnot(nlevels(b@idents)==9,ncol(b@data.signaling)==1000)
input <- b; input@net <- list(); input@netP <- list()
saveRDS(input,file.path(outdir,'representative_prepared_K9.rds'))
pars <- b@options$parameter
invisible(gc()); t0 <- proc.time()[['elapsed']]
e <- CellChatESS::computeCommunProbESS(input,algorithm='ess',type=pars$type.mean,
  raw.use=pars$raw.use,population.size=pars$population.size,
  nboot=pars$nboot,seed.use=pars$seed.use,Kh=pars$Kh,n=pars$n)
elapsed <- proc.time()[['elapsed']]-t0
e <- CellChat::filterCommunication(e,min.cells=10)
e <- CellChat::computeCommunProbPathway(e)
e <- CellChat::aggregateNet(e)
saveRDS(e,file.path(outdir,'representative_ESS_0.1.1_K9.rds'))
rr <- do.call(rbind,lapply(list(list(b,'new_ESS_vs_original'),list(s,'new_ESS_vs_saved_sparse_stream')),function(z){
  r <- cmp(e@net,z[[1]]@net,z[[2]])
  r$count_unequal_elements <- sum(e@net$count!=z[[1]]@net$count)
  r$count_elementwise_exact <- identical(e@net$count,z[[1]]@net$count)
  r
}))
rr$inference_call_sec <- elapsed;rr$CellChat <- '1.6.1';rr$CellChatESS <- '0.1.1'
rr$cells <- ncol(input@data.signaling);rr$K <- nlevels(input@idents);rr$nLR <- nrow(input@LR$LRsig)
write.table(rr,file.path(outdir,'P2_representative_accuracy.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(outdir,'P1_P2_sessionInfo.txt'))
writeLines(capture.output(print(CellChat::aggregateNet)),file.path(outdir,'aggregateNet_count_source.R'))
print(summary);print(rr)
