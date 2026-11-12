#!/usr/bin/env Rscript
# Fresh-process, common external inference timer; comparison is a separate process.
a <- commandArgs(TRUE)
if (length(a) != 10L) stop('METHOD INPUT OUTDIR CELLCHAT_VERSION ESS_LIB CELLCHAT_LIB AUTO_LIB THREADS SOURCE_ID INPUT_HASH')
method<-a[1]; input<-a[2]; outdir<-a[3]; cc_version<-a[4]
threads<-as.integer(a[8]); stopifnot(method %in% c('original','ess','autozyme'), threads>=1L)
.libPaths(unique(c(a[5],a[6],a[7],.libPaths())))
Sys.setenv(OMP_NUM_THREADS=as.character(threads),OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
 BLIS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',RCPP_PARALLEL_NUM_THREADS='1',
 AUTOZYME_THREADS='1',ZYME_THREADS='1')
Sys.setenv(AUTOZYME_DISABLED=if(method=='autozyme')'false' else 'true')
suppressPackageStartupMessages(library(CellChat))
stopifnot(!isNamespaceLoaded('autozyme'))
stopifnot(as.character(packageVersion('CellChat'))==cc_version)
stopifnot(normalizePath(find.package('CellChat'))==normalizePath(file.path(a[6],'CellChat')))
if (method=='ess') {
 suppressPackageStartupMessages(library(CellChatESS))
 stopifnot(as.character(packageVersion('CellChatESS'))=='0.2.2')
 stopifnot(identical(readLines(system.file('ESS_SOURCE_ID',package='CellChatESS'),warn=FALSE),Sys.getenv('CELLCHATESS_EXPECTED_SOURCE_ID')))
 writeLines(readLines(system.file('ESS_SOURCE_ID',package='CellChatESS'),warn=FALSE),file.path(outdir,'ESS_SOURCE_ID.txt'))
 stopifnot(normalizePath(find.package('CellChatESS'))==normalizePath(file.path(a[5],'CellChatESS')))
}
if (method=='autozyme') {
 stopifnot(cc_version=='2.2.0.9001')
 suppressPackageStartupMessages(library(autozyme))
 stopifnot(as.character(packageVersion('autozyme'))=='0.3.1')
 stopifnot(normalizePath(find.package('autozyme'))==normalizePath(file.path(a[7],'autozyme')))
 autozyme::set_threads(1L);autozyme::activate('cellchat')
 stopifnot(identical(unname(autozyme::status()[['cellchat']]),'active'))
}
writeLines(deparse(CellChat::computeCommunProb),file.path(outdir,'actual_computeCommunProb.R'))
writeLines(deparse(CellChat::triMean),file.path(outdir,'actual_triMean.R'))
future::plan('sequential')
obj<-readRDS(input)
stopifnot(identical(obj@options$datatype,'RNA'),is.factor(obj@idents))
obj@idents<-droplevels(obj@idents)
writeLines(capture.output(sessionInfo()),file.path(outdir,'sessionInfo_before.txt'))
writeLines(capture.output(str(.libPaths())),file.path(outdir,'library_paths.txt'))
gc()
started<-proc.time()[['elapsed']]
result<-if (method=='ess') {
 CellChatESS::computeCommunProbESS(obj,type='triMean',raw.use=TRUE,population.size=FALSE,
  nboot=100L,seed.use=1L,Kh=0.5,n=1,algorithm='ess',threads=threads)
} else {
 CellChat::computeCommunProb(obj,type='triMean',raw.use=TRUE,population.size=FALSE,
  nboot=100L,seed.use=1L,Kh=0.5,n=1)
}
elapsed<-proc.time()[['elapsed']]-started
prob<-result@net$prob;pval<-result@net$pval
stopifnot(identical(dim(prob),dim(pval)),all(is.finite(prob)),all(is.finite(pval)))
# Same downstream count definition as CellChat aggregateNet default, outside inference timer.
filtered<-prob;filtered[pval>=0.05]<-0
count<-apply(filtered>0,c(1,2),sum)
rm(filtered)
metadata<-if(method=='ess') result@options$CellChatESS else list()
used<-if(is.null(metadata$threads_used)) NA_integer_ else as.integer(metadata$threads_used)
saveRDS(list(prob=prob,pval=pval,count=count,ess_metadata=metadata),file.path(outdir,'result.rds'),compress=FALSE)
getver<-function(p) if(p %in% loadedNamespaces()) as.character(packageVersion(p)) else NA_character_
row<-data.frame(method=method,inference_call_sec=elapsed,cells=ncol(obj@data.signaling),
 signaling_genes=nrow(obj@data.signaling),K=nlevels(obj@idents),nLR=nrow(obj@LR$LRsig),
 nboot=100L,seed_use=1L,threads_requested=threads,threads_used=used,
 CellChat_version=getver('CellChat'),ESS_loaded_version=getver('CellChatESS'),AutoZyme_loaded_version=getver('autozyme'),
 candidate_source_id=a[9],input_sha256=a[10],input_file=normalizePath(input),
 probability_payload_bytes=8*length(prob),pvalue_payload_bytes=8*length(pval),
 probability_R_object_bytes=as.numeric(object.size(prob)),pvalue_R_object_bytes=as.numeric(object.size(pval)),
 timing_scope='external inference function call; excludes startup/load/validation/output/downstream count',
 stringsAsFactors=FALSE)
write.table(row,file.path(outdir,'measurement.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
writeLines(capture.output(sessionInfo()),file.path(outdir,'sessionInfo_after.txt'))
dlls<-getLoadedDLLs()
dllpaths<-vapply(dlls,function(x)x[['path']],character(1))
write.table(data.frame(path=dllpaths,md5=unname(tools::md5sum(dllpaths))),
 file.path(outdir,'loaded_dlls.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
