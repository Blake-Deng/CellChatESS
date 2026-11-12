a <- commandArgs(TRUE)
if(length(a)!=9L) stop('METHOD INPUT OUTDIR ESS_LIB CC_LIB AUTO_LIB DEPS_LIB CC_VERSION RAW_USE')
method<-a[1];stopifnot(method %in% c('original','autozyme','ess','dense'))
dir.create(a[3],recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(a[3],'result.rds')))stop('Choose a fresh output directory')
.libPaths(unique(c(a[5],a[4],a[6],a[7],.libPaths())))
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',
 VECLIB_MAXIMUM_THREADS='1',RCPP_PARALLEL_NUM_THREADS='1',AUTOZYME_THREADS='1',ZYME_THREADS='1',
 AUTOZYME_DISABLED=if(method=='autozyme')'false' else 'true')
suppressPackageStartupMessages(library(CellChat))
stopifnot(as.character(packageVersion('CellChat'))==a[8],
 normalizePath(find.package('CellChat'))==normalizePath(file.path(a[5],'CellChat')),
 !isNamespaceLoaded('autozyme'))
writeLines(deparse(CellChat::computeCommunProb),file.path(a[3],'CellChat_computeCommunProb_before.R'))
writeLines(deparse(CellChat::triMean),file.path(a[3],'CellChat_triMean_before.R'))
if(method=='autozyme') {
 suppressPackageStartupMessages(library(autozyme))
 stopifnot(as.character(packageVersion('autozyme'))=='0.3.1',normalizePath(find.package('autozyme'))==normalizePath(file.path(a[6],'autozyme')))
 autozyme::set_threads(1L);autozyme::activate('cellchat')
 stopifnot(identical(unname(autozyme::status()[['cellchat']]),'active'))
 writeLines(capture.output(autozyme::status()),file.path(a[3],'autozyme_status.txt'))
}
if(method %in% c('ess','dense')) {
 suppressPackageStartupMessages(library(CellChatESS))
 stopifnot(as.character(packageVersion('CellChatESS'))=='0.2.2',normalizePath(find.package('CellChatESS'))==normalizePath(file.path(a[4],'CellChatESS')))
}
writeLines(deparse(CellChat::computeCommunProb),file.path(a[3],'CellChat_computeCommunProb_active.R'))
writeLines(deparse(CellChat::triMean),file.path(a[3],'CellChat_triMean_active.R'))
future::plan('sequential')
x<-readRDS(a[2]);x@idents<-droplevels(x@idents)
settings<-list(type='triMean',raw.use=identical(a[9],'TRUE'),population.size=FALSE,nboot=100L,seed.use=1L,Kh=.5,n=1)
z<-if(method %in% c('ess','dense')) do.call(CellChatESS::computeCommunProbESS,c(list(object=x,threads=1L,algorithm=if(method=='dense')'dense' else 'ess'),settings)) else do.call(CellChat::computeCommunProb,c(list(object=x),settings))
rng<-.Random.seed
# Keep inference and downstream checks separate; both downstream paths see the same object.
filtered<-CellChat::filterCommunication(z,min.cells=10)
p<-if(method %in% c('ess','dense')) CellChatESS::computeCommunProbPathwayESS(filtered) else CellChat::computeCommunProbPathway(filtered)
nw<-if(method %in% c('ess','dense')) CellChatESS::aggregateNetESS(p) else CellChat::aggregateNet(p)
refp<-CellChat::computeCommunProbPathway(filtered);refnw<-CellChat::aggregateNet(refp)
res<-list(prob=z@net$prob,pval=z@net$pval,count=nw@net$count,weight=nw@net$weight,
 pathway_prob=p@netP$prob,pathways=p@netP$pathways,LRs=p@net$LRs,rng=rng,options=z@options,
 native_reference=list(count=refnw@net$count,weight=refnw@net$weight,pathway_prob=refp@netP$prob,pathways=refp@netP$pathways,LRs=refp@net$LRs))
saveRDS(res,file.path(a[3],'result.rds'))
writeLines(capture.output(sessionInfo()),file.path(a[3],'sessionInfo.txt'))
writeLines(capture.output(.libPaths()),file.path(a[3],'library_paths.txt'))
write.table(data.frame(method=method,CellChat=a[8],CellChatESS=if(method %in% c('ess','dense'))'0.2.2' else NA,
 autozyme=if(method=='autozyme')'0.3.1' else NA,cells=ncol(x@data.signaling),genes=nrow(x@data.signaling),K=nlevels(x@idents),nLR=nrow(x@LR$LRsig),nboot=100,seed=1,raw.use=settings$raw.use),file.path(a[3],'run_identity.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
cat('Smoke worker completed:',method,'\n')
