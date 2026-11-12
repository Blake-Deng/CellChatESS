# STAGE INPUT OUTDIR ESS_OR_VARIANT_LIBRARY CELLCHAT22_LIBRARY
# Wrap the same whole inference call for A0-A5. A1 uses its full-gene derived R wrapper.
a<-commandArgs(TRUE);stopifnot(length(a)==5L,a[1]%in%paste0('A',0:5))
.libPaths(c(a[4],a[5],.libPaths()))
Sys.setenv(AUTOZYME_DISABLED='true')
suppressPackageStartupMessages(library(CellChat));stopifnot(as.character(packageVersion('CellChat'))=='2.2.0.9001')
if(a[1]!='A0') {suppressPackageStartupMessages(library(CellChatESS));stopifnot(as.character(packageVersion('CellChatESS'))=='0.2.2')}
if(dir.exists(a[3]))stop('Choose a new stage output directory');dir.create(a[3],recursive=TRUE)
stopifnot(!isNamespaceLoaded('autozyme'),normalizePath(find.package('CellChat'))==normalizePath(file.path(a[5],'CellChat')))
if(a[1]!='A0')stopifnot(normalizePath(find.package('CellChatESS'))==normalizePath(file.path(a[4],'CellChatESS')))
future::plan('sequential');x<-readRDS(a[2]);x@idents<-droplevels(x@idents)
p<-list(type='triMean',raw.use=TRUE,population.size=FALSE,nboot=100L,seed.use=1L,Kh=.5,n=1)
gc();start<-proc.time()[['elapsed']]
z<-if(a[1]=='A0')do.call(CellChat::computeCommunProb,c(list(object=x),p))else do.call(CellChatESS::computeCommunProbESS,c(list(object=x,threads=1L,algorithm=if(a[1]%in%c('A1','A2'))'dense'else'ess'),p))
elapsed<-proc.time()[['elapsed']]-start
stopifnot(all(is.finite(z@net$prob)),all(is.finite(z@net$pval)))
count<-apply(z@net$prob>0 & z@net$pval<.05,c(1,2),sum)
saveRDS(list(prob=z@net$prob,pval=z@net$pval,count=count,metadata=z@options$CellChatESS),file.path(a[3],'result.rds'),compress=FALSE)
build<-if(a[1]=='A0')'upstream_original' else as.character(utils::packageDescription('CellChatESS')[['Config/CellChatESS/Build']])
if(!length(build))build<-'final_candidate'
write.table(data.frame(stage=a[1],inference_call_sec=elapsed,cells=ncol(x@data.signaling),K=nlevels(x@idents),nLR=nrow(x@LR$LRsig),CellChat='2.2.0.9001',ESS=if(a[1]=='A0')NA_character_ else '0.2.2',build=build,library=a[4],input=a[2]),file.path(a[3],'measurement.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
m<-z@options$CellChatESS
if(length(m)) write.table(as.data.frame(m),file.path(a[3],'kernel_metadata.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
writeLines(capture.output(sessionInfo()),file.path(a[3],'sessionInfo.txt'))
