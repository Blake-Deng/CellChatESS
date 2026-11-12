# INPUT OUTDIR ESS_LIB CELLCHAT_LIB CC_VERSION THREADS
a <- commandArgs(TRUE);stopifnot(length(a)==6L)
.libPaths(c(a[3],a[4],.libPaths()))
Sys.setenv(AUTOZYME_DISABLED='true')
suppressPackageStartupMessages({library(CellChat);library(CellChatESS)})
stopifnot(as.character(packageVersion('CellChat'))==a[5],
          as.character(packageVersion('CellChatESS'))=='0.2.2',
          normalizePath(find.package('CellChatESS'))==normalizePath(file.path(a[3],'CellChatESS')),
          !isNamespaceLoaded('autozyme'))
if(dir.exists(a[2]))stop('Fresh output directory required')
dir.create(a[2],recursive=TRUE)
x<-readRDS(a[1]);x@idents<-droplevels(x@idents)
future::plan('sequential');threads<-as.integer(a[6])
z<-CellChatESS::computeCommunProbESS(x,type='triMean',raw.use=TRUE,
  population.size=FALSE,nboot=100L,seed.use=1L,Kh=.5,n=1,threads=threads)
stopifnot(z@options$CellChatESS$threads_used==threads)
count<-apply(z@net$prob>0 & z@net$pval<.05,c(1,2),sum)
saveRDS(list(prob=z@net$prob,pval=z@net$pval,count=count),file.path(a[2],'result.rds'))
write.table(data.frame(CellChat=a[5],CellChatESS='0.2.2',requested=threads,
  threads_used=z@options$CellChatESS$threads_used,
  dll=getLoadedDLLs()[['CellChatESS']][['path']]),
  file.path(a[2],'threads.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
writeLines(capture.output(sessionInfo()),file.path(a[2],'sessionInfo.txt'))
