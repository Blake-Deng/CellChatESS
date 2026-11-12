a<-commandArgs(TRUE);stopifnot(length(a)==4L)
.libPaths(c(a[3],a[4],.libPaths()))
suppressPackageStartupMessages({library(CellChat);library(CellChatESS)})
future::plan('sequential')
stopifnot(as.character(packageVersion('CellChat'))=='2.2.0.9001',as.character(packageVersion('CellChatESS'))=='0.2.2')
dir.create(a[2],recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(a[2],'checks.tsv')))stop('Refusing to overwrite')
x<-readRDS(a[1]);x@idents<-droplevels(x@idents)
settings<-list(type='triMean',raw.use=TRUE,population.size=FALSE,nboot=100L,seed.use=1L,Kh=0.5,n=1)
original<-do.call(CellChat::computeCommunProb,c(list(object=x),settings))
suppressPackageStartupMessages(library(autozyme))
stopifnot(as.character(packageVersion('autozyme'))=='0.3.1')
autozyme::set_threads(1L);autozyme::activate('cellchat')
stopifnot(identical(unname(autozyme::status()[['cellchat']]),'active'))
auto<-do.call(CellChat::computeCommunProb,c(list(object=x),settings))
ess<-do.call(CellChatESS::computeCommunProbESS,c(list(object=x,threads=1L),settings))
checks<-list()
for(method in c('autozyme','ess')) {
 z<-if(method=='autozyme')auto else ess
 for(n in c('prob','pval')) {
  u<-original@net[[n]];v<-z@net[[n]]
  stopifnot(identical(dim(u),dim(v)),identical(dimnames(u),dimnames(v)),all(is.finite(u)),all(is.finite(v)))
  checks[[length(checks)+1L]]<-data.frame(method=method,output=n,max_absolute_difference=max(abs(u-v)),different_elements=sum(u!=v),elementwise_exact=all(u==v))
  if(n=='pval')stopifnot(all(u==v))else stopifnot(max(abs(u-v))<=1e-12)
 }
 count<-function(z)apply(z@net$prob>0 & z@net$pval<.05,c(1,2),sum)
 stopifnot(identical(count(original),count(z)))
 checks[[length(checks)+1L]]<-data.frame(method=method,output='count',max_absolute_difference=0,different_elements=0,elementwise_exact=TRUE)
}
write.table(do.call(rbind,checks),file.path(a[2],'checks.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
writeLines(capture.output(sessionInfo()),file.path(a[2],'sessionInfo.txt'))
writeLines(capture.output(autozyme::status()),file.path(a[2],'autozyme_status.txt'))
cat('THREE-METHOD SMOKE PASSED; numerical test only, no performance claim.\n')
