a<-commandArgs(TRUE);stopifnot(length(a)==2L)
.libPaths(c(a[2],.libPaths()))
suppressPackageStartupMessages({library(CellChat);library(Matrix)})
stopifnot(as.character(packageVersion('CellChat'))=='2.2.0.9001',!file.exists(a[1]))
set.seed(31022)
genes<-c('L','R','L1','L2','R1','R2','CA1','CI1','AG1','AN1','ANCHOR')
data<-matrix(runif(length(genes)*120,0.05,0.9),nrow=length(genes),dimnames=list(genes,paste0('cell',1:120)))
data['ANCHOR',1]<-2
meta<-data.frame(group=factor(rep(c('A','B','C'),each=40)),row.names=colnames(data))
x<-CellChat::createCellChat(Matrix::Matrix(data,sparse=TRUE),meta=meta,group.by='group')
complex<-data.frame(subunit_1=c('L1','R1'),subunit_2=c('L2','R2'),row.names=c('CL','CR'))
cofactor<-data.frame(cofactor1=c('CA1','CI1','AG1','AN1'),cofactor2=rep('',4),row.names=c('CA','CI','AG','AN'))
lr<-data.frame(interaction_name=paste0('LR',1:7),pathway_name=paste0('P',c(1,2,3,4,4,5,5)),
 ligand=c('L','CL','L','L','L','L','L'),receptor=c('R','R','CR','R','R','R','R'),
 agonist=c('','','','','','AG',''),antagonist=c('','','','','','','AN'),
 co_A_receptor=c('','','','CA','','',''),co_I_receptor=c('','','','','CI','',''),
 annotation=rep('Secreted Signaling',7),row.names=paste0('LR',1:7))
x@DB<-list(interaction=lr,complex=complex,cofactor=cofactor,geneInfo=data.frame(Symbol=genes))
x@LR<-list(LRsig=lr);x@data.signaling<-Matrix::Matrix(data,sparse=TRUE)
saveRDS(x,a[1],compress=FALSE)
cat('Synthetic numerical fixture only: all seven expression/modifier structures, 120 cells and three groups.\n')
