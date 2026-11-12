a<-commandArgs(TRUE);stopifnot(length(a)==4L)
.libPaths(c(a[2],a[3],a[4],.libPaths()))
suppressPackageStartupMessages({library(CellChat);library(CellChatESS);library(Matrix)})
stopifnot(as.character(packageVersion('CellChat'))=='2.2.0.9001')
dir.create(a[1],recursive=TRUE,showWarnings=FALSE)
rows<-list()
record<-function(name,u,v,exact=FALSE) {
 shape<-identical(dim(u),dim(v))&&identical(dimnames(u),dimnames(v))&&length(u)==length(v)
 d<-if(shape&&length(u))max(abs(u-v))else if(shape)0 else Inf
 ok<-shape&&is.finite(d)&&if(exact)all(u==v)else d<=1e-12
 rows[[length(rows)+1L]]<<-data.frame(check=name,elements=length(u),max_abs_difference=d,elementwise_equal=shape&&all(u==v),shape_match=shape,status=if(ok)'PASS' else 'FAIL')
}
reject<-function(name,expr) {
 e<-tryCatch({force(expr);NULL},error=function(e)e)
 rows[[length(rows)+1L]]<<-data.frame(check=name,elements=0,max_abs_difference=0,elementwise_equal=NA,shape_match=NA,status=if(inherits(e,'error'))'PASS' else 'FAIL')
}
base<-readRDS('validation/local_smoke/data/synthetic_structures_CC22.rds')
settings<-list(type='triMean',raw.use=TRUE,population.size=FALSE,nboot=100L,seed.use=1L,Kh=.5,n=1)
compare_inference<-function(x,label) {
 y<-do.call(CellChat::computeCommunProb,c(list(object=x),settings))
 z<-do.call(CellChatESS::computeCommunProbESS,c(list(object=x,threads=1L),settings))
 d<-do.call(CellChatESS::computeCommunProbESS,c(list(object=x,algorithm='dense'),settings))
 record(paste0(label,'_prob'),y@net$prob,z@net$prob)
 record(paste0(label,'_pval'),y@net$pval,z@net$pval,TRUE)
 record(paste0(label,'_dense_prob'),y@net$prob,d@net$prob)
 record(paste0(label,'_dense_pval'),y@net$pval,d@net$pval,TRUE)
 invisible(z)
}
x<-base
x@DB$complex<-rbind(x@DB$complex,data.frame(subunit_1=c('L1','L2'),subunit_2=c('',''),row.names=c('SINGLE_COMPLEX','L')))
x@LR$LRsig<-x@LR$LRsig[c(1,1),];rownames(x@LR$LRsig)<-c('single_complex','name_collision')
x@LR$LRsig$ligand<-c('SINGLE_COMPLEX','L')
compare_inference(x,'single_complex_and_gene_precedence')
bad<-x;bad@data.signaling<-bad@data.signaling[rownames(bad@data.signaling)!='L1',]
reject('missing_required_subunit_ESS_rejected',CellChatESS::computeCommunProbESS(bad))
reject('missing_required_subunit_reference_rejected',CellChat::computeCommunProb(bad))
x<-base;x@LR$LRsig<-x@LR$LRsig[1,,drop=FALSE];x@LR$LRsig$agonist<-'AG1'
# AG1 is a measured gene, but is not a cofactor row; reference treats it as neutral.
compare_inference(x,'unknown_modifier_not_gene_fallback')
x<-base;x@data.smooth<-as.matrix(x@data.signaling)*.8
settings$raw.use<-FALSE;compare_inference(x,'raw_use_false');settings$raw.use<-TRUE
x<-base
z<-CellChatESS::computeCommunProbESS(x);xd<-x;xd@data.signaling<-as.matrix(x@data.signaling)
zd<-CellChatESS::computeCommunProbESS(xd)
record('dense_sparse_probability',z@net$prob,zd@net$prob,TRUE)
record('dense_sparse_pval',z@net$pval,zd@net$pval,TRUE)
# Exercise actual public downstream APIs at p=.05 and one/zero retained pathways.
x<-base;x@LR$LRsig<-x@LR$LRsig[1:2,];x@LR$LRsig$pathway_name<-c('ONE','OTHER')
x@net<-list(prob=array(.25,c(3,3,2),dimnames=list(levels(x@idents),levels(x@idents),rownames(x@LR$LRsig))),pval=array(.05,c(3,3,2)))
x@net$pval[,,2]<-1
u<-CellChat::computeCommunProbPathway(x);v<-CellChatESS::computeCommunProbPathwayESS(x)
record('one_pathway_shape_and_values',u@netP$prob,v@netP$prob,TRUE)
u<-CellChat::aggregateNet(x);v<-CellChatESS::aggregateNetESS(x)
record('aggregate_excludes_p_equal_threshold',u@net$count,v@net$count,TRUE)
x@net$pval[]<-1
u<-CellChat::computeCommunProbPathway(x);v<-CellChatESS::computeCommunProbPathwayESS(x)
record('zero_pathway_shape_and_values',u@netP$prob,v@netP$prob,TRUE)
# Check both triMean implementations across the native sorting/multiselect boundary.
f<-getFromNamespace('group_tri_mean_cpp','CellChatESS')
set.seed(72)
for(nc in c(1,2,3,4,5,10,17,4095,4096,4097)) {
 m<-matrix(runif(4*nc),nrow=4);m[1,]<-0;m[2,seq_len(max(1L,floor(nc/2)))]<-0
 u<-vapply(seq_len(nrow(m)),function(i)CellChat::triMean(m[i,]),numeric(1))
 v<-as.numeric(f(m,rep(1L,nc),1L));record(paste0('triMean_group_size_',nc),u,v)
}
x<-base;x@idents<-x@idents[-1];reject('cell_count_mismatch',CellChatESS::computeCommunProbESS(x))
x<-base;names(x@idents)<-rev(colnames(x@data.signaling));reject('cell_order_mismatch',CellChatESS::computeCommunProbESS(x))
x<-base;reject('nboot_zero',CellChatESS::computeCommunProbESS(x,nboot=0))
reject('fractional_threads',CellChatESS::computeCommunProbESS(x,threads=1.5))
reject('nonpositive_hill_exponent',CellChatESS::computeCommunProbESS(x,n=0))
reject('unsupported_mean',CellChatESS::computeCommunProbESS(x,type='median'))
reject('unsupported_population_weight',CellChatESS::computeCommunProbESS(x,population.size=TRUE))
x@options$datatype<-'spatial';reject('unsupported_spatial',CellChatESS::computeCommunProbESS(x))
tab<-do.call(rbind,rows);write.table(tab,file.path(a[1],'checks.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
writeLines(capture.output(sessionInfo()),file.path(a[1],'sessionInfo.txt'))
if(any(tab$status!='PASS'))stop('Alignment checks failed; inspect checks.tsv')
cat('All targeted CellChat 2.2.0.9001 alignment checks passed:',nrow(tab),'\n')
