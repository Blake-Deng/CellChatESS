# MANIFEST OUT_TSV FINAL_ESS_LIBRARY CELLCHAT22_LIBRARY
# Direct strict positive L/Radj support; f_LR and f_union use different denominators.
a<-commandArgs(TRUE);stopifnot(length(a)==4L,!file.exists(a[2]))
.libPaths(c(a[3],a[4],.libPaths()))
suppressPackageStartupMessages({library(CellChat);library(CellChatESS);library(Matrix)})
stopifnot(as.character(packageVersion('CellChat'))=='2.2.0.9001',as.character(packageVersion('CellChatESS'))=='0.2.2')
m<-read.delim(a[1],stringsAsFactors=FALSE);stopifnot(all(c('experiment_id','dataset_id','input_file','input_sha256')%in%names(m)))
tri<-getFromNamespace('group_tri_mean_cpp','CellChatESS');lrfun<-getFromNamespace('computeExpr_LR','CellChat');cofun<-getFromNamespace('computeExpr_coreceptor','CellChat')
records<-list()
for(i in seq_len(nrow(m))) {
 stopifnot(identical(digest::digest(file=m$input_file[i],algo='sha256'),m$input_sha256[i]))
 x<-readRDS(m$input_file[i]);x@idents<-droplevels(x@idents);d<-x@data.signaling
 getFromNamespace('validate_nonnegative_finite_expression','CellChatESS')(d)
 maximum<-max(d);stopifnot(maximum>0,is.finite(maximum))
 g<-x@idents;K<-nlevels(g);pairs<-x@LR$LRsig;nLR<-nrow(pairs)
 avg<-tri(as.matrix(d/maximum),as.integer(g),K);rownames(avg)<-rownames(d);colnames(avg)<-levels(g)
 L<-lrfun(as.character(pairs$ligand),avg,x@DB$complex)
 R<-lrfun(as.character(pairs$receptor),avg,x@DB$complex)
 R<-R*cofun(x@DB$cofactor,avg,pairs,type='A')/cofun(x@DB$cofactor,avg,pairs,type='I')
 stopifnot(all(is.finite(L)),all(is.finite(R)))
 y<-CellChatESS::computeCommunProbESS(x,nboot=1L,seed.use=1L,threads=1L)
 union<-matrix(FALSE,K,K);count<-0;positive_zero<-0;nonzero_without<-0
 for(j in seq_len(nLR)) {
  support<-outer(L[j,]>0,R[j,]>0,'&');nonzero<-y@net$prob[,,j]!=0
  union<-union|support;count<-count+sum(support)
  positive_zero<-positive_zero+sum(support & !nonzero);nonzero_without<-nonzero_without+sum(nonzero & !support)
 }
 nnz<-if(inherits(d,'Matrix'))Matrix::nnzero(d)else sum(d!=0)
 sz<-as.numeric(table(g));r<-m[i,,drop=FALSE]
 r$cells<-ncol(d);r$K<-K;r$nLR<-nLR;r$signaling_genes<-nrow(d)
 r$matrix_nonzero<-as.double(nnz);r$signaling_gene_cell_elements<-as.double(nrow(d))*ncol(d)
 r$matrix_sparsity<-1-nnz/(as.double(nrow(d))*ncol(d))
 r$min_group_cells<-min(sz);r$median_group_cells<-median(sz);r$max_group_cells<-max(sz)
 r$positive_L_Radj_triples<-count;r$positive_group_pair_union<-sum(union)
 r$f_LR<-count/(as.double(nLR)*K*K);r$f_union<-sum(union)/(as.double(K)*K)
 r$positive_support_but_zero_probability<-positive_zero;r$nonzero_probability_without_support<-nonzero_without
 r$input_R_object_bytes<-as.numeric(object.size(d));r$ESS_nonzero_active_fraction<-y@options$CellChatESS$active_fraction
 r$CellChat<-'2.2.0.9001';r$CellChatESS<-'0.2.2'
 records[[i]]<-r
 rm(x,d,avg,L,R,y);gc(FALSE)
}
write.table(do.call(rbind,records),a[2],sep='\t',row.names=FALSE,quote=TRUE)
writeLines(capture.output(sessionInfo()),paste0(a[2],'.sessionInfo.txt'))
