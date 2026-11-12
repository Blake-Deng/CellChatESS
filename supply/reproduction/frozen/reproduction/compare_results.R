a<-commandArgs(TRUE);stopifnot(length(a)==3L)
x<-readRDS(a[1]);y<-readRDS(a[2]);rows<-list()
# A real benchmark difference must be recorded without aborting the remaining
# methods or blocks. Validation/ablation comparisons retain strict exit status.
benchmark_block <- grepl('/results/new/(primary216|autozyme192)/[^/]+/block[0-9]+/[^/]+_vs_original\\.tsv$', a[3])
for(n in c('prob','pval','count')) {
 u<-x[[n]];v<-y[[n]]
 present<-!is.null(u)&&!is.null(v)
 dims_equal<-present&&identical(dim(u),dim(v))
 names_equal<-present&&identical(dimnames(u),dimnames(v))
 finite<-present&&all(is.finite(u))&&all(is.finite(v))
 comparable<-dims_equal&&finite
 rows[[n]]<-data.frame(output=n,max_absolute_difference=if(comparable) max(abs(u-v)) else NA_real_,
  different_elements=if(comparable) sum(u!=v) else NA_real_,elementwise_exact=if(comparable) all(u==v) else NA,
  dimensions_equal=dims_equal,dimnames_equal=names_equal,finite_values=finite)
}
write.table(do.call(rbind,rows),a[3],sep='\t',row.names=FALSE,quote=TRUE)

checks<-do.call(rbind,rows)
issue<-any(!checks$dimensions_equal | !checks$dimnames_equal | !checks$finite_values) ||
 any(checks$different_elements[checks$output %in% c('pval','count')] != 0,na.rm=TRUE) ||
 any(checks$max_absolute_difference[checks$output=='prob']>1e-12,na.rm=TRUE)
if(issue) {
 message('Comparison differs or is not comparable; inspect the saved table before making accuracy claims.')
 if(!benchmark_block) stop('Validation comparison failed')
}
