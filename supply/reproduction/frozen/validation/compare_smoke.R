a<-commandArgs(TRUE);stopifnot(length(a) %in% c(1L,2L))
legacy16<-length(a)==2L && a[2]=='legacy16'
root<-a[1];paths<-list.dirs(root,recursive=FALSE,full.names=TRUE)
paths<-paths[file.exists(file.path(paths,'result.rds'))]
ref<-readRDS(file.path(root,'original','result.rds'));rows<-list()
compare<-function(u,v,method,output,exact) {
 shape<-identical(dim(u),dim(v))&&identical(dimnames(u),dimnames(v))&&identical(names(u),names(v))&&length(u)==length(v)
 finite<-is.numeric(u)&&is.numeric(v)&&all(is.finite(u))&&all(is.finite(v))
 md<-if(shape&&finite&&length(u))max(abs(u-v))else if(shape&&finite)0 else NA_real_
 equal<-shape&&if(finite)all(u==v)else identical(u,v)
 ok<-shape&&if(exact)equal else finite&&md<=1e-12
 rows[[length(rows)+1L]]<<-data.frame(method=method,output=output,elements=length(u),shape_and_names_match=shape,max_abs_difference=md,different_elements=if(shape)sum(u!=v)else NA,elementwise_equal=equal,require_exact=exact,status=if(ok)'PASS' else 'FAIL')
}
for(p in paths) {
 method<-basename(p);z<-readRDS(file.path(p,'result.rds'))
 for(n in c('prob','pval','count','weight','pathway_prob','pathways','LRs','rng')) compare(ref[[n]],z[[n]],method,n,n %in% c('pval','count','pathways','LRs','rng'))
 if(method %in% c('ess','dense')) for(n in names(z$native_reference)) compare(z$native_reference[[n]],z[[n]],paste0(method,'_downstream_API'),n,n %in% c('count','pathways','LRs'))
}
tab<-do.call(rbind,rows);write.table(tab,file.path(root,'checks.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
if(legacy16) {
 known<-tab$output=='LRs' & is.null(ref$LRs)
 tab$status[known]<-'KNOWN_VERSION_DIFFERENCE: CellChat 1.6.1 has no net$LRs field'
 write.table(tab,file.path(root,'checks.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
}
if(any(tab$status=='FAIL'))stop('Smoke comparison failed; see checks.tsv')
cat('All smoke comparisons passed:',nrow(tab),'checks\n')
