a<-commandArgs(TRUE);stopifnot(length(a)==3L)
x<-readRDS(a[1]);y<-readRDS(a[2]);rows<-list()
cnt<-function(p,pv){v<-p;v[pv>=.05]<-0;apply(v>0,c(1,2),sum)}
x$count<-cnt(x$prob,x$pval);y$count<-cnt(y$prob,y$pval)
for (nm in c("prob","pval","count")){
 u<-x[[nm]];v<-y[[nm]];ok<-identical(dim(u),dim(v))&&identical(dimnames(u),dimnames(v))&&all(is.finite(u))&&all(is.finite(v))
 rows[[nm]]<-data.frame(field=nm,dimensions_names_finite=ok,max_abs=if(ok)max(abs(u-v))else NA,changed=if(ok)sum(u!=v)else NA)
}
write.table(do.call(rbind,rows),a[3],sep="\t",row.names=FALSE,quote=FALSE)
