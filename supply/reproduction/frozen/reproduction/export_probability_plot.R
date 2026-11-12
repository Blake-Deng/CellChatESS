a<-commandArgs(TRUE);stopifnot(length(a)==3L,!file.exists(a[3]))
x<-readRDS(a[1]);y<-readRDS(a[2]);stopifnot(identical(dim(x$prob),dim(y$prob)),identical(dimnames(x$prob),dimnames(y$prob)),all(is.finite(x$prob)),all(is.finite(y$prob)))
write.table(data.frame(original=as.vector(x$prob),ESS=as.vector(y$prob),difference=as.vector(y$prob-x$prob)),a[3],sep='\t',row.names=FALSE,quote=FALSE)
cat('Representative input only; max absolute difference:',max(abs(x$prob-y$prob)),'\n')
