a <- commandArgs(TRUE); stopifnot(length(a)==2)
j <- read.delim(a[1],colClasses='character',check.names=FALSE)
out <- a[2]
if (file.exists(out)) stop('existing audit output')
for (i in seq_len(nrow(j))) {
 z <- j[i,]
 if(!file.exists(z$archive_file)) next
 if(identical(normalizePath(z$fresh_file),normalizePath(z$archive_file))) stop("same fresh/reference path")
 x <- readRDS(z$archive_file); y <- readRDS(z$fresh_file)
 for (field in c("prob","pval","count")) {
   u <- x[[field]]; v <- y[[field]]
   dims <- !is.null(u) && !is.null(v) && identical(dim(u),dim(v))
   labels <- dims && identical(dimnames(u),dimnames(v))
   finite <- dims && all(is.finite(u)) && all(is.finite(v))
   d <- if (finite) max(abs(u-v)) else NA_real_
   count <- if (finite) sum(u!=v) else NA_real_
   row <- data.frame(queue=z$queue,experiment_id=z$experiment_id,block=z$block,method=z$method,
      field=field,dimensions_match=dims,names_and_order_match=labels,finite=finite,
      elements=length(v),different_elements=count,max_abs=if(is.finite(d)) sprintf('%.17g',d) else NA_character_,
      value_identical=isTRUE(finite&&labels&&count==0))
   header <- !file.exists(out)
   write.table(row,out,sep='\t',quote=FALSE,row.names=FALSE,col.names=header,append=!header)
 }
 rm(x,y,u,v); if(i%%25==0)gc(FALSE)
}
