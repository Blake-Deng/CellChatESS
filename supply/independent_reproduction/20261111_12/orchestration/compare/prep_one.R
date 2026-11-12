a<-commandArgs(TRUE);stopifnot(length(a)==7)
suppressPackageStartupMessages({library(Seurat);library(CellChat);library(Matrix)})
stopifnot(as.character(packageVersion('CellChat'))=='1.6.1',as.character(packageVersion('Seurat'))=='5.5.1')
# Load only the definitions of the archived original benchmark. Never execute its main().
code<-parse(a[1]);stopifnot(identical(code[[length(code)]],quote(main())))
for (e in head(as.list(code),-1)) eval(e,.GlobalEnv)
src<-a[2];target<-a[3];seed<-as.integer(a[4]);expected<-a[5];out<-a[6];eid<-a[7]
dir.create(out,recursive=TRUE,showWarnings=FALSE);stopifnot(!file.exists(file.path(out,'prep_validation.tsv')))
set.seed(seed)
seu<-load_seurat(src)
seu<-downsample_cells(seu,target,seed)
chosen<-choose_label_col(seu,'auto')
seu<-prepare_groups(seu,chosen,seed=seed)
regen<-make_cellchat(seu);rm(seu);gc(FALSE)
regfile<-file.path(out,'regenerated_cellchat.rds');saveRDS(regen,regfile);ref<-readRDS(expected)
mat<-function(u,v) {
 dm<-identical(dim(u),dim(v));names<-dm&&identical(dimnames(u),dimnames(v))
 if(!names)return(c(identical=FALSE,shape_names=FALSE,max_abs=NA,n_diff=NA))
 dd<-u-v
 c(identical=identical(u,v),shape_names=TRUE,max_abs=if(length(dd)) max(abs(dd)) else 0,n_diff=Matrix::nnzero(dd))
}
checks<-list()
for(name in c('data','data.signaling')){
 u<-slot(ref,name);v<-slot(regen,name)
 checks[[name]]<-mat(u,v)
}
exact<-function(u,v) c(identical=identical(u,v),shape_names=NA,max_abs=NA,n_diff=NA)
checks$identities<-exact(ref@idents,regen@idents)
checks$cellchat_group_metadata<-exact(ref@meta$cellchat_group,regen@meta$cellchat_group)
checks$database<-exact(ref@DB,regen@DB)
checks$evaluated_LR_table<-exact(ref@LR$LRsig,regen@LR$LRsig)
z<-data.frame(id=eid,slot=names(checks),do.call(rbind,checks),branch=if(is.na(chosen))'Seurat de novo clustering' else 'source metadata',selected_label=if(is.na(chosen))'' else chosen,source_cells=ncol(ref@data),fresh_K=nlevels(regen@idents),expected_K=nlevels(ref@idents),stringsAsFactors=FALSE)
write.table(z,file.path(out,'prep_validation.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
