suppressPackageStartupMessages({library(Seurat);library(CellChat)})
stopifnot(as.character(packageVersion('Seurat'))=='5.5.1',as.character(packageVersion('CellChat'))=='1.6.1')
g <- c("cell_type","celltype","cell.type","cell_type_refined","annotation","annot","seurat_clusters","cluster","CellType")
root='/data/dzf/home_migrated_20261003/cellchat_acceleration/results/runs/checkpoints'
for (grp in c('Group1','Group2')) {
 src=file.path('/data/shared/3CA_data/PRJNA543322',grp,paste0(grp,'.rds')); y=readRDS(src)
 cat(sprintf('source=%s cells=%d recognized_columns=%s annotation_selection=%s\n',src,ncol(y),paste(intersect(g,names(y@meta.data)),collapse=','),if(length(intersect(g,names(y@meta.data)))) 'labels' else 'recluster'))
 for (s in c('25000','all')) for(r in 1:3) {
  id=paste0('3CA_data__',grp,'__cells-',s,'__rep-',r,'__engine-both__ablation-full')
  x=readRDS(file.path(root,id,'prepared_cellchat.rds'))
  cat(sprintf('id=%s cells=%d K=%d idents=%s\n',id,ncol(x@data.signaling),nlevels(x@idents),paste(levels(x@idents),collapse=',')))
 }
 rm(x,y);gc()
}
