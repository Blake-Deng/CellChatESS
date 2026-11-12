args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=5L)stop('usage: trace_xenium_original_bootstraps.R INPUT_RDS DIFF_TSV ORIGINAL_RDS ESS_RDS OUT_TSV')
input_path<-args[[1L]];diff_path<-args[[2L]];orig_path<-args[[3L]];ess_path<-args[[4L]];out_path<-args[[5L]]
archived_lib <- Sys.getenv("CELLCHAT_TRACE_R_LIB", "environment/historical_trace/library")
.libPaths(c(archived_lib,.libPaths()))
suppressPackageStartupMessages({library(CellChat);library(CellChatESS);library(Matrix)})
obj<-readRDS(input_path);diffs<-read.delim(diff_path,check.names=FALSE,stringsAsFactors=FALSE)
orig<-readRDS(orig_path);ess<-readRDS(ess_path)
group_levels<-dimnames(orig@net$prob)[[1L]];lr_names<-dimnames(orig@net$prob)[[3L]];K<-length(group_levels)
targets<-data.frame(lr_name=as.character(diffs$lr_name),source_group=as.character(diffs$source_group),target_group=as.character(diffs$target_group),
  source_index=vapply(diffs$source_group,function(x)match(x,group_levels),integer(1)),
  target_index=vapply(diffs$target_group,function(x)match(x,group_levels),integer(1)),stringsAsFactors=FALSE)
targets$lr_index<-vapply(targets$lr_name,function(x)match(x,lr_names),integer(1))
if(anyNA(targets$lr_index)||anyNA(targets$source_index)||anyNA(targets$target_index))stop('target labels are not present in reference result')
target_rows<-match(unique(targets$lr_name),rownames(obj@LR$LRsig))
if(anyNA(target_rows))stop('target LR rows are missing from the input')
pairs<-obj@LR$LRsig[target_rows,,drop=FALSE]
complex<-obj@DB$complex;cofactor<-obj@DB$cofactor;cf_cols<-grep('cofactor',colnames(cofactor),value=TRUE)
items<-unique(unlist(lapply(seq_len(nrow(pairs)),function(i)as.character(pairs[i,c('ligand','receptor','co_A_receptor','co_I_receptor','agonist','antagonist')])),use.names=FALSE))
items<-items[!is.na(items)&nzchar(items)]
needed<-character()
for(x in items){
  if(x%in%rownames(complex)){
    z<-unlist(complex[x,grep('^subunit',colnames(complex)),drop=FALSE],use.names=FALSE);needed<-c(needed,z[!is.na(z)&nzchar(z)])
  }else if(x%in%rownames(cofactor)){
    z<-unlist(cofactor[x,cf_cols,drop=FALSE],use.names=FALSE);needed<-c(needed,z[!is.na(z)&nzchar(z)])
  }else needed<-c(needed,x)
}
needed<-intersect(unique(needed),rownames(obj@data.signaling))
global_max<-max(obj@data.signaling);max_pos<-which(obj@data.signaling==global_max,arr.ind=TRUE)[1L,1L]
obj@data.signaling<-obj@data.signaling[unique(c(needed,rownames(obj@data.signaling)[max_pos])),,drop=FALSE]
if(max(obj@data.signaling)!=global_max)stop('gene projection changed global normalization')
obj@LR$LRsig<-pairs
fun<-CellChat::computeCommunProb
src<-deparse(fun,width.cutoff=500L)
idx<-grep('nReject <- rowSums\\(Pboot - Pnull > 0\\)',src)
if(length(idx)!=1L)stop('expected exactly one P-value comparison statement; found ',length(idx))
hook<-c(
 'nReject <- rowSums(Pboot - Pnull > 0)',
 'if (i %in% trace_lr_index) {',
 '  hit <- which(trace_lr_index == i)',
 '  for (hh in hit) {',
 '    ss <- trace_source_index[[hh]]; rr <- trace_target_index[[hh]]',
 '    lin <- ss + numCluster * (rr - 1L)',
 '    trace_rows[[length(trace_rows) + 1L]] <<- data.frame(discrepancy_id=hh, bootstrap=seq_len(nboot),',
 '      pnull_original=Pnull[[lin]], pboot_original=Pboot[lin,], delta_original=Pboot[lin,]-Pnull[[lin]],',
 '      reject_original=(Pboot[lin,]-Pnull[[lin]]) > 0, stringsAsFactors=FALSE)',
 '  }',
 '}')
src<-append(src,hook,after=idx-1L)
trace_env<-new.env(parent=asNamespace('CellChat'))
targets$trace_index<-match(targets$lr_name,rownames(pairs))
trace_env$trace_lr_index<-targets$trace_index
trace_env$trace_source_index<-targets$source_index
trace_env$trace_target_index<-targets$target_index
trace_env$trace_rows<-list()
trace_fun<-eval(parse(text=src),envir=trace_env)
cat('Running instrumented original CellChat for ',nrow(obj@LR$LRsig),' LR pairs, ',nrow(obj@data.signaling),' genes and ',nrow(targets),' target cells\n',sep='')
invisible(trace_fun(obj,LR.use=obj@LR$LRsig,type='triMean',raw.use=TRUE,population.size=FALSE,nboot=100L,seed.use=1L,Kh=0.5,n=1))
captured<-trace_env$trace_rows
if(length(captured)!=nrow(targets))stop('captured ',length(captured),' target rows; expected ',nrow(targets))
boot<-do.call(rbind,captured)
targets$discrepancy_id<-seq_len(nrow(targets))
boot<-merge(boot,targets,by='discrepancy_id',all.x=TRUE,sort=FALSE)
boot<-boot[order(boot$discrepancy_id,boot$bootstrap),]
reference_counts<-vapply(seq_len(nrow(targets)),function(j){
  si<-targets$source_index[[j]];ri<-targets$target_index[[j]];li<-targets$lr_index[[j]]
  as.integer(round(orig@net$pval[si,ri,li]*100))
},integer(1))
ess_counts<-vapply(seq_len(nrow(targets)),function(j){
  si<-targets$source_index[[j]];ri<-targets$target_index[[j]];li<-targets$lr_index[[j]]
  as.integer(round(ess@net$pval[si,ri,li]*100))
},integer(1))
boot$original_saved_pval<-orig@net$pval[cbind(boot$source_index,boot$target_index,boot$lr_index)]
boot$ess_saved_pval<-ess@net$pval[cbind(boot$source_index,boot$target_index,boot$lr_index)]
boot$original_saved_rejects<-reference_counts[boot$discrepancy_id]
boot$ess_saved_rejects<-ess_counts[boot$discrepancy_id]
summaries<-do.call(rbind,lapply(seq_len(nrow(targets)),function(j){
  z<-boot[boot$discrepancy_id==j,]
  data.frame(discrepancy_id=j,lr_name=targets$lr_name[[j]],source_group=targets$source_group[[j]],target_group=targets$target_group[[j]],
    original_reconstructed_rejects=sum(z$reject_original),original_saved_rejects=reference_counts[[j]],ess_saved_rejects=ess_counts[[j]],
    n_bootstrap_values=length(z$bootstrap),min_abs_delta=min(abs(z$delta_original)),min_delta=min(z$delta_original),max_delta=max(z$delta_original),
    stringsAsFactors=FALSE)
}))
write.table(boot,out_path,sep='\t',row.names=FALSE,quote=TRUE)
write.table(summaries,sub('\\.tsv$','_summary.tsv',out_path),sep='\t',row.names=FALSE,quote=TRUE)
writeLines(capture.output(sessionInfo()),sub('\\.tsv$','_sessionInfo.txt',out_path))
cat('Reference saved reject counts: ',paste(reference_counts,collapse=','),'\n',sep='')
cat('Reconstructed reject counts: ',paste(summaries$original_reconstructed_rejects,collapse=','),'\n',sep='')
cat('ESS saved reject counts: ',paste(ess_counts,collapse=','),'\n',sep='')
