# INPUT TARGETS_TSV OUTDIR TRACE_ESS_LIBRARY CELLCHAT22_LIBRARY
# Targets columns: lr_name, source_group, target_group. All cells/groups are retained.
a<-commandArgs(TRUE);stopifnot(length(a)==5L)
options(digits=17)
.libPaths(c(a[4],a[5],.libPaths()))
suppressPackageStartupMessages({library(CellChat);library(CellChatESS);library(Matrix)})
expected_version<-Sys.getenv('CELLCHAT_REFERENCE_VERSION','2.2.0.9001')
stopifnot(expected_version %in% c('1.6.1','2.2.0.9001'),as.character(packageVersion('CellChat'))==expected_version,as.character(packageVersion('CellChatESS'))=='0.2.2')
if(dir.exists(a[3]))stop('Choose a fresh trace output directory')
dir.create(a[3],recursive=TRUE);outdir<-normalizePath(a[3])
stopifnot(!isNamespaceLoaded('autozyme'),normalizePath(find.package('CellChatESS'))==normalizePath(file.path(a[4],'CellChatESS')),
 as.character(utils::packageDescription('CellChatESS')[['Config/CellChatESS/Build']])=='TRACE')
future::plan('sequential')
x<-readRDS(a[1]);targets<-read.delim(a[2],stringsAsFactors=FALSE)
stopifnot(all(c('lr_name','source_group','target_group')%in%names(targets)))
x@idents<-droplevels(x@idents);all_cells<-ncol(x@data.signaling);all_K<-nlevels(x@idents)
lrs<-unique(targets$lr_name);rows<-match(lrs,rownames(x@LR$LRsig));stopifnot(!anyNA(rows))
x@LR$LRsig<-x@LR$LRsig[rows,,drop=FALSE]
data<-x@data.signaling;normalizer<-max(data);stopifnot(is.finite(normalizer),normalizer>0)
if(inherits(data,'dgCMatrix')) {anchor<-rownames(data)[data@i[which.max(data@x)]+1L]} else {anchor<-rownames(data)[which(data==normalizer,arr.ind=TRUE)[1,1]]}
make_index<-getFromNamespace('make_index_matrix','CellChatESS')
indices<-list(
 make_index(as.character(x@LR$LRsig$ligand),rownames(data),complex_input=x@DB$complex,keep_duplicates=TRUE),
 make_index(as.character(x@LR$LRsig$receptor),rownames(data),complex_input=x@DB$complex,keep_duplicates=TRUE))
for(field in c('co_A_receptor','co_I_receptor','agonist','antagonist')) {
 indices[[length(indices)+1L]]<-make_index(as.character(x@LR$LRsig[[field]]),rownames(data),cofactor_input=x@DB$cofactor,cofactor_cols=grepl('cofactor',colnames(x@DB$cofactor)))
}
genes<-unique(unlist(indices,use.names=FALSE));needed<-rownames(data)[genes[genes>0L]]
keep<-unique(c(intersect(needed,rownames(data)),anchor));x@data.signaling<-data[keep,,drop=FALSE]
rm(data);gc()
stopifnot(max(x@data.signaling)==normalizer,ncol(x@data.signaling)==all_cells,nlevels(x@idents)==all_K)
targets$lr1<-match(targets$lr_name,rownames(x@LR$LRsig));s<-match(targets$source_group,levels(x@idents));r<-match(targets$target_group,levels(x@idents))
stopifnot(!anyNA(s),!anyNA(r));targets$pair0<-(s-1L)+all_K*(r-1L)
stopifnot(!anyDuplicated(paste(targets$lr1,targets$pair0)))
map<-file.path(outdir,'trace_targets_0based.tsv')
write.table(data.frame(lr0=targets$lr1-1L,pair0=targets$pair0),map,sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
trace_path<-file.path(outdir,'ess_direct_trace.tsv')
writeLines('lr_index_1based\tpair_index_0based\tbootstrap\tP_obs\tP_boot\tdelta\treject\tP_obs_hex\tP_boot_hex',trace_path)
# Instrument copies of function bodies in child environments; installed namespaces remain untouched.
original_records<-list();perm_hashes<-list()
record_original<-function(lr1,observed,boot) {
 for(pair0 in targets$pair0[targets$lr1==lr1]) {
  p<-observed[pair0+1L];v<-boot[pair0+1L,]
  original_records[[length(original_records)+1L]]<<-data.frame(lr_index_1based=lr1,pair_index_0based=pair0,bootstrap=seq_along(v),P_obs=p,P_boot=v,delta=v-p,reject=v-p>0,P_obs_hex=sprintf('%a',p),P_boot_hex=sprintf('%a',v))
 }
}
record_permutation<-function(permutation,group) {
 labels<-matrix(as.integer(group)[permutation],nrow=nrow(permutation),ncol=ncol(permutation))
 perm_hashes$original<<-digest::digest(labels,algo='sha256')
}
record_labels<-function(labels)perm_hashes$ess<<-digest::digest(labels,algo='sha256')
counts<-c(nReject=0L,permutation=0L,group_boot=0L)
walk<-function(expr,kind) {
 if(!is.call(expr))return(expr)
 if(identical(expr[[1]],as.name('<-')) && is.symbol(expr[[2]])) {
  lhs<-as.character(expr[[2]])
  if(kind=='original' && lhs=='nReject') {counts['nReject']<<-counts['nReject']+1L;return(substitute({.record_original(i,Pnull,Pboot); OLD},list(OLD=expr)))}
  if(kind=='original' && lhs=='permutation') {counts['permutation']<<-counts['permutation']+1L;return(substitute({OLD; .record_permutation(permutation,group)},list(OLD=expr)))}
  if(kind=='ess' && lhs=='group_boot' && is.call(expr[[3]]) && identical(expr[[3]][[1]],as.name('matrix'))) {counts['group_boot']<<-counts['group_boot']+1L;return(substitute({OLD; .record_labels(group_boot)},list(OLD=expr)))}
 }
 # Missing formal arguments in function calls are preserved by index assignment.
 for(i in seq_along(expr)[-1]) {
  if(identical(expr[[i]],quote(expr=)))next
  expr[i]<-list(walk(expr[[i]],kind))
 }
 expr
}
orig<-CellChat::computeCommunProb;ess<-CellChatESS::computeCommunProbESS
orig_env<-new.env(parent=environment(orig));orig_env$.record_original<-record_original;orig_env$.record_permutation<-record_permutation
body(orig)<-walk(body(orig),'original');environment(orig)<-orig_env
ess_env<-new.env(parent=environment(ess));ess_env$.record_labels<-record_labels
body(ess)<-walk(body(ess),'ess');environment(ess)<-ess_env
stopifnot(all(counts==1L))
writeLines(deparse(orig),file.path(outdir,'original_diagnostic_function.R'));writeLines(deparse(ess),file.path(outdir,'ess_diagnostic_function.R'))
params<-list(type='triMean',raw.use=TRUE,population.size=FALSE,nboot=100L,seed.use=1L,Kh=.5,n=1)
original<-do.call(orig,c(list(object=x),params))
Sys.setenv(CELLCHATESS_TRACE_TARGETS=map,CELLCHATESS_TRACE_OUTPUT=trace_path)
corrected<-do.call(ess,c(list(object=x,algorithm='ess',threads=1L),params))
Sys.unsetenv(c('CELLCHATESS_TRACE_TARGETS','CELLCHATESS_TRACE_OUTPUT'))
o<-do.call(rbind,original_records);e<-read.delim(trace_path)
key<-function(z)paste(z$lr_index_1based,z$pair_index_0based,z$bootstrap,sep=':')
stopifnot(nrow(o)==nrow(targets)*100L,nrow(e)==nrow(o),!anyDuplicated(key(o)),!anyDuplicated(key(e)))
e<-e[match(key(o),key(e)),];stopifnot(identical(key(o),key(e)))
comparison<-data.frame(lr_index_1based=o$lr_index_1based,pair_index_0based=o$pair_index_0based,bootstrap=o$bootstrap,original_Pobs=o$P_obs,original_Pboot=o$P_boot,ESS_Pobs=e$P_obs,ESS_Pboot=e$P_boot,original_reject=as.logical(o$reject),ESS_reject=as.logical(e$reject),decision_equal=as.logical(o$reject)==as.logical(e$reject),original_Pobs_hex=o$P_obs_hex,original_Pboot_hex=o$P_boot_hex,ESS_Pobs_hex=e$P_obs_hex,ESS_Pboot_hex=e$P_boot_hex)
o_export<-o;for(n in c('P_obs','P_boot','delta'))o_export[[n]]<-sprintf('%.17g',o[[n]])
write.table(o_export,file.path(outdir,'original_direct_trace.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
comparison$original_delta<-o$delta;comparison$ESS_delta<-e$delta
comparison$original_delta_hex<-sprintf('%a',o$delta);comparison$ESS_delta_hex<-sprintf('%a',e$delta)
comparison_export<-comparison
for(n in c('original_Pobs','original_Pboot','ESS_Pobs','ESS_Pboot','original_delta','ESS_delta'))comparison_export[[n]]<-sprintf('%.17g',comparison[[n]])
write.table(comparison_export,file.path(outdir,'comparison_all_decisions.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
ix<-cbind(s,r,targets$lr1)
summary<-data.frame(targets,original_pvalue=original@net$pval[ix],ESS_pvalue=corrected@net$pval[ix],original_count=round(original@net$pval[ix]*100),ESS_count=round(corrected@net$pval[ix]*100))
summary$original_Pobs<-sprintf('%.17g',original@net$prob[ix])
summary$ESS_Pobs<-sprintf('%.17g',corrected@net$prob[ix])
summary$original_Pobs_hex<-sprintf('%a',original@net$prob[ix])
summary$ESS_Pobs_hex<-sprintf('%a',corrected@net$prob[ix])
summary$original_reject_count<-vapply(seq_len(nrow(targets)),function(i)sum(o$reject[o$lr_index_1based==targets$lr1[i] & o$pair_index_0based==targets$pair0[i]]),numeric(1))
summary$ESS_reject_count<-vapply(seq_len(nrow(targets)),function(i)sum(e$reject[e$lr_index_1based==targets$lr1[i] & e$pair_index_0based==targets$pair0[i]]),numeric(1))
write.table(summary,file.path(outdir,'target_summary.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
write.table(data.frame(cells=all_cells,K=all_K,target_LR=nrow(x@LR$LRsig),retained_genes=length(keep),normalization_max=normalizer,anchor_gene=anchor,comparisons=nrow(o),different_decisions=sum(!comparison$decision_equal),permutation_hash_original=perm_hashes$original,permutation_hash_ESS=perm_hashes$ess,permutations_equal=identical(perm_hashes$original,perm_hashes$ess)),file.path(outdir,'scope_and_checks.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
writeLines(capture.output(sessionInfo()),file.path(outdir,'sessionInfo.txt'))
saveRDS(list(original_prob=original@net$prob,original_pval=original@net$pval,ess_prob=corrected@net$prob,ess_pval=corrected@net$pval),file.path(outdir,'target_LR_arrays.rds'),compress=FALSE)
trace_status <- data.frame(comparisons=nrow(comparison), different_decisions=sum(!comparison$decision_equal),
  pvalue_different_targets=sum(summary$original_pvalue != summary$ESS_pvalue),
  permutations_equal=identical(perm_hashes$original,perm_hashes$ess),
  conclusion=if(all(comparison$decision_equal) && identical(perm_hashes$original,perm_hashes$ess))
    'all target decisions equal' else 'differences retained for diagnosis')
write.table(trace_status,file.path(outdir,'trace_status.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
cat('DIRECT TRACE COMPLETED:',nrow(o),'decisions; differences are retained rather than hidden.\n')
