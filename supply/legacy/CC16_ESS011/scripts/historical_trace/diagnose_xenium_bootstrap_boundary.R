args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=4L)stop('usage: diagnose_xenium_bootstrap_boundary.R INPUT_RDS DIFF_TSV RAW_OUTPUTS_RDS OUT_TSV')
input_path<-args[[1L]];diff_path<-args[[2L]];raw_path<-args[[3L]];out_path<-args[[4L]]
archived_lib <- Sys.getenv("CELLCHAT_TRACE_R_LIB", "environment/historical_trace/library");.libPaths(c(archived_lib,.libPaths()))
suppressPackageStartupMessages({library(CellChat);library(CellChatESS);library(Matrix)})
obj<-readRDS(input_path);raw<-readRDS(raw_path);diffs<-read.delim(diff_path,check.names=FALSE,stringsAsFactors=FALSE)
set.seed(1L);perm_original<-replicate(100L,sample.int(ncol(obj@data.signaling),ncol(obj@data.signaling)))
set.seed(1L);perm_ess<-replicate(100L,sample.int(ncol(obj@data.signaling),ncol(obj@data.signaling)))
permutations_identical<-identical(perm_original,perm_ess)
if(!permutations_identical)stop('permutation sequences unexpectedly differ')
group<-obj@idents; group_int<-as.integer(group); group_levels<-levels(group); K<-nlevels(group); nC<-ncol(obj@data.signaling)
orig_prob<-raw$original_prob;ess_prob<-raw$ess_prob;orig_pval<-raw$original_pval;ess_pval<-raw$ess_pval
lr_names<-raw$lr_names;group_names<-raw$groups
get_ids<-function(x,db){if(is.na(x)||!nzchar(x))return(character());if(x%in%rownames(db)){z<-unlist(db[x,grep('^subunit',colnames(db)),drop=FALSE],use.names=FALSE);return(z[!is.na(z)&nzchar(z)])};x}
complex<-obj@DB$complex; cofactor<-obj@DB$cofactor
modifier_cols<-grep('cofactor',colnames(cofactor),value=TRUE)
# Restrict expression to the ligand/receptor and any declared complexes/modifiers used by discrepant LR rows.
rows<-lapply(diffs$lr_name,function(x)match(x,rownames(obj@LR$LRsig)))
if(anyNA(unlist(rows)))stop('discrepant LR absent from prepared object')
pairs<-obj@LR$LRsig[unique(unlist(rows)),,drop=FALSE]
items<-unique(unlist(pairs[,intersect(c('ligand','receptor','co_A_receptor','co_I_receptor','agonist','antagonist'),colnames(pairs)),drop=FALSE],use.names=FALSE))
needed<-unique(unlist(lapply(items,get_ids,db=complex),use.names=FALSE))
needed<-unique(unlist(lapply(needed,function(x)if(x%in%rownames(cofactor)){z<-unlist(cofactor[x,modifier_cols,drop=FALSE],use.names=FALSE);z[!is.na(z)&nzchar(z)]}else x),use.names=FALSE))
needed<-intersect(needed,rownames(obj@data.signaling));if(!length(needed))stop('no interaction genes found')
# Only group columns referenced by the five discrepant sender/receiver cells are needed for the comparison.
needed_groups<-unique(c(match(diffs$source_group,group_names),match(diffs$target_group,group_names)))
if(anyNA(needed_groups))stop('group label mismatch')
map<-integer(K);map[needed_groups]<-seq_along(needed_groups)
data<-as.matrix(obj@data.signaling[needed,,drop=FALSE])/max(obj@data.signaling)
boot_group<-matrix(group_int[perm_original],nrow=nC,ncol=100L)
boot_group_small<-matrix(map[boot_group],nrow=nC,ncol=100L)
cpp<-CellChatESS:::group_tri_mean_boot_cpp(data,boot_group_small,length(needed_groups))$avg_boot
dim(cpp)<-c(length(needed),length(needed_groups),100L)
# R reference uses the exact CellChat triMean helper and the same sampled cell indices.
r_means<-array(NA_real_,dim=c(length(needed),length(needed_groups),100L),dimnames=list(needed,NULL,NULL))
for(b in seq_len(100L)){
  gb<-boot_group[,b]
  for(gj in seq_along(needed_groups)){
    ix<-which(gb==needed_groups[[gj]])
    if(length(ix)) for(gi in seq_along(needed)) r_means[gi,gj,b]<-CellChat::triMean(data[gi,ix])
  }
  if(b%%10L==0L)cat('bootstrap ',b,'/100\n',sep='')
}
base_rows<-list();boundary_rows<-list()
for(j in seq_len(nrow(diffs))){
  lr_name<-diffs$lr_name[[j]];lr_idx<-match(lr_name,lr_names);pair_idx<-match(lr_name,rownames(obj@LR$LRsig));pair<-obj@LR$LRsig[pair_idx,,drop=FALSE]
  s<-match(diffs$source_group[[j]],group_names);r<-match(diffs$target_group[[j]],group_names);sj<-match(s,needed_groups);rj<-match(r,needed_groups)
  orig_null<-orig_prob[s,r,lr_idx];ess_null<-ess_prob[s,r,lr_idx]
  # Confirm average-expression arithmetic on the observed (non-bootstrap) data.
  avg_R<-CellChatESS:::group_tri_mean_cpp(data,group_int,K)
  colnames(avg_R)<-group_levels;rownames(avg_R)<-needed
  lrL<-as.character(pair$ligand);lrR<-as.character(pair$receptor)
  L0<-CellChat:::computeExpr_LR(lrL,avg_R,complex)[[s]]
  R0<-CellChat:::computeExpr_LR(lrR,avg_R,complex)[[r]]
  coA0<-CellChat:::computeExpr_coreceptor(cofactor,avg_R,pair,type='A')[[r]]
  coI0<-CellChat:::computeExpr_coreceptor(cofactor,avg_R,pair,type='I')[[r]]
  R0<-R0*coA0/coI0
  # Reconstruct Pboot for the exact sender-receiver cell and compare the strict > 0 decision.
  for(b in seq_len(100L)){
    avg_Rboot<-matrix(r_means[,,b],nrow=length(needed),ncol=length(needed_groups),dimnames=list(needed,group_levels[needed_groups]))
    avg_Cboot<-matrix(cpp[,,b],nrow=length(needed),ncol=length(needed_groups),dimnames=list(needed,group_levels[needed_groups]))
    calc<-function(a){
      Lb<-CellChat:::computeExpr_LR(lrL,a,complex)[[sj]]
      Rb<-CellChat:::computeExpr_LR(lrR,a,complex)[[rj]]
      coA<-CellChat:::computeExpr_coreceptor(cofactor,a,pair,type='A')[[rj]]
      coI<-CellChat:::computeExpr_coreceptor(cofactor,a,pair,type='I')[[rj]]
      Rb<-Rb*coA/coI
      dataLR<-Lb*Rb;P1<-dataLR^1/(0.5^1+dataLR^1)
      P2<-P3<-1
      if(!is.na(pair$agonist)&&nzchar(as.character(pair$agonist)))P2<-CellChat:::computeExpr_agonist(a,pair,cofactor,index.agonist=1L,Kh=0.5,n=1)[[sj]]*CellChat:::computeExpr_agonist(a,pair,cofactor,index.agonist=1L,Kh=0.5,n=1)[[rj]]
      if(!is.na(pair$antagonist)&&nzchar(as.character(pair$antagonist)))P3<-CellChat:::computeExpr_antagonist(a,pair,cofactor,index.antagonist=1L,Kh=0.5,n=1)[[sj]]*CellChat:::computeExpr_antagonist(a,pair,cofactor,index.antagonist=1L,Kh=0.5,n=1)[[rj]]
      P1*P2*P3
    }
    pb_R<-calc(avg_Rboot);pb_C<-calc(avg_Cboot)
    # Arithmetic margin relative to each implementation's observed null.
    delta_R<-pb_R-orig_null;delta_C<-pb_C-ess_null
    boundary_rows[[length(boundary_rows)+1L]]<-data.frame(discrepancy_id=j,lr_name=lr_name,source_group=diffs$source_group[[j]],target_group=diffs$target_group[[j]],bootstrap=b,
      pnull_original=orig_null,pnull_ess=ess_null,pboot_original_reconstructed=pb_R,pboot_ess_reconstructed=pb_C,
      delta_original=delta_R,delta_ess=delta_C,reject_original=delta_R>0,reject_ess=delta_C>0,
      pboot_abs_difference=abs(pb_R-pb_C),min_abs_boundary_distance=min(abs(delta_R),abs(delta_C)),
      stringsAsFactors=FALSE)
  }
  reject_R<-vapply(Filter(function(x)x$discrepancy_id==j,boundary_rows),function(x)x$reject_original,logical(1))
  reject_C<-vapply(Filter(function(x)x$discrepancy_id==j,boundary_rows),function(x)x$reject_ess,logical(1))
  base_rows[[j]]<-data.frame(discrepancy_id=j,lr_name=lr_name,source_group=diffs$source_group[[j]],target_group=diffs$target_group[[j]],
    permutations_identical=permutations_identical,original_saved_pval=orig_pval[s,r,lr_idx],ess_saved_pval=ess_pval[s,r,lr_idx],
    original_saved_rejects=round(orig_pval[s,r,lr_idx]*100),ess_saved_rejects=round(ess_pval[s,r,lr_idx]*100),
    original_reconstructed_rejects=sum(reject_R),ess_reconstructed_rejects=sum(reject_C),
    original_null=orig_null,ess_null=ess_null,original_base_L=L0,original_base_R=R0,
    stringsAsFactors=FALSE)
}
write.table(do.call(rbind,base_rows),out_path,sep='\t',row.names=FALSE,quote=TRUE)
write.table(do.call(rbind,boundary_rows),sub('\\.tsv$','_per_bootstrap.tsv',out_path),sep='\t',row.names=FALSE,quote=TRUE)
writeLines(capture.output(sessionInfo()),sub('\\.tsv$','_sessionInfo.txt',out_path))
cat('permutations_identical=',permutations_identical,' genes=',length(needed),' groups=',length(needed_groups),' discrepancies=',nrow(diffs),'\n',sep='')
