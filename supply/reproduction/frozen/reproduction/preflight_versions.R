# Rscript --vanilla preflight_versions.R RUNTIME_JSON EXTRACTED_CELLCHAT22_SOURCE OUTDIR
# Runtime JSON parsing uses jsonlite, already required by the supplied CellChat environment.
a<-commandArgs(TRUE);stopifnot(length(a)==3L)
rt<-jsonlite::fromJSON(a[1]);dir.create(a[3],recursive=TRUE,showWarnings=FALSE)
.libPaths(unique(c(rt$cellchat_22_library,rt$ess_library,rt$autozyme_library,rt$dependency_library,.libPaths())))
suppressPackageStartupMessages({library(CellChat);library(CellChatESS)})
stopifnot(as.character(packageVersion('CellChat'))=='2.2.0.9001',as.character(packageVersion('CellChatESS'))=='0.2.2',!isNamespaceLoaded('autozyme'))
expected_id<-if(!is.null(rt$ess_source_id))rt$ess_source_id else readLines('software/CellChatESS/inst/ESS_SOURCE_ID',warn=FALSE)
stopifnot(identical(readLines(system.file('ESS_SOURCE_ID',package='CellChatESS'),warn=FALSE),expected_id))
for(z in list(c('CellChat',rt$cellchat_22_library),c('CellChatESS',rt$ess_library)))stopifnot(normalizePath(find.package(z[1]))==normalizePath(file.path(z[2],z[1])))
# Prove the loaded original functions match the pinned upstream R source before activation.
e<-new.env(parent=asNamespace('CellChat'));sys.source(file.path(a[2],'R/modeling.R'),envir=e)
rows<-list()
for(n in c('computeCommunProb','triMean','computeExpr_LR','computeExpr_complex','computeExpr_coreceptor','computeExpr_agonist','computeExpr_antagonist','aggregateNet','computeCommunProbPathway')) {
 u<-get(n,envir=e);v<-getFromNamespace(n,'CellChat')
 ok<-identical(deparse(body(u)),deparse(body(v)))&&identical(formals(u),formals(v))
 rows[[length(rows)+1L]]<-data.frame(function_name=n,unmodified_source_match=ok)
}
tab<-do.call(rbind,rows);write.table(tab,file.path(a[3],'upstream_function_identity.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
stopifnot(all(tab$unmodified_source_match))
before<-lapply(c('computeCommunProb','triMean'),function(n)deparse(getFromNamespace(n,'CellChat')))
Sys.setenv(AUTOZYME_DISABLED='false',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
suppressPackageStartupMessages(library(autozyme));stopifnot(as.character(packageVersion('autozyme'))=='0.3.1',normalizePath(find.package('autozyme'))==normalizePath(file.path(rt$autozyme_library,'autozyme')))
autozyme::set_threads(1L);autozyme::activate('cellchat');stopifnot(identical(unname(autozyme::status()[['cellchat']]),'active'))
changed<-vapply(seq_along(before),function(i)!identical(before[[i]],deparse(getFromNamespace(c('computeCommunProb','triMean')[i],'CellChat'))),logical(1))
stopifnot(all(changed));autozyme::deactivate('cellchat')
stopifnot(all(vapply(seq_along(before),function(i)identical(before[[i]],deparse(getFromNamespace(c('computeCommunProb','triMean')[i],'CellChat'))),logical(1))))
writeLines(capture.output(sessionInfo()),file.path(a[3],'sessionInfo.txt'))
writeLines('PASS: exact upstream R functions, exact package paths/versions, and two AutoZyme replacements verified.',file.path(a[3],'PASS.txt'))
