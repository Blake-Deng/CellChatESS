args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
rt <- jsonlite::fromJSON(args[1])
out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
.libPaths(unique(c(rt$cellchat_22_library, rt$ess_library, rt$autozyme_library,
                  rt$dependency_library, .libPaths())))
Sys.setenv(AUTOZYME_DISABLED='true')
suppressPackageStartupMessages({library(CellChat); library(CellChatESS)})
ip <- as.data.frame(installed.packages(noCache=TRUE), stringsAsFactors=FALSE)
write.table(ip, file.path(out,'installed_R_packages.tsv'), sep='\t', quote=TRUE,
            row.names=FALSE, na='')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines(.libPaths(),file.path(out,'library_paths.txt'))
writeLines(capture.output(Sys.getenv(c('PATH','R_LIBS','LD_LIBRARY_PATH',
  'OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'))),
  file.path(out,'runtime_environment.txt'))
packages <- c('CellChat','CellChatESS','autozyme')
recursive <- tools::package_dependencies(packages, db=installed.packages(),
                                         recursive=TRUE)
dependencies <- unique(c(packages, unlist(recursive, use.names=FALSE)))
paths <- vapply(dependencies,find.package,character(1),quiet=FALSE)
write.table(data.frame(package=dependencies,
  version=vapply(dependencies,function(p)as.character(packageVersion(p)),character(1)),
  path=paths),file.path(out,'recursive_dependency_paths.tsv'),sep='\t',
  row.names=FALSE,quote=TRUE)
write.table(data.frame(package=names(getLoadedDLLs()),
  path=vapply(getLoadedDLLs(),function(z)z[['path']],character(1))),
  file.path(out,'loaded_dll_paths.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
cat('Recorded server package inventory and recursive dependencies.\n')
