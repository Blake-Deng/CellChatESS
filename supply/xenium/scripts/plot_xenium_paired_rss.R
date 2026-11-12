#!/usr/bin/env Rscript
# Plot actual 1-second RSS observations; lines connect samples without interpolation.
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 3L) stop('usage: ORIGINAL_DIR ESS_DIR OUTPUT_PREFIX')
orig_dir <- normalizePath(a[[1L]],mustWork=TRUE)
ess_dir <- normalizePath(a[[2L]],mustWork=TRUE)
prefix <- a[[3L]]
read_one <- function(d,method) {
  x <- read.delim(file.path(d,'rss_samples.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
  if (!nrow(x)) stop('No RSS samples for ',method)
  x$elapsed_min <- (as.numeric(x$sample_start_epoch_ms)-as.numeric(x$sample_start_epoch_ms[1L]))/60000
  x$RSS_GiB <- as.numeric(x$VmRSS_kB)/1024^2
  x$segment <- cumsum(c(TRUE, diff(as.numeric(x$sample_start_epoch_ms)) > 2500))
  x$method <- method
  run <- jsonlite::fromJSON(file.path(d,'run_metadata.json'),simplifyVector=FALSE)
  x$inference_start_min <- (as.numeric(run$function_call_start_epoch)*1000-as.numeric(x$sample_start_epoch_ms[1L]))/60000
  x$inference_end_min <- (as.numeric(run$function_call_end_epoch)*1000-as.numeric(x$sample_start_epoch_ms[1L]))/60000
  x
}
dat <- rbind(read_one(orig_dir,'Original CellChat 2.2.0.9001'),read_one(ess_dir,'CellChatESS 0.2.2'))
dat$method <- factor(dat$method,levels=c('Original CellChat 2.2.0.9001','CellChatESS 0.2.2'))
stage_lines <- unique(dat[,c('method','inference_start_min','inference_end_min')])
write.table(dat,file=paste0(prefix,'_rss_plot_data.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
library(ggplot2)
p <- ggplot(dat,aes(x=elapsed_min,y=RSS_GiB,colour=method,group=interaction(method,segment)))+
  geom_line(linewidth=.55,na.rm=TRUE)+geom_point(size=.65,alpha=.7,na.rm=TRUE)+
  geom_vline(data=stage_lines,aes(xintercept=inference_start_min),linetype='dotted',linewidth=.35,inherit.aes=FALSE)+
  geom_vline(data=stage_lines,aes(xintercept=inference_end_min),linetype='dashed',linewidth=.35,inherit.aes=FALSE)+
  facet_wrap(~method,ncol=1,scales='free_x')+
  scale_colour_manual(values=c('#315A8A','#B4474D'))+
  labs(x='Elapsed time from first RSS sample (min)',y='Process RSS (GiB)',
       title='Xenium full 627-LR run: sampled process memory',
       subtitle='Dots are raw /proc samples; dotted and dashed rules mark inference call start and return',
       colour=NULL)+
  theme_classic(base_size=11)+theme(legend.position='none',strip.background=element_blank(),
    strip.text=element_text(face='bold'),plot.title=element_text(face='bold'))
ggsave(paste0(prefix,'.pdf'),p,width=8,height=7,device=cairo_pdf)
ggsave(paste0(prefix,'.svg'),p,width=8,height=7,device=svglite::svglite)
ggsave(paste0(prefix,'.png'),p,width=8,height=7,dpi=300,bg='white')
