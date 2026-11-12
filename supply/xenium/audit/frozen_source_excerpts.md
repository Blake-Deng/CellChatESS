# Frozen source excerpts; no software edits

## server_install/extracted/CellChat_2.2.0.9001/CellChat-75253cd0c9e68410e6e721a6d3a0419a1d7e358f/R/modeling.R
SHA256 edf75f668a1527207c0d1a89316416becbbf7977a28f7226621d661366bd2cfe

```text
62: #'
63: computeCommunProb <- function(object, type = c("triMean", "truncatedMean","thresholdedMean", "median"), trim = 0.1, LR.use = NULL, raw.use = TRUE, population.size = FALSE,
64:                               distance.use = TRUE, interaction.range = 250, scale.distance = 0.01, k.min = 10, contact.dependent = TRUE, contact.range = NULL, contact.knn.k = NULL, contact.dependent.forced = FALSE, do.symmetric = TRUE,
65:                               nboot = 100, seed.use = 1L, Kh = 0.5, n = 1) {
66:   type <- match.arg(type)
67:   cat(type, "is used for calculating the average gene expression per cell group.", "\n")
68:   FunMean <- switch(type,
69:                     triMean = triMean,
70:                     truncatedMean = function(x) mean(x, trim = trim, na.rm = TRUE),
71:                     thresholdedMean = function(x) thresholdedMean(x, trim = trim, na.rm = TRUE),
72:                     median = function(x) median(x, na.rm = TRUE))
73: 
74:   if (raw.use) {
75:     data <- as.matrix(object@data.signaling)
76:   } else {
77:     data <- as.matrix(object@data.smooth)
78:   }
79:   if (is.null(LR.use)) {
80:     pairLR.use <- object@LR$LRsig
81:   } else {
82:     if (length(unique(LR.use$annotation)) > 1) {
83:       LR.use$annotation <- factor(LR.use$annotation, levels = c("Secreted Signaling","ECM-Receptor", "Non-protein Signaling", "Cell-Cell Contact"))
84:       LR.use <- LR.use[order(LR.use$annotation), , drop = FALSE]
85:       LR.use$annotation <- as.character(LR.use$annotation)
86:     }
87:     pairLR.use <- LR.use
88:   }
89:   complex_input <- object@DB$complex
90:   cofactor_input <- object@DB$cofactor
91:   my.sapply <- ifelse(
92:     test = future::nbrOfWorkers() == 1,
93:     yes = sapply,
94:     no = future.apply::future_sapply
95:   )
96: 
97:   ptm = Sys.time()
98: 
99:   pairLRsig <- pairLR.use
100:   group <- object@idents
101:   geneL <- as.character(pairLRsig$ligand)
102:   geneR <- as.character(pairLRsig$receptor)
103:   nLR <- nrow(pairLRsig)
104:   numCluster <- nlevels(group)
105:   if (numCluster != length(unique(group))) {
106:     stop("Please check `unique(object@idents)` and ensure that the factor levels are correct!
107:          You may need to drop unused levels using 'droplevels' function. e.g.,
108:          `meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`")
109:   }
110: 
111:   data.use <- data/max(data)
```

```text
203:   Prob <- array(0, dim = c(numCluster,numCluster,nLR))
204:   Pval <- array(0, dim = c(numCluster,numCluster,nLR))
205: 
206:   set.seed(seed.use)
207:   permutation <- replicate(nboot, sample.int(nC, size = nC))
208:   data.use.avg.boot <- my.sapply(
209:     X = 1:nboot,
210:     FUN = function(nE) {
211:       groupboot <- group[permutation[, nE]]
212:       data.use.avgB <- aggregate(t(data.use), list(groupboot), FUN = FunMean)
213:       data.use.avgB <- t(data.use.avgB[,-1])
214:       return(data.use.avgB)
215:     },
216:     simplify = FALSE
217:   )
```

```text
220:   for (i in 1:nLR) {
221:     # ligand/receptor
222:     dataLR <- Matrix::crossprod(matrix(dataLavg[i,], nrow = 1), matrix(dataRavg[i,], nrow = 1))
223:     P1 <- dataLR^n/(Kh^n + dataLR^n)
224:     P1_Pspatial <- P1*P.spatial
225:     if (sum(P1_Pspatial) == 0) {
```

```text
298:           Pboot = P1.boot*P2.boot*P3.boot*P4.boot*P.spatial
299:           return(as.vector(Pboot))
300:         }
301:       )
302:       Pboot <- matrix(unlist(Pboot), nrow=length(Pnull), ncol = nboot, byrow = FALSE)
303:       nReject <- rowSums(Pboot - Pnull > 0)
304:       p = nReject/nboot
305:       Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
306:     }
307:     setTxtProgressBar(pb = pb, value = i)
308:   }
309:   close(con = pb)
310:   Pval[Prob == 0] <- 1
311:   dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
312:   dimnames(Pval) <- dimnames(Prob)
313:   net <- list("prob" = Prob, "pval" = Pval)
314:   execution.time = Sys.time() - ptm
315:   object@options$run.time <- as.numeric(execution.time, units = "secs")
316: 
317:   object@options$parameter <- list(type.mean = type, trim = trim, raw.use = raw.use, population.size = population.size,  nboot = nboot, seed.use = seed.use, Kh = Kh, n = n,
318:                                    distance.use = distance.use, interaction.range = interaction.range, ratio = ratio, tol = tol, k.min = k.min,
319:                                    contact.dependent = contact.dependent, contact.range = contact.range, contact.knn.k = contact.knn.k, contact.dependent.forced = contact.dependent.forced
320:                                    )
321:   if (object@options$datatype != "RNA") {
322:     object@images$distance <- d.spatial
323:   }
324:   object@net <- net
325:   print(paste0('>>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [', Sys.time(),']'))
326:   return(object)
327: }
328: 
329: 
330: #' Compute the communication probability on signaling pathway level by summarizing all related ligands/receptors
331: #'
332: #' @param object CellChat object
333: #' @param net A list from object@net; If net = NULL, net = object@net
334: #' @param pairLR.use A dataframe giving the ligand-receptor interactions; If pairLR.use = NULL, pairLR.use = object@LR$LRsig
335: #' @param thresh threshold of the p-value for determining significant interaction
```

```text
403: aggregateNet <- function(object, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, remove.isolate = TRUE, thresh = 0.05, return.object = TRUE) {
404:   net <- object@net
405:   if (is.null(sources.use) & is.null(targets.use) & is.null(signaling) & is.null(pairLR.use)) {
406:     prob <- net$prob
407:     pval <- net$pval
408:     pval[prob == 0] <- 1
409:     prob[pval >= thresh] <- 0
410:     net$count <- apply(prob > 0, c(1,2), sum)
411:     net$weight <- apply(prob, c(1,2), sum)
412:     net$weight[is.na(net$weight)] <- 0
413:     net$count[is.na(net$count)] <- 0
414:   } else {
415:     df.net <- subsetCommunication(object, slot.name = "net",
416:                                   sources.use = sources.use, targets.use = targets.use,
```

## software/CellChatESS/src/cellchat_prob_cpp.cpp
SHA256 952e3c02c4c52f9190f1ef6b2345cec795f8c9d5fd3e9970fc15231375f59b87

```text
149: static inline double avg_expr_matrix(const NumericMatrix& avg, const IntegerMatrix& idx, int row, int k, bool direct_gene = true) {
150:   double sum_log = 0.0;
151:   double single_value = NA_REAL;
152:   int n = 0;
153:   for (int j = 0; j < idx.ncol(); ++j) {
154:     const int g = idx(row, j) - 1;
155:     if (g < 0) continue;
156:     const double x = avg(g, k);
157:     if (NumericVector::is_na(x)) continue;
158:     // A single gene uses the group summary directly, as in CellChat.
159:     // For multiple components, retain the original log-sum order.
160:     if (n == 0) {
161:       single_value = x;
162:     } else {
163:       if (n == 1) sum_log += std::log(single_value);
164:       sum_log += std::log(x);
165:     }
166:     ++n;
167:   }
168:   if (n == 0) return NA_REAL;
169:   if (n == 1) return direct_gene ? single_value : std::exp(std::log(single_value));
170:   return std::exp(sum_log / static_cast<double>(n));
171: }
```

```text
176:   double single_value = NA_REAL;
177:   int n = 0;
178:   const int boot_offset = G * K * b;
179:   for (int j = 0; j < idx.ncol(); ++j) {
180:     const int g = idx(row, j) - 1;
181:     if (g < 0) continue;
182:     const double x = avgBoot[g + G * k + boot_offset];
183:     if (NumericVector::is_na(x)) continue;
184:     // A single gene uses the group summary directly, as in CellChat.
185:     // For multiple components, retain the original log-sum order.
186:     if (n == 0) {
187:       single_value = x;
188:     } else {
189:       if (n == 1) sum_log += std::log(single_value);
190:       sum_log += std::log(x);
191:     }
192:     ++n;
193:   }
194:   if (n == 0) return NA_REAL;
195:   if (n == 1) return direct_gene ? single_value : std::exp(std::log(single_value));
196:   return std::exp(sum_log / static_cast<double>(n));
197: }
```

```text
472:   const int G = data.nrow();
473:   const int C = data.ncol();
474:   const int K = Rcpp::max(group);
475:   const int nboot = groupBoot.ncol();
476:   const int nLR = ligandIdx.nrow();
477:   const auto ligand_direct = direct_gene_flags(ligandIdx);
478:   const auto receptor_direct = direct_gene_flags(receptorIdx);
479:   const int KK = K * K;
480:   const R_xlen_t KK_x = static_cast<R_xlen_t>(KK);
481:   const double Kh_n = std::pow(Kh, n_power);
482:   const int requested_threads = std::max(1, threads);
483:   int threads_used = 1;
484: #ifdef _OPENMP
485:   threads_used = std::min(requested_threads, std::max(1, nboot));
486: #endif
487: 
488:   std::vector< std::vector<int> > boot_cells(nboot * K);
489:   for (int b = 0; b < nboot; ++b) {
```

```text
512:   // Layout: [bootstrap][gene][group]. The innermost group dimension is
513:   // contiguous, matching the ESS access pattern for one LR gene.
514:   std::vector<double> tri_cache(static_cast<size_t>(U) * K * nboot, NA_REAL);
515:   std::vector<unsigned char> tri_seen(static_cast<size_t>(U) * K * nboot, 0);
516: 
517:   NumericVector prob(Dimension(K, K, nLR));
518:   NumericVector pval(Dimension(K, K, nLR));
519:   std::vector<double> L(K), R(K), ago(K), ant(K);
520:   // Reuse serial scratch buffers across LR interactions.
521:   std::vector<double> Lb_serial(K), Rb_serial(K), agoB_serial(K), antB_serial(K);
522:   std::vector<int> gene_to_local(G, -1);
523: 
524:   long long active_pairs = 0;
525:   long long skipped_pairs = 0;
526:   long long boot_tri_mean_evals = 0;
527:   long long cache_hits = 0;
528:   long long streamed_lr = 0;
529:   long long max_lr_genes = 0;
```

```text
653:         Rb[k] = avg_expr_local(bootAvg, gene_to_local, receptorIdx, lr, k, K, receptor_direct[lr]);
654:         if (simple_lr) {
655:           agoB[k] = 1.0;
656:           antB[k] = 1.0;
657:         } else {
658:           const double coA = coreceptor_factor_local(bootAvg, gene_to_local, coAIdx, lr, k, K);
659:           const double coI = coreceptor_factor_local(bootAvg, gene_to_local, coIIdx, lr, k, K);
660:           Rb[k] = Rb[k] * coA / coI;
661:           agoB[k] = hasAgonist[lr] ? cofactor_factor_local(bootAvg, gene_to_local, agonistIdx, lr, k, K, true, Kh_n, n_power) : 1.0;
662:           antB[k] = hasAntagonist[lr] ? cofactor_factor_local(bootAvg, gene_to_local, antagonistIdx, lr, k, K, false, Kh_n, n_power) : 1.0;
663:         }
664:       }
665: 
666:       for (int j = 0; j < static_cast<int>(active_idx.size()); ++j) {
667:         const int idx = active_idx[j];
668:         const int s = idx % K;
669:         const int r = idx / K;
670:         const double dataLR = Lb[s] * Rb[r];
671:         double p1 = 0.0;
672:         if (dataLR > 0.0) {
673:           const double dataLR_n = std::pow(dataLR, n_power);
674:           p1 = dataLR_n / (Kh_n + dataLR_n);
675:         }
676:         const double val = p1 * (agoB[s] * agoB[r]) * (antB[s] * antB[r]);
677:         if (val - pnull[idx] > 0.0) reject_local[j] += 1;
678:       }
679:     };
680: 
681:     if (threads_used == 1) {
```

```text
737:   }
738: 
739:   const long long total_pairs = static_cast<long long>(nLR) * static_cast<long long>(KK);
740:   return List::create(
741:     Named("prob") = prob,
742:     Named("pval") = pval,
743:     Named("active_pairs") = static_cast<double>(active_pairs),
744:     Named("skipped_pairs") = static_cast<double>(skipped_pairs),
745:     Named("total_pairs") = static_cast<double>(total_pairs),
746:     Named("active_fraction") = total_pairs > 0 ? static_cast<double>(active_pairs) / static_cast<double>(total_pairs) : NA_REAL,
747:     Named("boot_tri_mean_evals") = static_cast<double>(boot_tri_mean_evals),
748:     Named("cache_hits") = static_cast<double>(cache_hits),
749:     Named("cache_genes") = static_cast<double>(U),
750:     Named("cache_slots") = static_cast<double>(static_cast<long long>(U) * static_cast<long long>(K) * static_cast<long long>(nboot)),
751:     Named("streamed_lr") = static_cast<double>(streamed_lr),
752:     Named("max_lr_genes") = static_cast<double>(max_lr_genes),
753:     Named("threads_used") = threads_used
754:   );
755: }
```

## software/CellChatESS/R/accel_cellchat.R
SHA256 ef6a07988481361a83df5fd8f3db2c693c59681877d9391e0853d03e69a4acda

```text
145:       Kh = Kh,
146:       n_power = n
147:     )
148:   }
149:   res$n_genes_kernel <- n_genes_kernel
150: 
151:   Prob <- res$prob
152:   Pval <- res$pval
153:   dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
154:   dimnames(Pval) <- dimnames(Prob)
155:   object@net <- list(prob = Prob, pval = Pval)
156:   object@options$run.time <- as.numeric(Sys.time() - ptm, units = "secs")
157:   object@options$parameter <- list(type.mean = type, trim = NULL, raw.use = raw.use,
158:                                    population.size = population.size, nboot = nboot,
159:                                    seed.use = seed.use, Kh = Kh, n = n,
160:                                    distance.use = distance.use,
161:                                    interaction.range = NULL, ratio = NULL, tol = NULL,
162:                                    k.min = NULL, contact.dependent = FALSE,
163:                                    contact.range = NULL, contact.knn.k = NULL,
164:                                    contact.dependent.forced = FALSE)
165:   object@options$CellChatESS <- list(
166:     algorithm = algorithm,
167:     active_pairs = res$active_pairs %||% NA_real_,
168:     skipped_pairs = res$skipped_pairs %||% NA_real_,
169:     total_pairs = res$total_pairs %||% NA_real_,
170:     active_fraction = res$active_fraction %||% NA_real_,
171:     n_genes_full = n_genes_full,
172:     n_genes_kernel = res$n_genes_kernel %||% n_genes_kernel,
173:     boot_tri_mean_evals = res$boot_tri_mean_evals %||% NA_real_,
174:     cache_hits = res$cache_hits %||% NA_real_,
175:     cache_genes = res$cache_genes %||% NA_real_,
176:     cache_slots = res$cache_slots %||% NA_real_,
177:     streamed_lr = res$streamed_lr %||% NA_real_,
178:     max_lr_genes = res$max_lr_genes %||% NA_real_,
179:     threads_requested = threads,
180:     threads_used = res$threads_used %||% 1L
181:   )
```

## Official filterCommunication, default min.samples

```text
908: filterCommunication <- function(object, min.cells = 10, min.samples = NULL, rare.keep = FALSE, nonFilter.keep = FALSE) {
909:   net <- object@net
910:   if (nonFilter.keep == TRUE) {
911:     cat("The non-filtered cell-cell communication is stored in `object@net$prob.nonFilter` and `object@net$pval.nonFilter`. \n")
912:     object@net$prob.nonFilter <- net$prob
913:     object@net$pval.nonFilter <- net$pval
914:   }
915:   num.interaction0 <- sum(net$prob > 0)
916:   cell.excludes <- which(as.numeric(table(object@idents)) <= min.cells)
917:   if (length(cell.excludes) > 0) {
918:     cat("The cell-cell communication related with the following cell groups are excluded due to the few number of cells: ", toString(levels(object@idents)[cell.excludes]), "!",'\t')
919:     net$prob[cell.excludes,,] <- 0
920:     net$prob[,cell.excludes,] <- 0
921:     num.interaction1 <- sum(net$prob > 0)
922:     pct.dicrease <- scales::percent((num.interaction0-num.interaction1)/num.interaction0, accuracy = .1)
923:     cat(paste0(pct.dicrease, " interactions are removed!",'\n'))
924:   } else {
925:     num.interaction1 <- num.interaction0
926:   }
927: 
928:   sample.info <- object@meta$samples
929:   sample.id <- levels(sample.info)
930:   if (is.null(min.samples)) {
931:     min.samples <- 1
932:   } else if (min.samples > length(sample.id)) {
933:     stop(paste0("There are only ", length(sample.id), " samples in the data. Please change the value of `min.samples`! "))
934:   }
935:   if (length(sample.id) >= 2 & min.samples >= 2) {
936:     if (object@options$parameter$raw.use == TRUE) {
937:       data <- as.matrix(object@data.signaling)
938:     } else {
939:       if ("data.smooth" %in% methods::slotNames(object) == FALSE) {
940:         stop("`object@data.smooth` is missing. Please update the CellChat object via `updateCellChat`! \n")
```
