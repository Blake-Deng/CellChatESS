computeCommunProbESS <- function(object, type = "triMean", raw.use = TRUE, population.size = FALSE,
                                     nboot = 100, seed.use = 1L, Kh = 0.5, n = 1,
                                     distance.use = NULL,
                                     threads = 1L,
                                     algorithm = c("ess", "dense")) {
  algorithm <- match.arg(algorithm)
  if (!is.numeric(threads) || length(threads) != 1L || !is.finite(threads) ||
      threads < 1 || threads != floor(threads) || threads > .Machine$integer.max) {
    stop("threads must be one positive integer.", call. = FALSE)
  }
  threads <- as.integer(threads)
  if (length(threads) != 1L || is.na(threads) || threads < 1L) {
    stop("threads must be one positive integer.", call. = FALSE)
  }
  check_supported_object(object, type, population.size, distance.use)
  if (!is.logical(raw.use) || length(raw.use) != 1L || is.na(raw.use)) {
    stop("raw.use must be TRUE or FALSE.", call. = FALSE)
  }

  ptm <- Sys.time()
  data <- if (raw.use) {
    object@data.signaling
  } else {
    if (!"data.smooth" %in% methods::slotNames(object)) {
      stop("object@data.smooth is missing. Run CellChat::projectData() first or use raw.use = TRUE.", call. = FALSE)
    }
    object@data.smooth
  }

  if (!length(data)) stop("No signaling data found in object@data.signaling.", call. = FALSE)
  validate_inference_input(object, data, raw.use, nboot, seed.use, Kh, n)
  validate_nonnegative_finite_expression(data)
  max_data <- max(data)
  if (!is.finite(max_data) || max_data <= 0) {
    stop("Signaling expression matrix has no positive finite values.", call. = FALSE)
  }

  pairLRsig <- object@LR$LRsig
  if (is.null(pairLRsig) || !nrow(pairLRsig)) {
    stop("object@LR$LRsig is empty. Run identifyOverExpressedInteractions() first.", call. = FALSE)
  }

  complex_input <- object@DB$complex
  cofactor_input <- object@DB$cofactor
  if (is.null(complex_input) || is.null(cofactor_input)) {
    stop("object@DB must include CellChat complex and cofactor tables.", call. = FALSE)
  }

  group <- object@idents
  if (!is.factor(group)) {
    stop("object@idents must be a factor for single-dataset CellChat objects.", call. = FALSE)
  }
  K <- nlevels(group)
  if (K != length(unique(group))) {
    stop("Unused factor levels found in object@idents. Drop unused levels before running.", call. = FALSE)
  }

  data.use <- data / max_data
  nC <- ncol(data.use)
  group_int <- as.integer(group)

  data_genes <- rownames(data.use)
  cofactor_cols <- grepl("cofactor", colnames(cofactor_input))
  ligand_idx <- make_index_matrix(as.character(pairLRsig$ligand), data_genes,
                                  complex_input = complex_input, keep_duplicates = TRUE)
  receptor_idx <- make_index_matrix(as.character(pairLRsig$receptor), data_genes,
                                    complex_input = complex_input, keep_duplicates = TRUE)
  coA_idx <- make_index_matrix(as.character(pairLRsig$co_A_receptor), data_genes,
                               cofactor_input = cofactor_input, cofactor_cols = cofactor_cols)
  coI_idx <- make_index_matrix(as.character(pairLRsig$co_I_receptor), data_genes,
                               cofactor_input = cofactor_input, cofactor_cols = cofactor_cols)
  agonist_names <- as.character(pairLRsig$agonist)
  antagonist_names <- as.character(pairLRsig$antagonist)
  hasAgonist <- !is.na(agonist_names) & nzchar(agonist_names)
  hasAntagonist <- !is.na(antagonist_names) & nzchar(antagonist_names)
  agonist_idx <- make_index_matrix(agonist_names, data_genes,
                                   cofactor_input = cofactor_input, cofactor_cols = cofactor_cols)
  antagonist_idx <- make_index_matrix(antagonist_names, data_genes,
                                      cofactor_input = cofactor_input, cofactor_cols = cofactor_cols)

  agonist_idx[!hasAgonist, ] <- 0L
  antagonist_idx[!hasAntagonist, ] <- 0L
  n_genes_full <- nrow(data.use)
  kernel_input <- list(data = as.matrix(data.use), indices = list(
    ligand=ligand_idx, receptor=receptor_idx, coA=coA_idx, coI=coI_idx,
    agonist=agonist_idx, antagonist=antagonist_idx))
  data.use <- kernel_input$data
  ligand_idx <- kernel_input$indices$ligand
  receptor_idx <- kernel_input$indices$receptor
  coA_idx <- kernel_input$indices$coA
  coI_idx <- kernel_input$indices$coI
  agonist_idx <- kernel_input$indices$agonist
  antagonist_idx <- kernel_input$indices$antagonist

  n_genes_kernel <- nrow(data.use)

  set.seed(seed.use)
  group_boot <- replicate(nboot, group_int[sample.int(nC, size = nC)])
  group_boot <- matrix(group_boot, nrow = nC, ncol = nboot)

  data.use.avg <- group_tri_mean_cpp(data.use, group_int, K)
  colnames(data.use.avg) <- levels(group)
  rownames(data.use.avg) <- rownames(data.use)

  if (identical(algorithm, "ess")) {
    res <- cellchat_prob_ess_cpp(
      data = data.use,
      group = group_int,
      groupBoot = group_boot,
      avg = data.use.avg,
      ligandIdx = ligand_idx,
      receptorIdx = receptor_idx,
      coAIdx = coA_idx,
      coIIdx = coI_idx,
      agonistIdx = agonist_idx,
      antagonistIdx = antagonist_idx,
      hasAgonist = hasAgonist,
      hasAntagonist = hasAntagonist,
      Kh = Kh,
      n_power = n,
      threads = threads
    )
  } else {
    avg_boot <- group_tri_mean_boot_cpp(data.use, group_boot, K)$avg_boot
    res <- cellchat_prob_from_avg_cpp(
      avg = data.use.avg,
      avgBoot = avg_boot,
      ligandIdx = ligand_idx,
      receptorIdx = receptor_idx,
      coAIdx = coA_idx,
      coIIdx = coI_idx,
      agonistIdx = agonist_idx,
      antagonistIdx = antagonist_idx,
      hasAgonist = hasAgonist,
      hasAntagonist = hasAntagonist,
      Kh = Kh,
      n_power = n
    )
  }
  res$n_genes_kernel <- n_genes_kernel

  Prob <- res$prob
  Pval <- res$pval
  dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
  dimnames(Pval) <- dimnames(Prob)
  object@net <- list(prob = Prob, pval = Pval)
  object@options$run.time <- as.numeric(Sys.time() - ptm, units = "secs")
  object@options$parameter <- list(type.mean = type, trim = NULL, raw.use = raw.use,
                                   population.size = population.size, nboot = nboot,
                                   seed.use = seed.use, Kh = Kh, n = n,
                                   distance.use = distance.use,
                                   interaction.range = NULL, ratio = NULL, tol = NULL,
                                   k.min = NULL, contact.dependent = FALSE,
                                   contact.range = NULL, contact.knn.k = NULL,
                                   contact.dependent.forced = FALSE)
  object@options$CellChatESS <- list(
    algorithm = algorithm,
    active_pairs = res$active_pairs %||% NA_real_,
    skipped_pairs = res$skipped_pairs %||% NA_real_,
    total_pairs = res$total_pairs %||% NA_real_,
    active_fraction = res$active_fraction %||% NA_real_,
    n_genes_full = n_genes_full,
    n_genes_kernel = res$n_genes_kernel %||% n_genes_kernel,
    boot_tri_mean_evals = res$boot_tri_mean_evals %||% NA_real_,
    cache_hits = res$cache_hits %||% NA_real_,
    cache_genes = res$cache_genes %||% NA_real_,
    cache_slots = res$cache_slots %||% NA_real_,
    streamed_lr = res$streamed_lr %||% NA_real_,
    max_lr_genes = res$max_lr_genes %||% NA_real_,
    threads_requested = threads,
    threads_used = res$threads_used %||% 1L
  )
  object
}

computeCommunProbPathwayESS <- function(object = NULL, net = NULL, pairLR.use = NULL, thresh = 0.05) {
  if (is.null(net)) net <- object@net
  if (is.null(pairLR.use)) pairLR.use <- object@LR$LRsig
  pathways <- unique(pairLR.use$pathway_name)
  pathway_id <- as.integer(factor(pairLR.use$pathway_name, levels = pathways))
  res <- pathway_sum_cpp(net$prob, net$pval, pathway_id, length(pathways), thresh)
  prob.pathways <- res$prob_pathway
  dimnames(prob.pathways) <- list(dimnames(net$prob)[[1]], dimnames(net$prob)[[2]], pathways)

  LR <- dimnames(net$prob)[[3]]
  LR.sig <- LR[res$lr_sum != 0]
  pathways.sig <- pathways[res$pathway_sum != 0]
  prob.pathways.sig <- prob.pathways[, , pathways.sig, drop = FALSE]
  if (length(pathways.sig)) {
    idx <- order(apply(prob.pathways.sig, 3, sum), decreasing = TRUE)
    pathways.sig <- pathways.sig[idx]
    # Match CellChat 2.2.0.9001, including its single-pathway dimension drop.
    prob.pathways.sig <- prob.pathways.sig[, , idx]
  }

  if (is.null(object)) {
    list(pathways = pathways.sig, prob = prob.pathways.sig)
  } else {
    object@net$LRs <- LR.sig
    object@netP$pathways <- pathways.sig
    object@netP$prob <- prob.pathways.sig
    object
  }
}

aggregateNetESS <- function(object, thresh = 0.05) {
  res <- aggregate_net_cpp(object@net$prob, object@net$pval, thresh)
  dimnames(res$count) <- dimnames(object@net$prob)[1:2]
  dimnames(res$weight) <- dimnames(object@net$prob)[1:2]
  object@net$count <- res$count
  object@net$weight <- res$weight
  object
}

computeAveExprESS <- function(object, features = NULL, group.by = NULL,
                                  type = "triMean", slot.name = c("data.signaling", "data"),
                                  data.use = NULL) {
  if (!identical(type, "triMean")) {
    stop("CellChatESS currently supports only type = 'triMean'.", call. = FALSE)
  }
  slot.name <- match.arg(slot.name)
  if (is.null(data.use)) data.use <- methods::slot(object, slot.name)
  features.use <- if (is.null(features)) row.names(data.use) else intersect(features, row.names(data.use))
  data.use <- as.matrix(data.use[features.use, , drop = FALSE])
  labels <- if (is.null(group.by)) object@idents else object@meta[[group.by]]
  if (!is.factor(labels)) labels <- factor(labels)
  avg <- group_tri_mean_cpp(data.use, as.integer(labels), nlevels(labels))
  rownames(avg) <- features.use
  colnames(avg) <- levels(labels)
  avg
}
