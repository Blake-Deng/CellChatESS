function (object, type = c("triMean", "truncatedMean", "thresholdedMean", "median"), trim = 0.1, LR.use = NULL, raw.use = TRUE, population.size = FALSE, distance.use = TRUE, interaction.range = 250, scale.distance = 0.01, k.min = 10, contact.dependent = TRUE, contact.range = NULL, contact.knn.k = NULL, contact.dependent.forced = FALSE, do.symmetric = TRUE, nboot = 100, seed.use = 1L, Kh = 0.5, n = 1) 
{
    type <- match.arg(type)
    cat(type, "is used for calculating the average gene expression per cell group.", "\n")
    FunMean <- switch(type, triMean = triMean, truncatedMean = function(x) mean(x, trim = trim, na.rm = TRUE), thresholdedMean = function(x) thresholdedMean(x, trim = trim, na.rm = TRUE), median = function(x) median(x, na.rm = TRUE))
    if (raw.use) {
        data <- as.matrix(object@data.signaling)
    }
    else {
        data <- as.matrix(object@data.smooth)
    }
    if (is.null(LR.use)) {
        pairLR.use <- object@LR$LRsig
    }
    else {
        if (length(unique(LR.use$annotation)) > 1) {
            LR.use$annotation <- factor(LR.use$annotation, levels = c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling", "Cell-Cell Contact"))
            LR.use <- LR.use[order(LR.use$annotation), , drop = FALSE]
            LR.use$annotation <- as.character(LR.use$annotation)
        }
        pairLR.use <- LR.use
    }
    complex_input <- object@DB$complex
    cofactor_input <- object@DB$cofactor
    my.sapply <- ifelse(test = future::nbrOfWorkers() == 1, yes = sapply, no = future.apply::future_sapply)
    ptm = Sys.time()
    pairLRsig <- pairLR.use
    group <- object@idents
    geneL <- as.character(pairLRsig$ligand)
    geneR <- as.character(pairLRsig$receptor)
    nLR <- nrow(pairLRsig)
    numCluster <- nlevels(group)
    if (numCluster != length(unique(group))) {
        stop("Please check `unique(object@idents)` and ensure that the factor levels are correct!\n         You may need to drop unused levels using 'droplevels' function. e.g.,\n         `meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`")
    }
    data.use <- data/max(data)
    nC <- ncol(data.use)
    data.use.avg <- aggregate(t(data.use), list(group), FUN = FunMean)
    data.use.avg <- t(data.use.avg[, -1])
    colnames(data.use.avg) <- levels(group)
    dataLavg <- computeExpr_LR(geneL, data.use.avg, complex_input)
    dataRavg <- computeExpr_LR(geneR, data.use.avg, complex_input)
    dataRavg.co.A.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avg, pairLRsig, type = "A")
    dataRavg.co.I.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avg, pairLRsig, type = "I")
    dataRavg <- dataRavg * dataRavg.co.A.receptor/dataRavg.co.I.receptor
    dataLavg2 <- t(replicate(nrow(dataLavg), as.numeric(table(group))/nC))
    dataRavg2 <- dataLavg2
    index.agonist <- which(!is.na(pairLRsig$agonist) & pairLRsig$agonist != "")
    index.antagonist <- which(!is.na(pairLRsig$antagonist) & pairLRsig$antagonist != "")
    if (object@options$datatype != "RNA") {
        data.spatial <- object@images$coordinates
        if ("spatial.factors" %in% names(object@images)) {
            ratio <- object@images$spatial.factors$ratio
            tol <- object@images$spatial.factors$tol
        }
        else {
            stop("`object@images$spatial.factors` is missing. Please update the object via `updateCellChat`! \n")
        }
        meta.t = data.frame(group = group, samples = object@meta$samples, row.names = rownames(object@meta))
        res <- computeRegionDistance(coordinates = data.spatial, meta = meta.t, interaction.range = interaction.range, ratio = ratio, tol = tol, k.min = k.min, contact.dependent = contact.dependent, contact.range = contact.range, contact.knn.k = contact.knn.k)
        d.spatial <- res$d.spatial
        adj.contact <- res$adj.contact
        if (distance.use) {
            print(paste0(">>> Run CellChat on spatial transcriptomics data using distances as constraints of the computed communication probability <<< [", Sys.time(), "]"))
            d.spatial <- d.spatial * scale.distance
            diag(d.spatial) <- NaN
            d.min <- min(d.spatial, na.rm = TRUE)
            if (d.min < 1) {
                cat("The suggested minimum value of scaled distances is in [1,2], and the calculated value here is ", d.min, "\n")
                stop("Please increase the value of `scale.distance` and use a value that is slighly smaller than ", format(1/d.min, digits = 2), "\n")
            }
            P.spatial <- 1/d.spatial
            P.spatial[is.na(d.spatial)] <- 0
            diag(P.spatial) <- max(P.spatial)
            d.spatial <- d.spatial/scale.distance
        }
        else {
            print(paste0(">>> Run CellChat on spatial transcriptomics data without distance values as constraints of the computed communication probability <<< [", Sys.time(), "]"))
            P.spatial <- matrix(1, nrow = numCluster, ncol = numCluster)
            P.spatial[is.na(d.spatial)] <- 0
        }
    }
    else {
        print(paste0(">>> Run CellChat on sc/snRNA-seq data <<< [", Sys.time(), "]"))
        d.spatial <- matrix(NaN, nrow = numCluster, ncol = numCluster)
        P.spatial <- matrix(1, nrow = numCluster, ncol = numCluster)
        adj.contact <- matrix(1, nrow = numCluster, ncol = numCluster)
        contact.dependent = FALSE
        contact.dependent.forced = FALSE
        contact.range = NULL
        contact.knn.k = NULL
        distance.use = NULL
        interaction.range = NULL
        ratio = NULL
        tol = NULL
        k.min = NULL
    }
    if (object@options$datatype == "RNA") {
        nLR1 <- nLR
    }
    else {
        if (contact.dependent.forced == TRUE) {
            cat("Force to run CellChat in a `contact-dependent` manner for all L-R pairs including secreted signaling.\n")
            P.spatial <- P.spatial * adj.contact
            nLR1 <- nLR
        }
        else {
            if (contact.dependent == TRUE && length(unique(pairLRsig$annotation)) > 0) {
                if (all(unique(pairLRsig$annotation) %in% c("Cell-Cell Contact"))) {
                  cat("All the input L-R pairs are `Cell-Cell Contact` signaling. Run CellChat in a contact-dependent manner. \n")
                  P.spatial <- P.spatial * adj.contact
                  nLR1 <- nLR
                }
                else if (all(unique(pairLRsig$annotation) %in% c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling"))) {
                  cat("Molecules of the input L-R pairs are diffusible. Run CellChat in a diffusion manner based on the `interaction.range`.\n")
                  nLR1 <- nLR
                }
                else {
                  cat("The input L-R pairs have both secreted signaling and contact-dependent signaling. Run CellChat in a contact-dependent manner for `Cell-Cell Contact` signaling, and in a diffusion manner based on the `interaction.range` for other L-R pairs. \n")
                  nLR1 <- max(which(pairLRsig$annotation %in% c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling")))
                }
            }
            else {
                cat("Run CellChat in a diffusion manner based on the `interaction.range` for all L-R pairs. Setting `contact.dependent = TRUE` if preferring a contact-dependent manner for `Cell-Cell Contact` signaling. \n")
                nLR1 <- nLR
            }
        }
    }
    Prob <- array(0, dim = c(numCluster, numCluster, nLR))
    Pval <- array(0, dim = c(numCluster, numCluster, nLR))
    set.seed(seed.use)
    permutation <- replicate(nboot, sample.int(nC, size = nC))
    data.use.avg.boot <- my.sapply(X = 1:nboot, FUN = function(nE) {
        groupboot <- group[permutation[, nE]]
        data.use.avgB <- aggregate(t(data.use), list(groupboot), FUN = FunMean)
        data.use.avgB <- t(data.use.avgB[, -1])
        return(data.use.avgB)
    }, simplify = FALSE)
    pb <- txtProgressBar(min = 0, max = nLR, style = 3, file = stderr())
    for (i in 1:nLR) {
        dataLR <- Matrix::crossprod(matrix(dataLavg[i, ], nrow = 1), matrix(dataRavg[i, ], nrow = 1))
        P1 <- dataLR^n/(Kh^n + dataLR^n)
        P1_Pspatial <- P1 * P.spatial
        if (sum(P1_Pspatial) == 0) {
            Pnull = P1_Pspatial
            Prob[, , i] <- Pnull
            p = 1
            Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
        }
        else {
            if (i > nLR1) {
                P.spatial <- P.spatial * adj.contact
            }
            if (is.element(i, index.agonist)) {
                data.agonist <- computeExpr_agonist(data.use = data.use.avg, pairLRsig, cofactor_input, index.agonist = i, Kh = Kh, n = n)
                P2 <- Matrix::crossprod(matrix(data.agonist, nrow = 1))
            }
            else {
                P2 <- matrix(1, nrow = numCluster, ncol = numCluster)
            }
            if (is.element(i, index.antagonist)) {
                data.antagonist <- computeExpr_antagonist(data.use = data.use.avg, pairLRsig, cofactor_input, index.antagonist = i, Kh = Kh, n = n)
                P3 <- Matrix::crossprod(matrix(data.antagonist, nrow = 1))
            }
            else {
                P3 <- matrix(1, nrow = numCluster, ncol = numCluster)
            }
            if (population.size) {
                P4 <- Matrix::crossprod(matrix(dataLavg2[i, ], nrow = 1), matrix(dataRavg2[i, ], nrow = 1))
            }
            else {
                P4 <- matrix(1, nrow = numCluster, ncol = numCluster)
            }
            Pnull = P1 * P2 * P3 * P4 * P.spatial
            Prob[, , i] <- Pnull
            Pnull <- as.vector(Pnull)
            Pboot <- my.sapply(X = 1:nboot, FUN = function(nE) {
                data.use.avgB <- data.use.avg.boot[[nE]]
                dataLavgB <- computeExpr_LR(geneL[i], data.use.avgB, complex_input)
                dataRavgB <- computeExpr_LR(geneR[i], data.use.avgB, complex_input)
                dataRavgB.co.A.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avgB, pairLRsig[i, , drop = FALSE], type = "A")
                dataRavgB.co.I.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avgB, pairLRsig[i, , drop = FALSE], type = "I")
                dataRavgB <- dataRavgB * dataRavgB.co.A.receptor/dataRavgB.co.I.receptor
                dataLRB = Matrix::crossprod(dataLavgB, dataRavgB)
                P1.boot <- dataLRB^n/(Kh^n + dataLRB^n)
                if (is.element(i, index.agonist)) {
                  data.agonist <- computeExpr_agonist(data.use = data.use.avgB, pairLRsig, cofactor_input, index.agonist = i, Kh = Kh, n = n)
                  P2.boot <- Matrix::crossprod(matrix(data.agonist, nrow = 1))
                }
                else {
                  P2.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
                }
                if (is.element(i, index.antagonist)) {
                  data.antagonist <- computeExpr_antagonist(data.use = data.use.avgB, pairLRsig, cofactor_input, index.antagonist = i, Kh = Kh, n = n)
                  P3.boot <- Matrix::crossprod(matrix(data.antagonist, nrow = 1))
                }
                else {
                  P3.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
                }
                if (population.size) {
                  groupboot <- group[permutation[, nE]]
                  dataLavg2B <- as.numeric(table(groupboot))/nC
                  dataLavg2B <- matrix(dataLavg2B, nrow = 1)
                  dataRavg2B <- dataLavg2B
                  P4.boot = Matrix::crossprod(dataLavg2B, dataRavg2B)
                }
                else {
                  P4.boot = matrix(1, nrow = numCluster, ncol = numCluster)
                }
                Pboot = P1.boot * P2.boot * P3.boot * P4.boot * P.spatial
                return(as.vector(Pboot))
            })
            Pboot <- matrix(unlist(Pboot), nrow = length(Pnull), ncol = nboot, byrow = FALSE)
            nReject <- rowSums(Pboot - Pnull > 0)
            p = nReject/nboot
            Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
        }
        setTxtProgressBar(pb = pb, value = i)
    }
    close(con = pb)
    Pval[Prob == 0] <- 1
    dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
    dimnames(Pval) <- dimnames(Prob)
    net <- list(prob = Prob, pval = Pval)
    execution.time = Sys.time() - ptm
    object@options$run.time <- as.numeric(execution.time, units = "secs")
    object@options$parameter <- list(type.mean = type, trim = trim, raw.use = raw.use, population.size = population.size, nboot = nboot, seed.use = seed.use, Kh = Kh, n = n, distance.use = distance.use, interaction.range = interaction.range, ratio = ratio, tol = tol, k.min = k.min, contact.dependent = contact.dependent, contact.range = contact.range, contact.knn.k = contact.knn.k, contact.dependent.forced = contact.dependent.forced)
    if (object@options$datatype != "RNA") {
        object@images$distance <- d.spatial
    }
    object@net <- net
    print(paste0(">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [", Sys.time(), "]"))
    return(object)
}
