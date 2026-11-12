function (object, sources.use = NULL, targets.use = NULL, signaling = NULL, 
    pairLR.use = NULL, remove.isolate = TRUE, thresh = 0.05, 
    return.object = TRUE) 
{
    net <- object@net
    if (is.null(sources.use) & is.null(targets.use) & is.null(signaling) & 
        is.null(pairLR.use)) {
        prob <- net$prob
        pval <- net$pval
        pval[prob == 0] <- 1
        prob[pval >= thresh] <- 0
        net$count <- apply(prob > 0, c(1, 2), sum)
        net$weight <- apply(prob, c(1, 2), sum)
        net$weight[is.na(net$weight)] <- 0
        net$count[is.na(net$count)] <- 0
    }
    else {
        df.net <- subsetCommunication(object, slot.name = "net", 
            sources.use = sources.use, targets.use = targets.use, 
            signaling = signaling, pairLR.use = pairLR.use, thresh = thresh)
        df.net$source_target <- paste(df.net$source, df.net$target, 
            sep = "_")
        df.net2 <- df.net %>% group_by(source_target) %>% summarize(count = n(), 
            .groups = "drop")
        df.net3 <- df.net %>% group_by(source_target) %>% summarize(prob = sum(prob), 
            .groups = "drop")
        df.net2$prob <- df.net3$prob
        a <- stringr::str_split(df.net2$source_target, "_", simplify = T)
        df.net2$source <- as.character(a[, 1])
        df.net2$target <- as.character(a[, 2])
        cells.level <- levels(object@idents)
        if (remove.isolate) {
            message("Isolate cell groups without any interactions are removed. To block it, set `remove.isolate = FALSE`")
            df.net2$source <- factor(df.net2$source, levels = cells.level[cells.level %in% 
                unique(df.net2$source)])
            df.net2$target <- factor(df.net2$target, levels = cells.level[cells.level %in% 
                unique(df.net2$target)])
        }
        else {
            df.net2$source <- factor(df.net2$source, levels = cells.level)
            df.net2$target <- factor(df.net2$target, levels = cells.level)
        }
        count <- tapply(df.net2[["count"]], list(df.net2[["source"]], 
            df.net2[["target"]]), sum)
        prob <- tapply(df.net2[["prob"]], list(df.net2[["source"]], 
            df.net2[["target"]]), sum)
        net$count <- count
        net$weight <- prob
        net$weight[is.na(net$weight)] <- 0
        net$count[is.na(net$count)] <- 0
    }
    if (return.object) {
        object@net <- net
        return(object)
    }
    else {
        return(net)
    }
}
# <bytecode: 0x61153cd8dc90>
# <environment: namespace:CellChat>
