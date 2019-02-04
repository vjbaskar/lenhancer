
#' getpeaks
getPeaks <- function(metafile, gene, distance = 100000){
    x <- metafile [ abs(metafile$distance) <= distance & metafile$gene == gene, ]
    return(x)
}

getExpression <- function(expression, gene){
    x <- expression [ gene, ]
    return(x)
}

getRegulation <- function(regulation_data, region_id){
    x <- regulation_data [rownames(regulation_data) %in% region_id, ]
}

get_coCRE <- function(gene_reg_tfbs, gene_reg_signal, gene_meta, min_tfb_events = 2, coCRE_corr_cutoff = 0.5, singleton_cutoff = 20000){
    require(igraph)
    require(arules)

    # Confirming binary input
    gene_reg_tfbs [ gene_reg_tfbs > 0 ] <- 1

    # Keep only regions with more than 2 tf binding events
    gene_reg_tfbs = as.matrix(gene_reg_tfbs)
    x = rowSums(gene_reg_tfbs)
    remRows = names(which(x < 2))

    # Get dice correlation of tfbs
    gene_reg_tfbs = as.matrix(gene_reg_tfbs)
    dice.matx = 1 - dissimilarity(gene_reg_tfbs, method = "dice")
    dice.matx = as.matrix(dice.matx)

    # Get spearman correlation of reg signal
    gene_reg_signal = as.matrix(gene_reg_signal)
    scor.matx = cor(t(gene_reg_signal), method="spearman")

    # Get the final correlation matrix
    final.matx = dice.matx * scor.matx
    final.matx[ final.matx <= coCRE_corr_cutoff ] <- 0
    final.matx [ remRows, remRows] <- 0
    final.matx = round(final.matx, 2)
    diag(final.matx) <- 1

    # Get cliques
    temp = graph.adjacency(final.matx, mode="undirected", weighted=TRUE)
    g <- temp
    E(g)$label <- E(g)$weight
    co <- layout_with_kk(g, weights = E(g)$weight , kkconst = 1000000)

    g.all = g
    g=delete.vertices(g,which(degree(g)<1))
    wcomm = cluster_fast_greedy(g)
    wcomm.mem = membership(wcomm)
    g.vert = V(g)$name
    x = gene_reg_tfbs[g.vert,]
    wcomm.all = cluster_fast_greedy(g.all)
    wcomm.all.mem = membership(wcomm.all)

    rn =c()
    pred.3d = c()
    predictors = t(gene_reg_signal)

    regions_class <- data.frame(regions = gene_meta$peak, distance = gene_meta$distance, class = as.character("singleton"))
    regions_class$class <- as.character(regions_class$class)


    for(i in unique(wcomm.all.mem)){
        cres = names(which(wcomm.all.mem==i))
        #CREs.merge.matx[,]
        x = predictors[, cres]
        add_data = 0
        if(!is.null(dim(x))){
            #  message(i)
            regions_class [ regions_class$regions %in% colnames(x), "class"] <- "coCRE"
            pred.3d = cbind(pred.3d, rowMeans(x))
            add_data = 1
        } else {
            #message("Orphan node ...")
            if(regions_class [ regions_class$regions == cres, "distance" ] <= singleton_cutoff){
                pred.3d = cbind(pred.3d, x)
                add_data = 1
            } else {
                add_data = 0
            }
        }
        if(add_data == 1){
            temp = cres
            temp = paste(temp, collapse="_")
            #rn = c(rn,paste("merge.",temp,collapse = "_",sep=""))
            rn = c(rn, temp)
        }
    }
    colnames(pred.3d) <- rn

    #### Singleton CREs
    singletons <- regions_class [ regions_class$class == "singleton", ]
    coCREs <- regions_class [ regions_class$class == "singleton", ]

    return(pred.3d)
}



