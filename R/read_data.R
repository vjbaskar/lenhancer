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

