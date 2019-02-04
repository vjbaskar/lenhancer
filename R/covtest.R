
getPValue <- function(predictors, response){

    require(covTest)
    require(selectiveInference)

    # Predictors and response
    x = scale(predictors,T,T)
    cn = colnames(x)
    cn = cn [ ! grepl("temp", cn)]
    x = as.matrix(x[,cn])
    colnames(x) = cn
    n = nrow(x)
    if(ncol(x) == 1){
        temp = cbind(rep(0,n), rep(1,n))
        colnames(temp) <- c("temp1","temp2")
        x = cbind(x,temp)
    }
    x = as.matrix(x)
    y = as.matrix(response)
    dim(x)
    dim(y)
    # Estimate sigma - sd of noise. sqrt(sum((y-yhat)^2) / (n-df-1))
    ## x is scaled. Therefore not using standardisation or normalisation
    sigma = estimateSigma(x, y, standardize = FALSE)$sigmahat
    sigma
    # covtest using an estimated sigma
    a=lars(x,y, type = "lasso", normalize = FALSE)
    ctest = covTest(a,x,y, sigma.est = sigma)
    # Extract p-values
    res = ctest$results
    pvals = res[,3]
    delV = res[,2]
    cn = colnames(x)
    df = data.frame(dropInVariance = delV, pvalue = pvals)
    df[ is.na(df$pvalue), "pvalue"] <- 1
    df = df[ !is.na(df$pvalue), ]
    model.pvalue = min(df$pvalue)
    return( model.pvalue )
}
