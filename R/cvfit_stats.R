

getShortListed <- function(cvfit, lambda, predictors){
    shortlisted = coef(cvfit,s=lambda)
    colnames(shortlisted) <- "lassoval"
    shortlisted = as.data.frame(as.matrix(shortlisted))
    shortlisted

    reqCols = rownames(shortlisted)
    reqCols = reqCols[which(!grepl("Intercept",reqCols))]
    reqCols

    shortlisted = shortlisted[reqCols,]
    names(shortlisted) = reqCols
    shortlisted = shortlisted[shortlisted != 0]
    shlst = names(shortlisted)
    temp = c("expr", shlst)
    final.shortlisted = as.data.frame(predictors[,shlst])
    colnames(final.shortlisted) = shlst
    return(final.shortlisted)
}

getPredictedCor <- function(cvfit.obj, xmat, y, lambda = "lambda.1se"){
    lambdaVal = cvfit.obj[lambda]
    #    message("LambdaVal = ",lambdaVal )
    y.pred = predict.cv.glmnet(cvfit.obj,xmat , s=lambda)
    pcoef = cor.test(y,y.pred, method="spearman")
    temp = list(x = xmat, y.pred = y.pred, pcoef=pcoef, cvfit = cvfit.obj)
    return(temp)

}
