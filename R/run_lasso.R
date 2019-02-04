lasso.glmnet <- function(x.lasso, y.lasso, alphaVal=1, family = "gaussian", scale.predictors = TRUE,  nfold=10){

    lasso.mod=glmnet(x.lasso,y.lasso, alpha=alphaVal, family = family, standardize = scale.predictors)
    # Cross-validation = n fold
    cvfit = cv.glmnet(x.lasso, y.lasso, alpha=alphaVal, family = family, nfold=nfold, standardize = scale.predictors,  nlambda =10000, type.gaussian="naive", lambda = seq(0.0001, 10, 0.001))

    lmin = cvfit$lambda.min
    l1se = cvfit$lambda.1se

    df = list( lasso.mod = lasso.mod,
               cvfit = cvfit,
               lmin = lmin,
               l1se = l1se
    )
    return(df)
}

cvfit.pickOne <- function(cvfit){
    temp = cvfit$nzero
    usort = sort(unique(temp))
    usort = usort[ usort > 0]
    totalPred = usort[1]
    nzero.1 = which(temp==totalPred)[1]
    lambdaval = cvfit$lambda[nzero.1]
    return(lambdaval)
}

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

run_lasso <- function(gene_coCRE_signal, gene_expression, ltype= "l1se", alphaVal = 1, scale.predictors = TRUE, family = "gaussian", nfold = nfold){
    require(glmnet)
    predictors = x = as.matrix(gene_coCRE_signal)
    response = y = as.matrix(gene_expression)
    lasso.glmnet.outpt = lasso.glmnet(predictors, response, alphaVal = alphaVal, family, scale.predictors = scale.predictors, nfold =  nfold)

    lasso.mod = lasso.glmnet.outpt$lasso.mod
    cvfit = lasso.glmnet.outpt$cvfit
    lmin = lasso.glmnet.outpt$lmin
    l1se = lasso.glmnet.outpt$l1se
    lambda = lmin

    final.shortlisted.1se = getShortListed(cvfit,l1se, predictors )
    final.shortlisted.lmin = getShortListed(cvfit,lmin, predictors )

    pickOne.l1se = F; pickOne.lmin=F
    if(ncol(final.shortlisted.1se)==0){
        message("1se has zero predictors. Resorting to the best")
        l1se = cvfit.pickOne(cvfit)
        pickOne.l1se = T
    }

    if(ncol(final.shortlisted.lmin)==0){
        message("lmin has zero predictors. Resorting to the best")
        lmin = cvfit.pickOne(cvfit)
        pickOne.lmin = T
    }

    final.shortlisted.1se = getShortListed(cvfit,l1se, predictors )
    final.shortlisted.lmin = getShortListed(cvfit,lmin, predictors )

    cor.1se = getPredictedCor(cvfit, x, y, lambda = l1se)$pcoef
    cor.lmin = getPredictedCor(cvfit, x, y, lambda = lmin)$pcoef

    model.pvalue =getPValue(predictors, response)
    rlist = list()
    rlist[["l1se_pred"]] <- final.shortlisted.1se
    rlist[["lmin_pred"]] <- final.shortlisted.lmin
    rlist[["p"]] <- model.pvalue
    return(rlist)

}
