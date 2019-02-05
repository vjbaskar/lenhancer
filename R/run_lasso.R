
run_lasso <- function(gene_coCRE_signal, gene_expression, ltype= "l1se", alphaVal = 1, scale.predictors = TRUE, family = "gaussian", nfold = nfold){
   # require(glmnet)
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
