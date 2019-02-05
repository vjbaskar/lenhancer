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
