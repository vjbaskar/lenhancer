
cvfit.pickOne <- function(cvfit){
    temp = cvfit$nzero
    usort = sort(unique(temp))
    usort = usort[ usort > 0]
    totalPred = usort[1]
    nzero.1 = which(temp==totalPred)[1]
    lambdaval = cvfit$lambda[nzero.1]
    return(lambdaval)
}
