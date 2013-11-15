### Function to calculate perform GLARMA fit and compute LRT and Wald Tests of
### serial dependence when the alternative is a GLARMA process.  This function
### takes a GLARMA object and uses its attributes to set up a GLM fit that
### matches the GLARMA model regression structure.  This is done to ensure that
### the GLM object is the null hypothesis for testing against the GLARMA object.
### Written by W.T.M. Dunsmuir Last Modified April 7 2011.

likeTests <- likTests <- function(object) {

    ## Calculate Null Hypothesis values of likelihood.

    r <- object$r
    pq <- object$pq

    y <- object$y
    X <- object$X
    type <- object$type

    if (type == "Poi")
        (GLM <- glm(y ~ -1 + X, family = poisson, x = T))


    if (type == "Bin")
        (GLM <- glm(y ~ -1 + X, family = binomial(link = logit),
                na.action = na.omit, x = T))

    if (type == "NegBin")
        ## fix up neg bin case later
        (GLM <- glm.nb(y ~ -1 + X, init.theta = 1, x = T))


    ll.null <- r - GLM$aic/2


    ## Calculate LRT and WALD against alternative hypothesis

    LRT <- NA
    LRT.P <- NA
    Wald <- NA
    Wald.P <- NA

    if (object$errCode == 0 & object$WError == 0) {

        ## LRT
        ll.alt <- object$logLik
        LRT <- 2 * (ll.alt - ll.null)
        LRT.P <- 1 - pchisq(LRT, pq)

        ## Wald
        thetahat <- object$delta[r + 1:pq]
        Wald <- (thetahat) %*%
                solve(object$cov[r + 1:pq, r + 1:pq]) %*%
                thetahat
        Wald.P <- 1 - pchisq(Wald, pq)

    }
    likTests <- matrix(c(LRT, Wald, LRT.P, Wald.P), ncol = 2)
    colnames(likTests) <- c("Statistic", "p-value")
    rownames(likTests) <- c("LR Test", "Wald Test")
    class(likTests) <- "likTests"
    likTests
}

print.likTests <- function(x, ...){
  printCoefmat(x, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n")
}
