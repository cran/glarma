#####################################################################
#  Creates a GLARMA model object for simulation.
#
#
#EDIT:
#   User Manual:
#   Dunsmuir (2010) "R software for fitting observation driven
#   regression models for univariate time series"
#
#   Contributers: See Dunsmuir(2010).
#
#   Last Modified: March 29, 2011.
#
####################################################################
# NOTES on assumptions for manual entry:
# X is assumed to be at least a single column of ones of length sample size, n,
# for the intercept term.
#


glarmaSimModel <- function(X, beta, offset = NULL,
                           type = "Poi", mBin = NULL, alpha = NULL,
                           residType = "Pearson",
                           phiLags = NULL, phi = NULL,
                           thetaLags = NULL, theta = NULL)
{
  GSM <- list()
  n <- NROW(X)
  ## Test for correct values of arguments
  if ((!is.null(offset)) & (length(offset) != n))
    stop("When offset is NULL offset length must equal number of rows of X")
  if ((type == "Bin") & (length(mBin)!=n))
    stop("When type is 'Bin',  mBin length must equal number of rows of X")
  if (type == "NegBin") {
    if ((length(alpha) > 1) |(is.null(alpha)))
      stop("For negative binomial value of alpha is required")
    if (alpha < 0)
      stop("For negative binomial alpha must be positive")
  }
  if (length(phiLags) != length(phi))
    stop("lengths of phi and phiLags do not match")
  if (length(thetaLags) != length(theta))
    stop("lengths of theta and thetaLags do not match")
  if (NCOL(X) != length(beta))
    stop("Number of regression coefficients not equal to number of regressors")

  GSM$X <- X
  GSM$beta <- beta
  GSM$offset <- offset
  GSM$type <- type
  GSM$mBin <- mBin
  GSM$alpha <- alpha
  GSM$residType <- residType
  GSM$phiLags <- phiLags
  GSM$phi <- phi
  GSM$thetaLags <- thetaLags
  GSM$theta <- theta

  class(GSM)  <- "glarmaSimModel"
  GSM
}
