\name{likTests}
\alias{likTests}
\alias{likeTests}
\alias{print.likTests}

\title{
  Likelihood Ratio Test and Wald Test for \code{GLARMA} Fit
}

\description{

  Function to carry out the likelihood ratio and Wald tests of serial
  dependence when the alternative is a GLARMA process.  This function
  takes a \code{"glarma"} object and uses its attributes to set up a GLM
  fit that matches the GLARMA model regression structure.  This is done
  to ensure that the GLM object is the null hypothesis for testing
  against the \code{"glarma"} object.

}

\usage{
likTests(object)
likeTests(object)
\method{print}{likTests}(x, ...)
}

\arguments{
  \item{object}{An object of class \code{"glarma"}, obtained from a call to
    \code{\link{glarma}}.}

  \item{x}{An object of class \code{"likTests"}, a result of a call to
    \code{likTests}}

  \item{...}{Further arguments passed to or from other methods.}
 }

\details{
  This function carries out the likelihood ratio and Wald tests
  for comparing the null model and the alternative model.

  \code{likeTests} is an alias for \code{likTests}.
}

\value{
  \code{likTests} returns an object of class \code{"likTests"}. A matrix is
  shown with the statistics and p-value for each test. The significance
  stars alongside help to identify any probabilities less than 0.05 or
  0.01.
}


\author{
  "William T.M. Dunsmuir" <w.dunsmuir@unsw.edu.au>
}

\examples{
### Binomial response
data(Polio)
y <- Polio[, 2]
X <- as.matrix(Polio[, 3:8])
glarmamod <- glarma(y, X, thetaLags = c(1,2,5), type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 2.22e-16)
likTests(glarmamod)
likeTests(glarmamod)
### Negative binomial response
glarmamod <- glarma(y, X, thetaLags = c(1,2,5), type = "NegBin", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 2.22e-16)
glarmamod
summary(glarmamod)
likeTests(glarmamod)

### Negative Binomial Response with offset
glarmamod.offset <- glarma(y, X, thetaLags = c(1,2,5), type = "NegBin",
                           method = "FS", offset=log(X[,1]),
                           residuals = "Pearson", maxit = 100, grad = 2.22e-16)
glarmamod.offset
summary(glarmamod.offset)
likeTests(glarmamod.offset)
}

\keyword{Diagnostic}