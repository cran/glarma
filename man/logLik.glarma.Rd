\name{logLik.glarma}
\alias{logLik.glarma}

\title{Extract Log-Likelihood from GLARMA Models}

\description{
  \code{logLik} is a generic function which extracts the GLARMA model
  log-likelihood from objects returned by modeling functions.
}

\usage{
\method{logLik}{glarma}(object, deriv, ...)
}

\arguments{
  \item{object}{An object of class \code{"glarma"}, a result of a call
    to \code{\link{glarma}}.}
  \item{deriv}{Numeric; either "0", "1"
    or "2". It is used to choose and extract the log-likehood, its
    derivative or its second derivative respectively from the
    \code{"glarma"} object. The default is "0".}
  \item{...}{Further arguments passed to or from other methods.}

}

\details{
  This is an S3 generic function. \code{logLik} returns the
  log-likelihood, its derivative, or its second derivative from the
  object of class \code{glarma} based on the value of the argument
  \code{deriv}. "0" is for the log-likelihood, "1" is for the first
  derivative of log-likelihood and "2" is for the second derivative of
  the log-likelihood.
}

\value{
  The log-likelihood, the derivative of the log-likelihood or the second
  derivative of the log-likelihood extracted from the GLARMA model
  object \code{object}.
}

\seealso{
\code{\link{coef.glarma}}, \code{\link{residuals.glarma}},
\code{\link{fitted.glarma}}, \code{\link{glarma}}.
}

\examples{
data(Polio)
Y <- Polio[, 2]
X <- as.matrix(Polio[, 3:8])
glarmamod <- glarma(Y, X, thetaLags = c(1, 2, 5), type = "Poi", method ="FS",
                    residuals = "Pearson", maxit = 100 , grad = 1e-6)

logLik(glarmamod, deriv = 0)
logLik(glarmamod, deriv = 1)
logLik(glarmamod, deriv = 2)
}

\keyword{Accessor Functions}
