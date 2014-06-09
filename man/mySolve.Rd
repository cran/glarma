\name{mySolve}

\alias{mySolve}

\title{Matrix Inversion of the Hessian of the Log-Likelihood}

\description{
  Inverts the second derivative matrix of the log-likelihood to obtain
  the estimated covariance matrix of the parameters.
 }

\usage{
mySolve(A)
}

\arguments{
  \item{A}{Matrix; the negative second derivative of the
    log-likelihood}
}

\details{
  \code{mySolve} attempts to invert its matrix argument. If the matrix
  supplied is not invertible, \code{ErrCode} is set to 1.
}

\value{

  \item{Ainv}{inverse of the negative second derivative of the
    loglikelihood. If the inverse is unable to be obtained,
    returns the original negative second derivative of the
    log-likelihood.}

  \item{ErrCode}{Numeric; 0 if the inverse can be found, 1 if not.}

}

\author{
  "William T.M. Dunsmuir" <w.dunsmuir@unsw.edu.au>
}


\examples{
### Using the polio data
data(Polio)
y <- Polio[, 2]
X <- as.matrix(Polio[, 3:8])

## Construct the vectors of phi lags and theta lags
theta.lags <- c(1, 2, 5)
phi.lags <- rep(0, 0)
## Construct the initial delta vector
delta <- c("Intcpt" = 0.2069383, "Trend" = -4.7986615 ,
           "CosAnnual" = -0.1487333, "SinAnnual" = -0.5318768,
           "CosSemiAnnual" = 0.1690998, "SinSemiAnnual" = -0.4321435,
           "theta_1" = 0, "theta_2"= 0, "theta_5"= 0 )

## Calculate the second derivative of the loglikelihood
glarmamod <- glarmaPoissonPearson(y, X, delta = delta, phiLags = phi.lags,
                           thetaLags = theta.lags, method = "FS")

## estimate the covariance matrix of the estimators from the second
## derivative of the loglikelihood
mySolve(-glarmamod$ll.dd)
}

\keyword{GLARMA}