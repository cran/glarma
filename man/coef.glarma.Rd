\name{coef.glarma}
\alias{coef.glarma}
\alias{coefficients.glarma}

\title{Extract GLARMA Model Coefficients}

\description{
\code{coef} is a generic function which extracts GLARMA model coefficients
from objects returned by modeling functions. \code{coefficients} is an alias
for it.
}

\usage{
\method{coef}{glarma}(object, types = "all", ...)
}

\arguments{
  \item{object}{An object of class \code{"glarma"}, obtained from a call
    to \code{\link{glarma}}.}

  \item{types}{Character; which coefficients to extract, either
    \code{ARMA}), \code{beta}, \code{NB} or \code{all}. The default is
    \code{all}.}
    
  \item{...}{Further arguments passed to or from other methods.}
}

\details{  
  This is an S3 generic function. \code{coef} or \code{coefficients}
  return the requested coefficients from the object of class
  \code{"glarma"}. By changing the argument \code{type}, either the ARMA
  coefficients (\code{ARMA}), regression coefficients (\code{beta}) or
  all coefficients are returned. In the case of negative binomial
  counts, the negative binomial coefficient \eqn{\alpha}{alpha} is also
  returned if type is \code{all}, or if type is \code{NB}. The default
  is \code{all}.
}

\value{ 
  \code{ARMA} coefficients, \code{beta} coefficients, \code{NB}
  coefficients or all of these three types of coefficients are extracted
  from the \code{glarma} model object \code{object}.

  A named numeric vector or list of named numeric vectors is returned.
}


\seealso{
  \code{\link{fitted.glarma}} and \code{\link{residuals.glarma}} for
  related methods;
}

\examples{
data(Polio)
Y <- Polio[, 2]
X <- as.matrix(Polio[, 3:8])
glarmamod <- glarma(Y, X, thetaLags = c(1, 2, 5), type = "Poi",
                    method = "FS", residuals= "Pearson",
                    maxit = 100, grad = 1e-6)

coef(glarmamod, type = "ARMA")
coef(glarmamod, type = "beta")
coef(glarmamod, type = "all")
}

\keyword{Accessor Functions}

