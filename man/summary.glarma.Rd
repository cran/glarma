\name{summary.glarma}
\alias{summary.glarma}
\alias{print.summary.glarma}
\alias{glarmaModelEstimates}

\title{Summarize GLARMA Fit}

\description{
  \code{summary} method for class \code{glarma} and functions to
  generate the estimates for this \code{summary} method.
}

\usage{
\method{summary}{glarma}(object, tests = TRUE, ...)
\method{print}{summary.glarma}(x, digits = max(3L, getOption("digits") - 3L), ...)
glarmaModelEstimates(object)
}

\arguments{
 \item{object}{An object of class \code{"glarma"}, obtained from a call
   to \code{\link{glarma}}.}

 \item{x}{An object of class \code{"summary.glarma"}, obtained from a
   call to \code{summary.glarma}.}

 \item{digits}{Numeric; minimum number of significant digits to be used
 for most numbers.}

 \item{tests}{Logical; if \code{TRUE}, the likelihood-ratio test and
 the Wald test are shown in the summary. The default is \code{TRUE}.}

 \item{...}{Further arguments passed to or from other methods.}
}

\value{
  \code{summary.glarma} returns an object of class
  \code{"summary.glarma"}, a list with components
  \item{call}{the component from \code{object}}
  \item{null.deviance}{null deviance of the GLM with the same regression
    structure as the GLARMA model.}
  \item{df.null}{null degrees of freedom of the GLM with the same
    regression structure as the GLARMA model.}
  \item{phi.lags}{the component from \code{object}.}
  \item{theta.lags}{the component from \code{object}.}
  \item{pq}{the component from \code{object}.}
  \item{iter}{the component from \code{object}.}
  \item{deviance}{the deviance of the fitted model.}
  \item{df.residual}{the degrees of freedom of the fitted model.}
  \item{deviance.resid}{the component from \code{object}.}
  \item{aic}{the component from \code{object}.}
  \item{methods}{vector specifying the count distribution of the
    GLARMA model, the iteration method and the type of residual used.}
  \item{tests}{whether tests were asked for.}
  \item{likTests}{if \code{tests} is \code{TRUE}, the result of a call
    to \code{likTests}, \code{NULL} otherwise.}
  \item{coefficients1}{the matrix of beta coefficients, standard
    errors, z-ratio and p-values.}
  \item{coefficients2}{the matrix of ARMA coefficients, standard
    errors, z-ratio and p-values.}
  \item{coefficients3}{when the count distribution is negative binomial,
  a matrix with 1 row, giving the negative binomial parameter, its
  standard error, z-ratio and p-value.}
}

\author{
  "William T.M. Dunsmuir" <w.dunsmuir@unsw.edu.au> and
  "Cenanning Li" <cli113@aucklanduni.ac.nz>
}

\seealso{\code{\link{glarma}}, \code{\link{summary}}.}

\examples{
## For examples see example(glarma)
}

\keyword{Print}