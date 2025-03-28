\name{glarmaSimModel}
\alias{glarmaSimModel}

\title{
  Create a glarma simulation model
}
\description{
  Defines a glarma model for input to the function \code{glarmaSim}.
}
\usage{
glarmaSimModel(X, beta, offset = NULL, type = "Poi", mBin = NULL,
               alpha = NULL, residType = "Pearson",
               phiLags = NULL, phi = NULL,
               thetaLags = NULL, theta = NULL)
}

\arguments{
  \item{X}{Matrix; the explanatory variables. A vector of ones should be
    added to the data matrix as the first column for the \code{beta} of
    the intercept.}
  \item{beta}{Numeric vector; values of the regression coefficients.}
  \item{offset}{Either \code{NULL} or a numeric vector of length equal
    to the number of cases.}
  \item{type}{Character; the count distribution. Possible values are
    \code{"Poi"} (Poisson), \code{"Bin"} (binomial) and \code{"NegBin"}
    (negative binomial). The default is the Poisson distribution.}
  \item{mBin}{Numeric vector of length equal to the forecast horizon;
    only for the binomial case, the number of trials for each time point.}
  \item{alpha}{Numeric; for the negative binomial case, the shape
    parameter for \code{\link[MASS]{glm.nb}}}
  \item{residType}{Character; the type of residuals to be used. Possible
    values are \code{"Pearson"} and \code{"Score"}, and for the binomial
    distribution \code{"Identity"} is also allowed. The default is to
    use Pearson residuals.}
  \item{phiLags}{Numeric vector; AR orders.}
  \item{phi}{Numeric vector; values for the corresponding
    AR orders.}
  \item{thetaLags}{Numeric vector; MA orders.}
  \item{theta}{Numeric vector; values for the corresponding
    MA orders.}
}

\value{
  Creates a list object of type "glarmaSimModel" with the same list
  elements as in the function call.  
}

\author{
  "David J. Scott" <d.scott@auckland.ac.nz> and
   "William T.M. Dunsmuir" <w.dunsmuir@unsw.edu.au>
}

\seealso{
See \code{\link{glarmaSim}}.
}
\examples{
### Test glarmaSimModel
data("Polio")
y <- Polio[, 2]
X <- as.matrix(Polio[, 3:8])
glarmamod <- glarma(y, X, thetaLags = c(1,2,5), type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)
summary(glarmamod)
PolioModel <- glarmaSimModel(X, beta = coef.glarma(glarmamod, type = "beta"),
                             phiLags = glarmamod$phiLags,
                             phi = rep(0.1, length(glarmamod$phiLags)),
                             thetaLags = glarmamod$thetaLags,
                             theta = rep(0.1, length(glarmamod$thetaLags)),
                             type = glarmamod$type,
                             residType = glarmamod$residType)
str(PolioModel)
coef.glarma(glarmamod, type = "ARMA")
PolioModel <- extractGlarmaSimModel(glarmamod)
str(PolioModel)

}

\keyword{ts }
