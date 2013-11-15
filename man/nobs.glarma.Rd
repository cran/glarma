\name{nobs.glarma}

\alias{nobs.glarma}

\title{Extract the Number of Observations from a GLARMA Model Fit}

\description{
  An accessor function used to extract the number of observations from a
  \code{"glarma"} object.
}

\usage{
\method{nobs}{glarma}(object, ...)
}

\arguments{
  \item{object}{An object of class \code{"glarma"}, obtained from a call
    to \code{\link{glarma}}.}
  \item{...}{Further arguments passed to or from other methods.}
}

\value{
  The number of observations extracted from the object \code{object}.
}

\author{
  "Cenanning Li" <cli113@aucklanduni.ac.nz>
}

\seealso{
\code{\link{coef.glarma}}, \code{\link{residuals.glarma}},
\code{\link{fitted.glarma}}, \code{\link{glarma}}.
}
\examples{
### Example from Davis, Dunsmuir Wang (1999)
## MA(1,2,5), Pearson Residuals, Fisher Scoring
data(Polio)
y <- Polio[, 2]
X <- as.matrix(Polio[, 3:8])
glarmamod <- glarma(y, X, thetaLags = c(1,2,5), type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 2.22e-16)

nobs(glarmamod)
}

\keyword{Accessor Functions}
