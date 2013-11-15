\name{extractAIC.glarma}
\alias{extractAIC.glarma}

\title{Extract AIC from a GLARMA Model}

\description{
  \code{extractAIC} method for class \code{"glarma"}. Used to
  extract AIC from a \code{glarma} object.
}

\usage{
\method{extractAIC}{glarma}(fit, ...)
}

\arguments{
  \item{fit}{An object of class \code{"glarma"}, obtained from a call to
     \code{\link{glarma}}.}
  \item{...}{Further arguments passed to or from other methods.}
}

\value{
  AIC extracted from \code{object}
}

\seealso{
\code{\link{coef.glarma}}, \code{\link{residuals.glarma}},
\code{\link{glarma}}.
}

\keyword{Accessor Functions}
