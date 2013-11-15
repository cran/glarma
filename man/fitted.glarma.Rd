\name{fitted.glarma}
\alias{fitted.glarma}
\alias{fitted.values.glarma}

\title{Extract GLARMA Model Fitted Values}

\description{

  \code{fitted}  method for class \code{"glarma"}. \code{fitted.values}
  is an alias for \code{fitted}.
}

\usage{
\method{fitted}{glarma}(object, ...)
}

\arguments{
   \item{object}{An object of class \code{"glarma"}, obtained from a
     call to \code{\link{glarma}}.}
   \item{...}{Further arguments passed to or from other methods.}
}

\details{
  This is an S3 generic function. \code{fitted} or \code{fitted.values}
  return the required fitted values from an object of class
  \code{"glarma"}.
}

\value{
  Fitted values \code{mu} extracted from the object \code{object}.
}

\seealso{
\code{\link{coef.glarma}}, \code{\link{residuals.glarma}},
\code{\link{glarma}}.
}

\keyword{Accessor Functions}