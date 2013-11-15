\name{residuals.glarma}
\alias{residuals.glarma}
\alias{resid.glarma}

\title{Extract GLARMA Model Residuals}

\description{
  \code{residuals} is a generic function which extracts model residuals
  from objects returned by the modeling function
  \code{\link{glarma}}. \code{resid} is an alias for \code{residuals}.
}

\usage{
\method{residuals}{glarma}(object, ...)
}

\arguments{
 \item{object}{An object of class \code{"glarma"}, a result of a call
   to \code{\link{glarma}}.}
 \item{...}{Further arguments passed to or from other methods.}
}

\value{
  Residuals extracted from the object \code{object}.
}

\author{
  "William T.M. Dunsmuir" <w.dunsmuir@unsw.edu.au> and
  "Cenanning Li" <cli113@aucklanduni.ac.nz>
}

\seealso{
\code{\link{coefficients.glarma}}, \code{\link{fitted.glarma}},
\code{\link{glarma}}.
}

\keyword{Accessor Functions}
