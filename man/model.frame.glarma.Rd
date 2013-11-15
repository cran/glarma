\name{model.frame.glarma}
\alias{model.frame.glarma}

\title{Extracting the Model Frame of the GLARMA Model}

\description{
  \code{model.frame} (a generic function) and its methods return a
  data frame with the variables that are used in the
  \code{\link{glarma}} model.
}

\usage{
\method{model.frame}{glarma}(formula, ...)
}

\arguments{
  \item{formula}{An object of class \code{glarma}, obtained from a call
    to \code{\link{glarma}}.}
  \item{...}{Further arguments passed to or from other methods.}
}

\details{
  This is an S3 generic function. It extracts the response variable
  vector and the matrix of the explanatory variables from the object of
  class \code{"glarma"}, and combines them as a data frame.
}

\value{
  A data frame with the variables used in the fitted
  \code{\link{glarma}} model.
}

\author{
  Cenanning Li <cli113@aucklanduni.ac.nz>
}


\seealso{
\code{\link{coef.glarma}}, \code{\link{residuals.glarma}},
\code{\link{fitted.glarma}}, \code{\link{glarma}}.
}

\examples{
data(Polio)
print(y <- Polio[, 2])
X <- as.matrix(Polio[, 3:8])
str(X)
head(X)


glarmamod <- glarma(y, X, thetaLags = c(1, 2, 5), type = "Poi",
                    method = "FS", residuals = "Pearson",
                    maxit = 100, grad = 1e-6)

str(model.frame(glarmamod))
head(model.frame(glarmamod))
}

\keyword{Accessor Functions}
