\name{paramGen}
\alias{paramGen}
\alias{deltaGen}
\alias{thetaGen}
\alias{phiGen}

\title{Parameter Generators}

\description{
  Functions which use the arguments of a \code{\link{glarma}} call to
  generate the initial \code{delta}, \code{theta} and \code{phi}
  vectors.
}

\usage{
deltaGen(y, X, offset = NULL, phiInit, thetaInit, type, alpha,
         beta, alphaInit)
thetaGen(thetaLags, thetaInit)
phiGen(phiLags, phiInit)
}

\arguments{
  \item{y}{Numeric vector; response variable.}

  \item{X}{Matrix; the explanatory variables. A vector of ones should be
    added to the data matrix as the first column for the \code{beta} of
    the intercept.}

  \item{offset}{Either \code{NULL} or a numeric vector of length equal
    to the number of cases. Used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.}

  \item{phiInit}{Numeric vector; initial values for the
    corresponding AR orders.}

  \item{thetaInit}{Numeric vector; initial values for the
    corresponding orders.}

  \item{type}{Character; the count distribution. The default is the
    Poisson distribution.}

   \item{beta}{Numeric vector; initial values of the parameters of
    variables. It is for the user to construct the specific \code{delta}
    by themselves.}

  \item{alpha}{Numeric; an optional initial \code{theta} parameter for
    \code{\link[MASS]{glm.nb}}.}

  \item{alphaInit}{Numeric; an initial \code{theta} parameter for
    \code{glarma} for negative binomial counts.}

  \item{thetaLags}{Numeric vector; MA orders}

  \item{phiLags}{Numeric vector; AR orders}
}



\details{
  The \code{thetaGen} and \code{phiGen} functions take the arguments,
  \code{thetaLags}, \code{phiLags}, \code{thetaInit} and \code{phiInit},
  in a \code{\link{glarma}} call to generate and return the initial
  \code{theta} and \code{phi} vectors with orders corresponding to their
  names. Then the \code{deltaGen} function uses the values returned by
  \code{thetaGen}, \code{phiGen} and other arguments in the
  \code{\link{glarma}} call to generate and return the initial
  \code{delta} vector with correct names.
}

\value{
  \code{thetaGen} returns a list containing \code{thetaLags} and
  \code{thetaInit}. \code{thetaInit} is the initial \code{theta} vector
  with its corresponding MA orders as its names.

  \code{phiGen} returns a list containing \code{phiLags} and
  \code{phiInit}. \code{phiInit} is the initial \code{phi} vector
  with its corresponding MA orders as its names.

  \code{deltaGen} returns a named vector giving the values of
  \code{beta}, \code{phiInit}, \code{thetaInit} plus \code{alpha} in the
  negative binomial case.
}

\author{
  "Cenanning Li" <cli113@aucklanduni.ac.nz> and
  "William T.M. Dunsmuir" <w.dunsmuir@unsw.edu.au>
}

\examples{
### Using the polio data
data(Polio)
y <- Polio[, 2]
X <- as.matrix(Polio[, 3:8])

## generate the theta vector
theta.lags <- c(1, 2, 5)
theta.init <- c(0.0, 0.0, 0.0)

theta <- thetaGen(theta.lags, theta.init)
print(thetaLags <- theta[[1]])
print(theta.init <- theta[[2]])

## generate the vector of phi
phi.lags <- rep(0, 0)
phi.init <- rep(0, 0)
phi <- phiGen(phi.lags, phi.init)
print(phiLags <- phi[[1]])
print(phi.init <- phi[[2]])

## generate the delta vector
delta <- deltaGen(y = y, X = X, phiInit = phi.init,
                  thetaInit = theta.init, type = "Poi",
                  alpha = 1)
delta
}

\keyword{Initial Parameter Generators}