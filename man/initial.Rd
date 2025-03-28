\name{initial}
\alias{initial}

\title{Initial Parameter Generator for GLARMA from GLM}

\description{
   Function used to generate initial values of parameters for the
   GLARMA model from \code{\link{glm}} or \code{\link[MASS]{glm.nb}}.
}

\usage{
initial(y, X, offset = NULL, type = "Poi", alpha = 1)
}

\arguments{
  \item{y}{Numeric vector; response variable.}

  \item{X}{Matrix; explanatory variables. A vector of ones should be
    added to the data matrix as the first column for the
    \eqn{\beta}{beta} of the intercept.}

  \item{offset}{Either \code{NULL} or a numeric vector of length equal
    to the number of cases. Used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.}  

  \item{type}{Character; the distribution of the counts.}

  \item{alpha}{Numeric; an optional initial value for the
    \code{theta} parameter in the negative binomial distribution; the
    default value is 1.}

}

\details{
  Generates and returns the initial parameters for the GLARMA model
  under the specified distribution by fitting a generalized linear
  model.
}

\value{
  \item{beta}{A named numeric vector of initial coefficients.}
  \item{y}{If requested, the \eqn{y} vector used.}
  \item{X}{If requested, the model matrix.}
  \item{alpha}{The \code{theta} parameter in the negative binomial
    distribution returned by \code{\link[MASS]{glm.nb}}.  \code{NULL} if any
    other distribution is used.}
  \item{type}{The distribution of the counts in the GLARMA model.}
  \item{null.deviance}{Null deviance of the GLM with the same regression
    structure as the GLARMA model.}
  \item{df.null}{Null degrees of freedom of the GLM with the same
    regression structure as the GLARMA model.}

}

\author{
  "William T.M. Dunsmuir" <w.dunsmuir@unsw.edu.au>
}

\examples{
### Using the polio data
data(Polio)
y <- Polio[, 2]
X <- as.matrix(Polio[, 3:8])

glmMod <- initial(y, X, type = "Poi", alpha=1)
str(glmMod)
head(glmMod)
}

\keyword{Initial Parameter Generator}
