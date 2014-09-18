\name{PIT}
\alias{glarmaPredProb}
\alias{glarmaPIT}

\title{Non-randomized Probability Integral Transformation}

\description{
Functions to produce the non-randomized probability integral
transform (PIT) to check the adequacy of the distributional assumption of
the GLARMA model.
}


\usage{
glarmaPredProb(object)
glarmaPIT(object, bins = 10)
}


\arguments{
  \item{object}{An object of class \code{"glarma"}, obtained from a
    call to \code{\link{glarma}}.}

  \item{bins}{Numeric; the number of bins used in the
    PIT.}

}
\details{
These functions are used for the assessment of predictive distributions
in discrete data. They obtain the predictive probabilities
and the probability integral transformation for a fitted GLARMA
model.
}

\value{
  \code{glarmaPredProb} returns a list with values:

  \item{upper}{the predictive cumulative probabilities used
     as the upper bound for computing the non-randomized PIT.}

   \item{lower}{the predictive cumulative probabilities used
     as the lower bound for computing the non-randomized PIT.}

  \code{glarmaPIT} returns a list with values:

  \item{upper}{the predictive cumulative probabilities used
     as the upper bound for computing the non-randomized PIT.}

   \item{lower}{the predictive cumulative probabilities used
     as the lower bound for computing the non-randomized PIT.}

   \item{conditionalPIT}{the conditional probability integral
     transformation given the observed counts.}

   \item{PIT}{the probability integral transformation.}
}

\references{
Czado, Claudia and Gneiting, Tilmann and Held, Leonhard (2009)
Predictive model assessment for count data. \emph{Biometrics},
\bold{65}, 1254--1261.

Jung, Robert.C and Tremayne, A.R (2011) Useful models for time series
of counts or simply wrong ones? \emph{Advances in Statistical
Analysis}, \bold{95}, 59--91.
}

\author{
  "David J. Scott" <d.scott@auckland.ac.nz> and
  "Cenanning Li" <cli113@aucklanduni.ac.nz>
}

\examples{
### Example from Davis, Dunsmuir Wang (1999)
## MA(1,2,5), Pearson Residuals, Fisher Scoring
data(Polio)
y <- Polio[, 2]
X <- as.matrix(Polio[, 3:8])
glarmamod <- glarma(y, X, thetaLags = c(1,2,5), type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 2.22e-16)
glarmaPredProb(glarmamod)
glarmaPIT(glarmamod)
}

\keyword{ts}
