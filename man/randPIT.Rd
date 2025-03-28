\name{normRandPIT}
\alias{normRandPIT}

\title{
  Random normal probability integral transformation
}
\description{
  Function to create the normalized conditional (randomized) quantile
  residuals.
}
\usage{
normRandPIT(object)
}
\arguments{
  \item{object}{an object of class "glarma"}
}
\details{
  The function \code{\link{glarmaPredProb}} produces the non-randomized
  probability integral transformation (PIT). It returns estimates of the
  cumulative predictive probabilities as upper and lower bounds of a
  collection of intervals. If the model is correct, a histogram drawn
  using these estimated probabilities should resemble a histogram
  obtained from a sample from the uniform distribution. This function
  aims to produce observations which instead resemble a sample from a
  normal distribution. Such a sample can then be examined by the usual
  tools for checking normality, such as histograms, Q-Q normal plots and
  for checking independence, autocorrelation and partial autocorrelation
  plots, and associated portmanteau statistics.

  For each of the intervals produced by \code{\link{glarmaPredProb}}, a
  random uniform observation is generated, which is then converted to a
  normal observation by applying the inverse standard normal
  distribution function (that is \code{\link{qnorm}}). The vector of
  these values is returned by the function in the list element
  \code{rt}. In addition non-random observations which should appear
  similar to a sample from a normal distribution are obtained by
  applying \code{qnorm} to the mid-points of the predictive distribution
  intervals. The vector of these values is returned by the function in
  the list element \code{rtMid}.

}
\value{
  A list consisting of two elements:

  \item{rt}{the normalized conditional (randomized) quantile residuals}

  \item{rtMid}{the midpoints of the predictive probability intervals}

}
\references{
  Berkowitz, J. (2001) Testing density forecasts, with applications to
  risk management. \emph{Journal of Business & Economic Statistics},
  \bold{19}, 465--474.

  Dunn, Peter K. and Smyth, Gordon K. (1996) Randomized quantile
  residuals. \emph{Journal of Computational and Graphical Statistics},
  \bold{5}, 236--244.
}
\author{
  "William T.M. Dunsmuir" <w.dunsmuir@unsw.edu.au> and
  "David J Scott" <d.scott@auckland.ac.nz>
}
\seealso{
  See also as \code{\link{glarmaPredProb}}.
}
\examples{
data(DriverDeaths)
y <- DriverDeaths[, "Deaths"]
X <- as.matrix(DriverDeaths[, 2:5])
Population <- DriverDeaths[, "Population"]

### Offset included
glarmamodOffset <- glarma(y, X, offset = log(Population/100000),
                          phiLags = c(12),
                          type = "Poi", method = "FS",
                          residuals = "Pearson", maxit = 100, grad = 1e-6)
rt <- normRandPIT(glarmamodOffset)$rt
par(mfrow = c(2,2))
hist(rt, main = "Histogram of Randomized Residuals",
     xlab = expression(r[t]))
box()
qqnorm(rt, main = "Q-Q Plot of  Randomized Residuals" )
abline(0, 1, lty = 2)
acf(rt, main = "ACF of Randomized Residuals")
pacf(rt, main = "PACF of Randomized Residuals")
}
\keyword{ts }
