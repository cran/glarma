\name{plotPIT}
\alias{histPIT}
\alias{qqPIT}

\title{
PIT Plots for a \code{glarma} Object
}

\description{

Two plots for the non-randomized PIT are currently available for
checking the distributional assumption of the fitted GLARMA model: the
PIT histogram, and the uniform Q-Q plot for PIT.

}

\usage{
histPIT(object, bins = 10, line = TRUE, colLine = "red",
        colHist = "royal blue", lwdLine = 2, main = NULL, ...)
qqPIT(object, bins = 10, col1 = "red", col2 = "black",
      lty1 = 1, lty2 = 2, type = "l", main = NULL, ...)
}


\arguments{
  \item{object}{An object of class \code{"glarma"}, obtained from a
    call to \code{\link{glarma}}.}

  \item{bins}{Numeric; the number of bins shown in the PIT
    histogram or the PIT Q-Q plot. By default, it is 10.}

  \item{line}{Logical; if \code{TRUE}, the line for displaying the
    standard uniform distribution will be shown for the purpose of
    comparison. The default is \code{TRUE}.}

  \item{colLine}{Numeric or character; the colour of the line for
    comparison in PIT histogram.}

  \item{lwdLine}{Numeric; the line widths for the comparison line in
    PIT histogram.}

  \item{colHist}{Numeric or character; the colour of the histogram for
    PIT.}

  \item{col1}{Numeric or character; the colour of the sample uniform
    Q-Q plot in PIT.}

  \item{col2}{Numeric or character; the colour of the theoretical
    uniform Q-Q plot in PIT.}

  \item{lty1}{An integer or character string; the line types for the
    sample uniform Q-Q plot in PIT, see \code{\link{par}(lty = .)}.}

  \item{lty2}{An integer or character string; the line types for the
    theoretical uniform Q-Q plot in PIT, see \code{\link{par}(lty =
    .)}.}

  \item{type}{A 1-character string; the type of plot for the sample
    uniform Q-Q plot in PIT.}

  \item{main}{A character string giving a title. For each plot the
    default provides a useful title.}

  \item{...}{Further arguments passed to \code{\link{plot.default}}
    and \code{\link{plot.ts}}.}  }

\details{
The histogram and the Q-Q plot are used to compare the fitted profile
with \bold{U}(0, 1). If they match relatively well, it means the
distributional assumption is satisfied.
}


\references{
Czado, Claudia and Gneiting, Tilmann and Held, Leonhard (2009)
Predictive model assessment for count data. \emph{Biometrics},
\bold{65}, 1254--1261.

Jung, Robert.C and Tremayne, A.R (2011) Useful models for time series
of counts or simply wrong ones? \emph{AStA Advances in Statistical
Analysis}, \bold{95}, 59--91.
}

\author{
  "David J. Scott" <d.scott@auckland.ac.nz> and
  "Cenanning Li" <cli113@aucklanduni.ac.nz>
}

\examples{
## For examples see example(plot.glarma)
}

\keyword{Diagnostic}
