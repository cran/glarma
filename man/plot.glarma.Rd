\name{plot.glarma}
\alias{plot.glarma}

\title{Plot Diagnostics for a \code{glarma} Object}

\description{
  Six plots (selectable by \code{which}) are currently available: a
  time series plot with observed values of the dependent variable, fixed
  effects fit, and \code{GLARMA} fit; an ACF plot of residuals; a plot
  of residuals against time; a Normal Q-Q plot; the PIT histogram; and
  the uniform Q-Q plot for the PIT. By default, all plots are provided.
}

\usage{
\method{plot}{glarma}(x, which = c(1L:6L), fits = 1L:3L,
     ask = prod(par("mfcol")) < length(which) && dev.interactive(),
     lwdObs = 1, lwdFixed = 1, lwdGLARMA = 1,
     colObs = "black", colFixed = "blue", colGLARMA = "red",
     ltyObs = 2, ltyFixed = 1, ltyGLARMA = 1,
     pchObs = 1, legend = TRUE, residPlotType = "h", bins = 10,
     line = TRUE, colLine = "red", colHist = "royal blue",
     lwdLine = 2, colPIT1 = "red", colPIT2 = "black",
     ltyPIT1 = 1, ltyPIT2 = 2, typePIT = "l",...)
}

\arguments{

 \item{x}{An object of class \code{"glarma"}, obtained from a call to
   \code{\link{glarma}}.}

 \item{which}{Numeric; if a subset of the plots is required, specify a
   subset of the numbers 1:6. 1 is the time series plot with observed
   values of the dependent variable, fixed effects fit, and \code{GLARMA}
   fit. 2 is the ACF plot of residuals. 3 is a plot of residuals
   against time. 4 is the Normal Q-Q plot. 5 is the PIT histogram. 6
   is the uniform Q-Q plot for the PIT. By default, all plots are
   provided.}

 \item{fits}{Numeric; if a subset of fits on the time series plot is
   required, specify a subset of the numbers 1:3. 1 is the observed
   values of the dependent variable, 2 is the fixed effects fit, and 3 is
   \code{GLARMA} fit. By default, all fits are provided.}

 \item{ask}{Logical; if \code{TRUE}, the user is asked before each
   plot, see \code{\link{par}(ask = .)}.}

 \item{lwdObs}{Numeric; the line widths for lines of the observed values
   of the dependent variable appearing in the time series plot.}

 \item{lwdFixed}{Numeric; the line widths for lines of the fixed
   effects fit appearing in the time series plot.}

 \item{lwdGLARMA}{Numeric; the line widths for lines of \code{GLARMA}
   fit appearing in the time series plot.}

 \item{ltyObs}{An integer or character string; the line types for the
   line of the observed data of the dependent variable appearing in the
   time series plot, see \code{\link{par}(lty = .)}.}

 \item{ltyFixed}{An integer or character string; the line types for the
   line of the fixed effects fit appearing in the time series plot, see
   \code{\link{par}(lty = .)}.}

 \item{ltyGLARMA}{An integer or character string; the line types for the
   line of \code{GLARMA} fit appearing in the time series plot, see
   \code{\link{par}(lty = .)}.}

 \item{pchObs}{Numeric; the point type for the point of the observed
   data of the dependent variable appearing in the time series plot.}

 \item{colObs}{Numeric or character; the colour of lines or points of
   the observed data of the dependent variable appearing in the time
   series plot.}

 \item{colFixed}{Numeric or character; the colour of lines of the
   fixed effects fit appearing in the time series plot.}

 \item{colGLARMA}{Numeric or character; the colour of lines of
   \code{GLARMA} fit appearing in the time series plot.}

 \item{legend}{Logical; if \code{TRUE}, the legend for the fits in the
   time series plot would be shown. By default, it would be shown.}

 \item{residPlotType}{A 1-character string giving the type of plot
   desired. The following values are possible, for details, see
   \code{\link{plot}}: "\code{p}" for points, "\code{l}" for lines,
   "\code{b}" for both points and lines, "\code{c}" for empty points
   joined by lines, "\code{o}" for overplotted points and lines,
   "\code{s}" and "\code{S}" for stair steps and "\code{h}" for
   histogram-like vertical lines. Finally, "\code{n}" does not produce
   any points or lines.}

 \item{bins}{Numeric; the number of bins shown in the PIT
   histogram. By default, it is 10.}

 \item{line}{Logical; if \code{TRUE}, the line for displaying the
   standard uniform distribution will be shown for the purpose of
   comparison. The default is \code{TRUE}.}

 \item{colLine}{Numeric or character; the colour of the line for
   comparison in the PIT histogram.}

 \item{lwdLine}{Numeric; the line widths for the comparison line in
   the PIT histogram.}

 \item{colHist}{Numeric or character; the colour of the histogram for
   the PIT.}

 \item{colPIT1}{Numeric or character; the colour of the sample uniform
   Q-Q plot in the PIT.}

 \item{colPIT2}{Numeric or character; the colour of the theoretical
   uniform Q-Q plot in the PIT.}

 \item{ltyPIT1}{An integer or character string; the line types for the
   sample uniform Q-Q plot in the PIT, see \code{\link{par}(lty = .)}.}

 \item{ltyPIT2}{An integer or character string; the line types for the
   theoretical uniform Q-Q plot in the PIT, see \code{\link{par}(lty = .)}.}

 \item{typePIT}{A 1-character string; the type of plot for the sample
   uniform Q-Q plot in the PIT.}

 \item{...}{Further arguments passed to \code{\link{plot.default}} and
   \code{\link{plot.ts}}.}

}

\details{

  \code{plot.glarma} is an S3 generic function for objects of class
  \code{glarma}.

  The plots in this method display the fixed effects fit,
  GLARMA fit and various types of residuals for the GLARMA fit
  under the Poisson distribution, the binomial distribution
  or the negative binomial distribution. The observed values of the
  dependent variable shown in the time series plot are mainly used to
  compare with the two fits.

  The fixed effects fit is calculated from \eqn{\eta}{eta}, the
  multiplication of the data matrix \code{X} and \eqn{\beta}{beta}
  coefficients in GLARMA model. In contrast, the GLARMA fit is
  calculated from \eqn{W}{W}, the product of the data matrix \code{X}
  and \eqn{\delta}{delta} in the GLARMA model, which is the combination
  of both the \eqn{\beta}{beta} and ARMA coefficients.

  There are some major differences for computing the fixed effects fit
  from \eqn{\eta}{eta} and the GLARMA fit from \eqn{W} under
  different distributions.

  Under the Poisson distribution and negative binomial distribution,
   \deqn{\mathsf{fit}_{\mathsf{fixed}} = \exp{\eta}}{%
     fixed effects fit = exp(eta)}
   and
   \deqn{\mathsf{fit}_{\mathsf{glarma}} = \exp{W}.}{%
     glarma fit = exp(W).}

   Under the binomial distribution,
   \deqn{\mathsf{fit}_{\mathsf{fixed}} = \frac{1}{(1+e^{-\eta})}}{%
     fixed effects fit = 1 / (1 + exp(-eta))}
   and
   \deqn{\mathsf{fit}_{\mathsf{glarma}} = \frac{1}{(1+e^{-W})}.}{%
         glarma fit = 1 / (1 + exp(-W)).}

   The residuals are calculated from the observed data and GLARMA
   fit. The exact computation for the residuals depends on the type
   of residuals used. The details are given in
   \code{\link{glarma}}. The ACF plot, the residuals against time and
   the normal Q-Q plot are all based on these residuals. Further details
   about those three plots are passed to \code{\link{acf}} and
   \code{\link{qqnorm}}.

   The number of plots to be shown in the window depends on the value of
   the graphical parameter \code{mfrow} (or \code{mfcol}). If the
   displayed window is set to be large enough to show all six plots,
   they will be shown at one time. Otherwise, the required number of
   plots will appear each time in the displayed window, and the user
   will need to enter \code{return} to see subsequent plots.

   For the time series plot in the function, the fit displayed is
   specified by the argument \code{fits}. The legend will be shown if
   \code{legend} is \code{TRUE}. It will appear under the title
   of the time series plot. Also the legend and the title will alter
   automatically according to the fits shown in the plot.

}

\author{
  "Cenanning Li" <cli113@aucklanduni.ac.nz>
}


\seealso{
\code{\link{plot.ts}}, \code{\link{qqnorm}}, \code{\link{acf}},
\code{\link{plot.default}}.
}

\examples{
### A example from Davis, Dunsmuir Wang (1999)
## MA(1,2,5), Pearson Residuals, Fisher Scoring
data(Polio)
y <- Polio[, 2]
X <- as.matrix(Polio[, 3:8])
glarmamod <- glarma(y, X, thetaLags = c(1, 2, 5), type = "Poi",method = "FS",
                    residuals = "Pearson", maxit = 100 , grad = 1e-6)

## The default plots are shown
plot(glarmamod)

## The plots used only to compared GLARMA fit and the observed data
plot(glarmamod, which = 1L, fits = c(1, 3))
}

\keyword{Diagnostic}
