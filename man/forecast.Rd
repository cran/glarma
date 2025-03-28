\name{forecast}
\alias{forecast}
\alias{forecast.glarma}

\title{
  Forecasting GLARMA time series
}
\description{
\code{forecast} is a generic function for forecasting time series or
time series models. The function invokes particular \emph{methods} which
depend on the class of the first argument.

Currently the only method provided by the package is for objects of
class \code{"glarma"}.

}
\usage{

forecast(object, ...)
\method{forecast}{glarma}(object, n.ahead = 1, newdata = 0,
         newoffset = 0, newm = 1, ...)
}

\arguments{
  \item{object}{An object of class "glarma" obtained from a call to
    \code{glarma} }

  \item{n.ahead}{The number of periods ahead to be forecast.}

  \item{newdata}{The model matrix \eqn{X} comprising the values of the
    predictors for the times for which the series is to be
    predicted. Number of rows must be equal to \code{n.ahead} }

  \item{newoffset}{A vector containing the values of the offset for the
    times for which the series is to be predicted.  Length must be equal
    to \code{n.ahead} if the model includes an offset}

  \item{newm}{A vector containing the number of trials when forecasting
    binomial or binary time series. Defaults to the binary case. Length
    must be equal to \code{n.ahead}}

  \item{\dots}{Further arguments for the call, currently unused.  }

}
\details{
  Only one forecasting method is currently provided, for objects of
  class "glarma". This produces an object of class "glarmaForecast".

  When forecasting one step ahead, the values in the matrix
  \code{newdata} (and \code{offset} if there is an offset) in the GLARMA
  model are used along with the regression coefficients in the model to
  obtain the predicted value of \eqn{\eta}{eta}, the regression
  component of the state variable \eqn{W}. The predicted value of the
  ARMA component of the state variable is then added to this value to
  give the predicted value of \eqn{W}.

  When further predictions are required, since no data is available to
  calculate the predicted value of the state variable, an observation is
  generated from the predicted distribution and the methodology for one
  step ahead is then used on this generated data. This process is
  repeated until predictions are obtained for the required number of
  time periods (specified by \code{n.ahead}). Note that the value of
  \code{n.ahead} must equal the row dimension of \code{newdata} and if
  they are specified, of \code{newoffset} and \code{newm}.

  For completeness a randomly generated value of the time series is
  produced even for one step-ahead prediction.

  Note that the forecasted time series returned as the component
  \code{fitted} is then a randomly generated sample path for the
  predicted time series. If a sample of such paths is produced by
  repeated calls to \code{forecast} then sample predicted distributions
  can be obtained for the forecast series.

  In the case of binary or binomial time series in addition to values of
  the predictors in the regression component of the state variable and
  the values of any offset, the numbers of trials for the binomially
  distributed future observations are required. This information should
  be provided in the argument \code{newm}. If not, the number of trials
  defaults to 1, which is the case of binary responses.

}
\value{
  \code{forecast} currently has no default method.

  When \code{object} is of class \code{"glarma"}, \code{forecast}
  returns an object of class \code{"glarmaForecast"} with components:

  \item{eta}{the forecast values of the regression component of the
    state variable}

  \item{W}{the forecast values of the state variable \eqn{W}}

  \item{mu}{the conditional mean \eqn{\mu_t}{mu_t}}

  \item{Y}{the simulated series based on the fitted model}

  \item{n.ahead}{the number of steps ahead for which the forecasts were
    requested in the call to \code{forecast}}

  \item{newdata}{the model matrix \eqn{X} comprising the values of the
    predictors for the times for which the series is to be predicted}

  \item{newoffset}{the vector containing the values of the offset for
    the times for which the series is to be predicted}

  \item{newm}{the vector giving the number of trials when forecasting
    binomial or binary time series}

  \item{model}{the \code{"glarma"} object from the call to
    \code{forecast}}

}

\author{
  "William T.M. Dunsmuir" <w.dunsmuir@unsw.edu.au> and
  "David J Scott" <d.scott@auckland.ac.nz>
}

\examples{
require(zoo)
### Model number of deaths
data(DriverDeaths)
y <- DriverDeaths[, "Deaths"]
X <- as.matrix(DriverDeaths[, 2:5])
Population <- DriverDeaths[, "Population"]

### Offset included
glarmamod <- glarma(y, X, offset = log(Population/100000),
                    phiLags = c(12),
                    thetaLags = c(1),
                    type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)
print(summary(glarmamod))

XT1 <- matrix(X[72,], nrow = 1)
offsetT1 <- log(Population/100000)[72]

mu <- forecast(glarmamod, 1, XT1, offsetT1)$mu
print(mu)


### Save some values
allX <- X
allFits <- fitted(glarmamod)
ally <- y

### Look at a succession of forecasts
### Using actual values in forecasts
forecasts <- numeric(72)
for (i in (62:71)){
    y <- DriverDeaths[1:i, "Deaths"]
    X <- as.matrix(DriverDeaths[1:i, 2:5])
    Population <- DriverDeaths[1:i, "Population"]

    ## Offset included
    glarmamod <- glarma(y, X, offset = log(Population/100000),
                        phiLags = c(12),
                        thetaLags = c(1),
                        type = "Poi", method = "FS",
                        residuals = "Pearson", maxit = 100, grad = 1e-6)
    XT1 <- matrix(allX[i + 1, ], nrow = 1)
    offsetT1 <- log(DriverDeaths$Population[i + 1]/100000)
    mu <- forecast(glarmamod, 1, XT1, offsetT1)$mu
    if (i == 62){
        forecasts[1:62] <- fitted(glarmamod)
    }
    forecasts[i+1] <- mu
}
par(mfrow = c(1,1))
forecasts <- ts(forecasts[63:72], start = c(1985, 10), deltat = 1/12)
fitted <- ts(allFits, start = c(1980, 8), deltat = 1/12)
obs <- ts(DriverDeaths$Deaths, start = c(1980, 8), deltat = 1/12)
plot(obs, ylab = "Driver Deaths", lty = 2,
     main = "Single Vehicle Nighttime Driver Deaths in Utah")
points(obs)
lines(fitted, lwd = 2)
lines(forecasts, col = "red")
par(xpd = NA)
graph.param <-
    legend("top",
           legend = c("observations",expression(estimated~mu[t]),
                      expression(predicted~mu[t])),
           ncol = 3,
           cex = 0.7,
           bty = "n", plot = FALSE)
legend(graph.param$rect$left,
       graph.param$rect$top + graph.param$rect$h,
       legend = c("observations", expression(estimated~mu[t]),
                  expression(predicted~mu[t])),
       col = c("black","black","red"),
       lwd = c(1,2,1), lty = c(2,1,1),
       pch = c(1, NA_integer_, NA_integer_),
       ncol = 3,
       cex = 0.7,
       bty = "n",
       text.font = 4)
par(xpd = FALSE)

### Generate a sample of Y values 2 steps ahead and examine the distribution
data(DriverDeaths)
y <- DriverDeaths[, "Deaths"]
X <- as.matrix(DriverDeaths[, 2:5])
Population <- DriverDeaths[, "Population"]

### Fit the glarma model to the first 70 observations
glarmamod <- glarma(y[1:70], X[1:70, ],
                    offset = log(Population/100000)[1:70],
                    phiLags = c(12),
                    thetaLags = c(1),
                    type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)

nObs <- NROW(X)
n.ahead <- 2
### Specify the X matrix and offset for the times where predictions
### are required
XT1 <- as.matrix(X[(nObs - n.ahead + 1):nObs, ])
offsetT1 <- log(Population/100000)[(nObs - n.ahead + 1):nObs]
nSims <- 500
forecastY <- matrix(ncol = n.ahead, nrow = nSims)
forecastMu <- matrix(ncol = n.ahead, nrow = nSims)

### Generate sample predicted values
for(i in 1:nSims){
    temp <-  forecast(glarmamod, n.ahead, XT1, offsetT1)
    forecastY[i, ] <- temp$Y
    forecastMu[i, ] <- temp$mu
}
### Examine distribution of sample of Y values n.ahead
table(forecastY[, 2])
par(mfrow = c(2,1))
barplot(table(forecastY[, 2]),
        main = "Barplot of Sample Y Values 2 Steps Ahead")
hist(forecastY[, 2], xlab = "Sample Y values",
     main = paste0("Histogram of Sample Y Values 2 Steps Ahead",
                   "\nwith 0.025 and 0.975 Quantiles"))
abline(v = quantile(forecastY[, 2], c(0.025, 0.975)), col = "red")


}

\keyword{ts}
\keyword{methods }
