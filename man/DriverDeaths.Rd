\name{DriverDeaths}
\docType{data}
\alias{DriverDeaths}

\title{Single Vehicle Nighttime Driver Deaths in Utah}

\description{
  This data set gives the number of single vehicle nighttime
  driver deaths in the state of Utah by month over the period August
  1980 to July 1986, along with observations on a number of possible
  predictors. The aim of the study from which it was taken was to
  investigate the effect of the lowering of the legal blood alcohol
  concentration (BAC) while driving, from 0.1 to 0.08 units, and the
  simultaneous introduction of administrative license revocation. The
  time period for the observations is centred on the month of the
  intervention, August 1983.
}

\usage{data(DriverDeaths)}


\format{
  A data frame containing the following columns:
  \tabular{rll}{

  [, 1]\tab Deaths \tab Number of single vehicle nighttime driver deaths
  monthly.\cr

  [, 2]\tab Intercept \tab A vector of ones, providing the intercept in
  the model.\cr

  [, 3]\tab ReducedBAC \tab Indicator of before or after lowering of
  legal blood alcohol level.0 for months prior to August 1983, 1 for
  months on or after August 1983.\cr

  [, 4]\tab FriSat \tab Number of Friday and Saturday nights in the month.\cr

  [, 5]\tab lnOMVDRate \tab Log of the number of other motor vehicle
  deaths per 100,000 of population.\cr

  [, 6]\tab Population \tab Adult population of the State of Utah.\cr
 }
}


\source{
  Debra H. Bernat, William T.M. Dunsmuir, and Alexander C. Wagenaar (2004)
  Effects of lowering the legal BAC to 0.08 on single-vehicle-nighttime
  fatal traffic crashes in 19 jurisdictions. \emph{Accident Analysis
  & Prevention}, \bold{36}, 1089--1097.
}
\examples{
### Model number of deaths
data(DriverDeaths)
y <- DriverDeaths[, "Deaths"]
X <- as.matrix(DriverDeaths[, 2:5])
Population <- DriverDeaths[, "Population"]

### Offset included
glarmamodOffset <- glarma(y, X, offset = log(Population/100000),
                          phiLags = c(12),
                          type = "Poi", method = "FS",
                          residuals = "Pearson", maxit = 100, grad = 1e-6)
print(summary(glarmamodOffset))
par(mfrow =c(3,2))
plot(glarmamodOffset)


### No offset included
glarmamodNoOffset <- glarma(y, X, phiLags = c(12),
                            type = "Poi", method = "FS",
                            residuals = "Pearson", maxit = 100, grad = 1e-6)
print(summary(glarmamodNoOffset))
par(mfrow=c(3,2))
plot(glarmamodNoOffset)
}

\keyword{datasets}