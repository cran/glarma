\name{Asthma}
\docType{data}
\alias{Asthma}

\title{Daily Presentations of Asthma at Campbelltown Hospital}

\description{
The data arose from a single hospital (at Campbelltown) as part of a
larger (ongoing) study into the relationship between atmospheric
pollution and the number of asthma cases presenting themselves to
various emergency departments in local hospitals in the South West
region of Sydney, Australia.
}

\usage{data(Asthma)}


\format{
  A data frame containing the following columns:
  \tabular{rll}{

  [, 1]\tab Count \tab Daily counts of asthma at Campbelltown
  Hospital.\cr

  [, 2]\tab Intercept \tab A vector of ones, providing the intercept in
  the model.\cr

  [, 3]\tab Sunday \tab Takes value one for Sundays, otherwise
  zero.\cr

  [, 4]\tab Monday \tab Takes value one for Mondays, otherwise
  zero.\cr

  [, 5]\tab CosAnnual \tab cos((2*pi*t)/365), annual cosine term.\cr

  [, 6]\tab SinAnnual \tab sin((2*pi*t)/365), annual sine term.\cr

  [, 7]\tab H7 \tab Scaled lagged and smoothed humidity variable. \cr

  [, 8]\tab NO2max \tab Maximum daily nitrogen dioxide.\cr

  [, 9:16] \tab T1.1990 - T2.1993 \tab Smooth shapes to capture school
  terms in each year.\cr

 }
}


\source{
Davis, Richard A and Dunsmuir, William TM and Streett, Sarah B (2003)
Observation-driven models for Poisson counts. \emph{Biometrika},
\bold{90}, 777--790.
}
\examples{
### Example with asthma data
data(Asthma)
y <- Asthma[,1]
X <- as.matrix(Asthma[,2:16])

## Model in Davis, Dunsmuir and Streett (2003)

## MA(7) specification - see Davis, Dunsmuir and Streett (2003)

## Pearson Residuals, Fisher Scoring
glarmamod <- glarma(y, X, thetaLags = 7, type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)
glarmamod
summary(glarmamod)

likTests(glarmamod)
plot.glarma(glarmamod)

\dontrun{
## Example is specified as \dontrun because it takes too long
## for package inclusion on CRAN

## Pearson Residuals, Newton Raphson, Negative Binomial
## Initial value of the shape parameter take to be zero
glarmamod <- glarma(y, X, thetaLags = 7, type = "NegBin", method = "NR",
                    residuals = "Pearson", alphaInit = 0,
                    maxit = 100, grad = 1e-6)
glarmamod
summary(glarmamod)

likTests(glarmamod)
plot.glarma(glarmamod)
}

}
\keyword{datasets}