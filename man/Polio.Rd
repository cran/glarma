\name{Polio}
\docType{data}
\alias{Polio}
\title{Cases of Poliomyelitis in the U.S.}
\description{
  This data set gives the monthly number of cases of poliomyelitis in
  the U.S. for the years 1970--1983 as reported by the Center for
  Disease Control. The polio data frame has 168 rows and 8 columns.
}

\usage{data(Polio)}
\format{
  A data frame containing the following columns:
  \tabular{rll}{
  [, 1]\tab Cases \tab monthly number of cases of poliomyelitis.\cr
  [, 2]\tab Intcpt \tab a vector of ones, providing the intercept in
   the model.\cr
  [, 3]\tab Trend \tab a linear trend.\cr
  [, 4]\tab CosAnnual \tab cosine harmonics at periods of 12.\cr
  [, 5]\tab SinAnnual \tab sine harmonics at periods of 12.\cr
  [, 6]\tab CosSemiAnnual \tab cosine harmonics at periods of 6.\cr
  [, 7]\tab SinSemiAnnual \tab sine harmonics at periods of 6.\cr
  }
 }

\source{
  Zeger, S.L (1988) A regression model for time series of
  counts. \emph{Biometrika}, \bold{75}, 621--629.
}

\keyword{datasets}