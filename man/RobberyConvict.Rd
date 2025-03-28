\name{RobberyConvict}
\alias{RobberyConvict}
\docType{data}
\title{
  Court Convictions for Armed Robbery in New South Wales
}

\description{
  Monthly counts of charges laid and convictions made in Local Courts
  and Higher Court in armed robbery in New South Wales from 1995--2007.
}

\usage{data(RobberyConvict)}

\format{
  A data frame containing the following columns:
  \tabular{rll}{

  [, 1]\tab Date \tab Date in month/year format.\cr

  [, 2]\tab Incpt \tab A vector of ones, providing the intercept in the
  model.\cr

  [, 3]\tab Trend \tab Scaled time trend.\cr

  [, 4]\tab Step.2001 \tab Unit step change from 2001 onwards.\cr

  [, 5]\tab Trend.2001 \tab Change in trend term from 2001 onwards.\cr

  [, 6]\tab HC.N \tab Monthly number of cases for robbery (Higher
  Court).\cr

  [, 7]\tab HC.Y \tab Monthly number of convictions for robbery
  (Higher court).\cr

  [, 8]\tab HC.P \tab Proportion of convictions to charges for robbery
  (Higher court).\cr

  [, 9]\tab LC.N \tab Monthly number of cases for robbery (Lower
  court).\cr

  [, 10]\tab LC.Y \tab Monthly number of convictions for robbery
  (Lower court).\cr

  [, 11]\tab LC.P \tab Proportion of convictions to charges for
  robbery (Lower court).\cr
 }
}


\source{
Dunsmuir, William TM, Tran, Cuong, and Weatherburn, Don (2008)
\emph{Assessing the Impact of Mandatory DNA Testing of Prison Inmates
in NSW on Clearance, Charge and Conviction Rates for Selected Crime
Categories.}
}
\examples{
### Example with Robbery Convictions
data(RobberyConvict)
datalen <- dim(RobberyConvict)[1]
monthmat <- matrix(0, nrow = datalen, ncol = 12)
dimnames(monthmat) <- list(NULL, c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
months <- unique(months(strptime(RobberyConvict$Date, format = "\%m/\%d/\%Y"),
                        abbreviate=TRUE))


for (j in 1:12) {
  monthmat[months(strptime(RobberyConvict$Date,  "\%m/\%d/\%Y"),
                  abbreviate = TRUE) == months[j], j] <- 1
}

RobberyConvict <- cbind(rep(1, datalen), RobberyConvict, monthmat)
rm(monthmat)

## LOWER COURT ROBBERY
y1 <- RobberyConvict$LC.Y
n1 <- RobberyConvict$LC.N

Y <- cbind(y1, n1-y1)

glm.LCRobbery <- glm(Y ~ Step.2001 +
                        I(Feb + Mar + Apr + May + Jun + Jul) +
                        I(Aug + Sep + Oct + Nov + Dec),
                     data = RobberyConvict, family = binomial(link = logit),
                     na.action = na.omit, x = TRUE)

summary(glm.LCRobbery, corr = FALSE)

X <- glm.LCRobbery$x


## Newton Raphson
glarmamod <- glarma(Y, X, phiLags = c(1), type = "Bin", method = "NR",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)
glarmamod
summary(glarmamod)

## LRT, Wald tests.
likTests(glarmamod)

## Residuals and Fit Plots
plot.glarma(glarmamod)


## HIGHER COURT ROBBERY
y1 <- RobberyConvict$HC.Y
n1 <- RobberyConvict$HC.N

Y <- cbind(y1, n1-y1)

glm.HCRobbery <- glm(Y ~ Trend + Trend.2001 +
                       I(Feb + Mar + Apr + May + Jun) + Dec,
                     data = RobberyConvict, family = binomial(link = logit),
                     na.action = na.omit, x = TRUE)

summary(glm.HCRobbery,corr = FALSE)

X <- glm.HCRobbery$x


## Newton Raphson
glarmamod <- glarma(Y, X, phiLags = c(1, 2, 3), type = "Bin", method = "NR",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)
glarmamod
summary(glarmamod)


## LRT, Wald tests.
likTests(glarmamod)

## Residuals and Fit Plots
plot.glarma(glarmamod)
}

\keyword{datasets}
