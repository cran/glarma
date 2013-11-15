\name{OxBoatRace}
\docType{data}
\alias{OxBoatRace}

\title{Oxford-Cambridge Boat Race}

\description{
Results of the boat race between Oxford and Cambridge from
1829--2011.
}

\usage{data(OxBoatRace)}

\format{
  A data frame containing the following columns:
  \tabular{rll}{

  [, 1]\tab Year \tab Year in which the race occurred. Some years are
  missing when the race was not run.\cr

  [, 2]\tab Intercept \tab A vector of ones, providing the intercept in
  the model.\cr

  [, 3]\tab Camwin \tab A binary response, zero for an Oxford win, one
  for a Cambridge win.\cr

  [, 4]\tab WinnerWeight \tab Weight of winning team's crew.\cr

  [, 5]\tab LoserWeight \tab Weight of losing team's crew.\cr

  [, 6]\tab Diff \tab Difference between winning team's weight and
  losing team's weight.\cr
 }
}


\source{
Klingenberg, Bernhard (2008) Regression models for binary time series
with gaps. \emph{Computational Statistics & Data Analysis},
\bold{52}, 4076--4090.
}
\examples{
### Example with Oxford-Cambridge Boat Race
data(OxBoatRace)

y1 <- OxBoatRace$Camwin
n1 <- rep(1, length(OxBoatRace$Year))
Y <- cbind(y1, n1 - y1)
X <- cbind(OxBoatRace$Intercept, OxBoatRace$Diff)
colnames(X) <- c("Intercept", "Weight Diff")

oxcamglm <- glm(Y ~ Diff + I(Diff^2),
                data = OxBoatRace,
                family = binomial(link = "logit"), x = TRUE)
summary(oxcamglm)

X <- oxcamglm$x

glarmamod <- glarma(Y, X, thetaLags = c(1, 2), type = "Bin", method = "NR",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)

summary(glarmamod)
likTests(glarmamod)

## Plot Probability of Cambridge win versus Cambridge Weight advantage:
beta <- coef(glarmamod, "beta")
par(mfrow = c(1, 1))
plot(OxBoatRace$Diff, 1 / (1 + exp(-(beta[1] + beta[2] * OxBoatRace$Diff +
                                       beta[3] * OxBoatRace$Diff^2))),
     ylab = "Prob", xlab = "Weight Diff")
title("Probability of Cambridge win \n versus Cambridge weight advantage")

## Residuals and fit plots
par(mfrow=c(3, 2))
plot.glarma(glarmamod)
}
\keyword{datasets}