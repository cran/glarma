\name{glarmaSim}
\alias{glarmaSim}

\title{
  Simulate a glarma process.
}
\description{
  Simulate a single instance of a glarma process specified by an object
  of class \code{glarmaSim}.
}
\usage{
glarmaSim(object)
}
\arguments{
  \item{object}{An object of class "glarmaSim" either constructed to be
  of that class or derived from a fitted glarma model using the function
  \code{extractGlarmaSimModel}. }
}

\value{
  An object of class \code{"glarmaSimulation"}.
}

\author{
  "William T.M. Dunsmuir" <w.dunsmuir@unsw.edu.au> and
  "David J Scott" <d.scott@auckland.ac.nz>
}

\seealso{
  \code{\link{glarmaSim}}
}
\examples{
### Polio data
data("Polio")
y <- Polio[, 2]
X <- as.matrix(Polio[, 3:8])
glarmamod <- glarma(y, X, thetaLags = c(1,2,5), type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)
PolioModel <- extractGlarmaSimModel(glarmamod)
str(PolioModel)
par(mfrow = c(3,1))
sim <- glarmaSim(PolioModel)
ts.plot(sim$W)
ts.plot(sim$mu)
ts.plot(sim$Y)

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
str(glarmamod)
BoatRaceModel <- glarmaSimModel(X, beta = coef.glarma(glarmamod, type = "beta"),
                                phiLags = glarmamod$phiLags,
                                phi = rep(0.1, length(glarmamod$phiLags)),
                                thetaLags = glarmamod$thetaLags,
                                theta = c(0.34,0.56),
                                type = glarmamod$type, m = n1,
                                residType = glarmamod$residType)
str(BoatRaceModel)
BoatRaceModel <- extractGlarmaSimModel(glarmamod)
str(BoatRaceModel)
par(mfrow = c(3,1))
sim <- glarmaSim(BoatRaceModel)
ts.plot(sim$W)
ts.plot(sim$mu)
ts.plot(sim$Y)


### Example with asthma data, negative binomial
data(Asthma)
y <- Asthma[, 1]
X <- as.matrix(Asthma[, 2:16])

## Pearson Residuals, Newton Raphson, Negative Binomial
## Initial value of the shape parameter take to be zero
glarmamod <- glarma(y, X, thetaLags = 7, type = "NegBin", method = "NR",
                    residuals = "Pearson", alphaInit = 0,
                    maxit = 100, grad = 1e-6)
summary(glarmamod)
str(glarmamod)
AsthmaModel <- glarmaSimModel(X, beta = coef.glarma(glarmamod, type = "beta"),
                              phiLags = glarmamod$phiLags,
                              phi = rep(0.1, length(glarmamod$phiLags)),
                              thetaLags = glarmamod$thetaLags,
                              theta = c(0.044),
                              type = glarmamod$type, alpha = 37.19,
                              residType = glarmamod$residType)
str(AsthmaModel)
AsthmaModel <- extractGlarmaSimModel(glarmamod)
str(AsthmaModel)
par(mfrow = c(3,1))
sim <- glarmaSim(AsthmaModel)
ts.plot(sim$W)
ts.plot(sim$mu)
ts.plot(sim$Y)
}
\keyword{ ts }
\keyword{ method }
