require(glarma)
### Extract model from a glarma fit to use in simulation
### Example with asthma data
data(Asthma)
y <- Asthma[,1]
X <- as.matrix(Asthma[,2:16])
## Pearson Residuals, Fisher Scoring
glarmamod <- glarma(y, X, thetaLags = 7, type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)
mod <- extractGlarmaSimModel(glarmamod)
str(mod)

glarmamod <- glarma(y, X, thetaLags = 7, type = "NegBin", method = "NR",
                    residuals = "Pearson", alphaInit = 0,
                    maxit = 100, grad = 1e-6)
mod <- extractGlarmaSimModel(glarmamod)
str(mod)

data(DriverDeaths)
y <- DriverDeaths[, "Deaths"]
X <- as.matrix(DriverDeaths[, 2:5])
Population <- DriverDeaths[, "Population"]

### Offset included
glarmamodOffset <- glarma(y, X, offset = log(Population/100000),
                          phiLags = c(12),
                          type = "Poi", method = "FS",
                          residuals = "Pearson", maxit = 100, grad = 1e-6)
mod <- extractGlarmaSimModel(glarmamodOffset)
str(mod)

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
mod <- extractGlarmaSimModel(glarmamod)
str(mod)


data(RobberyConvict)
datalen <- dim(RobberyConvict)[1]
monthmat <- matrix(0, nrow = datalen, ncol = 12)
dimnames(monthmat) <- list(NULL, c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
months <- unique(months(strptime(RobberyConvict$Date, format = "%m/%d/%Y"),
                        abbreviate=TRUE))


for (j in 1:12) {
  monthmat[months(strptime(RobberyConvict$Date,  "%m/%d/%Y"),
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
mod <- extractGlarmaSimModel(glarmamod)
str(mod)
