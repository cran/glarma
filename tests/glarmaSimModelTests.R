require(glarma)
### Test glarmaSimModel
data("Polio")
y <- Polio[, 2]
X <- as.matrix(Polio[, 3:8])
glarmamod <- glarma(y, X, thetaLags = c(1,2,5), type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)
summary(glarmamod)
PolioModel <- glarmaSimModel(X, beta = coef.glarma(glarmamod, type = "beta"),
                             phiLags = glarmamod$phiLags,
                             phi = rep(0.1, length(glarmamod$phiLags)),
                             thetaLags = glarmamod$thetaLags,
                             theta = rep(0.1, length(glarmamod$thetaLags)),
                             type = glarmamod$type,
                             residType = glarmamod$residType)
str(PolioModel)
coef.glarma(glarmamod, type = "ARMA")
PolioModel <- extractGlarmaSimModel(glarmamod)
str(PolioModel)
