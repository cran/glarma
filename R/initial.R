initial <- function(y, X, type = "Poi", alpha = 1) {

    if (type == "Poi"){
        (GLM <- glm(y ~ -1 + X, family = poisson, x = TRUE))
        GLMInt <- glm(y ~ X, family = poisson, x = TRUE)
    }


    if (type == "NegBin"){
        (GLM <- glm.nb(y ~ -1 + X, init.theta = alpha, x = TRUE))
        GLMInt <- glm.nb(y ~ X, init.theta = alpha, x = TRUE)
    }

    if (type == "Bin"){
        (GLM <- glm(y ~ -1 + X, family = binomial(link = logit),
                    na.action = na.omit, x = TRUE))
        GLMInt <- glm(y ~ X, family = binomial(link = logit),
                   na.action = na.omit, x = TRUE)
    }

    list(beta = GLM$coefficients, y = y, X = X, alpha = GLM$theta, type = type,
         null.deviance = GLMInt$null.deviance, df.null = GLMInt$df.null)

}
