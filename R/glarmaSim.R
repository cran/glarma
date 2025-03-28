### Program to simulate a single sample path of length n
### based on a glarma model object of class glarmaModel
###

glarmaSim <- function(object){
  if(!is(object, "glarmaSimModel"))
    stop("Model not of glarmaSimModel class")
  n <- NROW(X)
  ## Regression and offset part of glarma simulation model
  beta <- object$beta
  X <- object$X
  offset <- object$offset
  if(is.null(offset)) eta <- X %*% beta else eta <- X %*% beta + offset
  ## for Binomial need number of trials
  m <- object$mBin
  ## glarma recursion part
  phiLags <- object$phiLags
  phi <- object$phi
  thetaLags <- object$thetaLags
  theta <- object$theta
  p <- length(phiLags)
  q <- length(thetaLags)

  ## initialize arrays of sample size plus initializing conditions
  ## assume zero initilizing for glarma part e_t and Z_t
  if ((p + q) > 0) {
    mpq <- max(phiLags[p], thetaLags[q])
  }

  nmpq <- n + mpq
  e <- array(0, nmpq)
  Z <- array(0, nmpq)
  W <- array(0, nmpq)
  mu <- array(0, nmpq)
  sig2 <- array(0, nmpq)
  Y <- array(0, nmpq)
  for (time in 1:n){
    tmpq <- time + mpq
    if (p > 0) {
      for (i in 1:p) {
        Z[tmpq] <- Z[tmpq] + phi[i] * (Z + e)[tmpq - phiLags[i]]
      }}

    if (q > 0) {
      for (i in 1:q) {
        Z[tmpq] <- Z[tmpq] + theta[i] * e[tmpq - thetaLags[i]]
      }
    }
    W[tmpq] <- eta[time] + Z[tmpq]
    if(object$type == "Bin") {
      pi[tmpq] <- 1/(1 + exp(-W[tmpq]))
      mu[tmpq] <- m[time]*pi[tmpq]
      sig2[tmpq] <- mu[tmpq]*(1 - pi[tmpq])
      Y[tmpq] <- rbinom(1, m[time], 1/(1 + exp(-W[tmpq])))
      if(object$residType == "Pearson"){
        e[tmpq] <- (Y[tmpq] - mu[tmpq])/sig2[tmpq]^0.5
      }
      if(object$residType == "Score"){
        e[tmpq] <- (Y[tmpq] - mu[tmpq])/sig2[tmpq]
      }
      if(object$residType == "Identity"){
        e[tmpq] <- (Y[tmpq] - mu[tmpq])
      }
    }
    if(object$type == "Poi") {
      mu[tmpq] <- exp(W[tmpq])
      Y[tmpq] <- rpois(1, exp(W[tmpq]))
      if(object$residType == "Pearson"){
        e[tmpq] <- (Y[tmpq] - mu[tmpq])/mu[tmpq]^0.5
      }
      if(object$residType == "Score"){
        e[tmpq] <- (Y[tmpq] - mu[tmpq])/mu[tmpq]
      }
      if(object$residType == "Identity"){
        e[tmpq] <- (Y[tmpq] - mu[tmpq])
      }
    }
    if(object$type == "NegBin") {
      alpha <- object$alpha
      mu[tmpq] <- exp(W[tmpq])
      sig2[tmpq] <-  mu[tmpq] + mu[tmpq]^2/(alpha^2)
      Y[tmpq] <- rnbinom(1, mu = mu[tmpq], size =  alpha)
      if(object$residType == "Pearson"){
        e[tmpq] <- (Y[tmpq] - mu[tmpq])/sig2[tmpq]^0.5
      }
      if(object$residType == "Score"){
        e[tmpq] <- (Y[tmpq] - mu[tmpq])/sig2[tmpq]
      }
      if(object$residType == "Identity"){
        e[tmpq] <- (Y[tmpq] - mu[tmpq])
      }
    }

  }
  out <- list(Y = Y[mpq + 1:n], eta = eta,  W = W[mpq + 1:n],
              e = e[mpq + 1:n], mu = mu[mpq + 1:n], model = object)
  class(out) <- "glarmaSimulation"
  return(out)
}


