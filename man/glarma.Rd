\name{glarma}
\alias{glarma}
\alias{glarmaPoissonPearson}
\alias{glarmaPoissonScore}
\alias{glarmaBinomialIdentity}
\alias{glarmaBinomialPearson}
\alias{glarmaBinomialScore}
\alias{glarmaNegBinPearson}
\alias{glarmaNegBinScore}
\alias{print.glarma}

\title{Generalized Linear Autoregressive Moving Average Models with
  Various Distributions}


\description{

  The function \code{glarma} is used to fit generalized linear
  autoregressive moving average models with various distributions
  (Poisson, binomial, negative binomial) using either Pearson residuals
  or score residuals, and for the binomial distribution, identity
  residuals. It also estimates the parameters of the GLARMA model with
  various distributions by using Fisher scoring or Newton-Raphson
  iteration.

  For Poisson and negative binomial response distributions the log link
  is currently used. For binomial responses the logit link is currently
  used.
}

\usage{
glarma(y, X, offset = NULL, type = "Poi", method = "FS", residuals = "Pearson",
       phiLags,  thetaLags, phiInit, thetaInit, beta, alphaInit,
       alpha = 1, maxit = 30, grad = 2.22e-16)

glarmaPoissonPearson(y, X, offset = NULL, delta, phiLags, thetaLags,
                     method =  "FS")

glarmaPoissonScore(y, X, offset = NULL, delta, phiLags, thetaLags,
                   method =  "FS")

glarmaBinomialIdentity(y, X, offset = NULL, delta, phiLags, thetaLags,
                       method =  "FS")

glarmaBinomialPearson(y, X, offset = NULL, delta, phiLags, thetaLags,
                      method =  "FS")

glarmaBinomialScore(y, X, offset = NULL, delta, phiLags, thetaLags,
                    method =  "FS")

glarmaNegBinPearson(y, X, offset = NULL, delta, phiLags, thetaLags,
                    method =  "FS")

glarmaNegBinScore(y, X, offset = NULL, delta, phiLags, thetaLags,
                  method =  "FS")
}

\arguments{

  \item{y}{Numeric vector; the response variable. If the
    response variable is for the model with the binomial
    distribution, it should be a n by 2 matrix, one column is the
    number of successes and another is the number of failures.}

  \item{X}{Matrix; the explanatory variables. A vector of ones should be
    added to the data matrix as the first column for the \code{beta} of
    the intercept.}

  \item{offset}{Either \code{NULL} or a numeric vector of length equal
    to the number of cases. Used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.}

  \item{beta}{Numeric vector; initial values of the regression
    coefficients.}

  \item{phiLags}{Numeric vector; AR orders.}

  \item{phiInit}{Numeric vector; initial values for the corresponding
    AR orders.}

  \item{thetaLags}{Numeric vector; MA orders.}

  \item{thetaInit}{Numeric vector; initial values for the corresponding
    MA orders.}

  \item{delta}{Numeric vector; initial values of the parameters for the
    GLARMA estimation procedure. It is a combination of the parameters of
    \code{beta}, the AR terms and the MA terms.}

  \item{alpha}{Numeric; an optional initial shape parameter for
    \code{\link[MASS]{glm.nb}}.}

  \item{alphaInit}{Numeric; an initial shape parameter for
    \code{glarma} for negative binomial counts.}

  \item{type}{Character; the count distribution. Possible values are
    \code{"Poi"} (Poisson), \code{"Bin"} (binomial) and \code{"NegBin"}
    (negative binomial). The default is the
    Poisson distribution.}

  \item{method}{Character; method of iteration to be used. Possible
    values are \code{"FS"} (Fisher scoring), and \code{"NR"}
    (Newton-Raphson). The default is to use Fisher scoring to estimate
    the parameters of a GLARMA model.}

  \item{residuals}{Character; the type of residuals to be used. Possible
    values are \code{"Pearson"} and \code{"Score"}, and for the binomial
    distribution \code{"Identity"} is also allowed. The default is to
    use Pearson residuals.}

  \item{maxit}{Numeric; the maximum number of iterations allowed.}

  \item{grad}{Numeric; the tolerance for recognizing numbers, which
    are smaller than the specified tolerance, as zero.}

}

\details{

  Models for \code{glarma} are specified symbolically. A typical model
  has the form \code{y} (response), \code{X} (terms) where \code{y} is
  the count or factor reponse vector, \code{X} is a series of terms
  which specifies a linear predictor for the response. It should be noted
  that the first column of \code{X} should be a vector of 1s as the
  intercept in the model. Four initial parameters that need to be
  estimated are combined into \eqn{\delta = (\beta, \phi, \theta,
  \alpha)}{delta = (beta, phi, theta, alpha)}, where \eqn{\alpha}{alpha}
  is an optional parameter to accomodate the negative binomial
  model. Note that in the function \code{\link[MASS]{glm.nb}} from the
  package \pkg{MASS}, this parameter is called \code{theta}.

  For Poisson and negative binomial response distributions the log link
  is currently used. For binomial responses the logit link is currently
  used.

  The generalized linear autoregressive moving average models are
  computed as follows.

  The linear predictor for the response is
  \deqn{\log{\mu_t} = W_t = X_t^T\beta + \mbox{offset} + Z_t.}{log(mu_t) = %
    W_t = transpose(X_t) * beta + offset + Z_t.}

  The infinite moving average from the linear predictor is
  \deqn{Z_t = \sum_{i=1}^\infty \gamma_i e_{t-i}.}{Z_t = sum(gamma_i *%
    residuals_(t-i)).}

  This infinite moving average, is computed using the autoregressive
  moving average recursions

  \deqn{Z_t = \phi_1 (Z_{t-1} + e_{t-1}) + ... +%
              \phi_p (Z_{t-p} + e_{t-p}) + \theta_1 e_{t-1}%
	      + ... + \theta_q e_{t-q}}{%
	      Z_t = phi_1 * (Z_(t-1) + e_(t-1)) + ... + %
	      phi_p * (Z_{t-p} + e_(t-p)) + theta_1 * e_{1} + ... + %
	      theta_q * e_{t-q}}

    where \eqn{p} and \eqn{q} are the orders of \eqn{\phi}{phi}
    and \eqn{\theta}{theta} respectively and the non-zero lags of the vectors
    \code{phi} and \code{theta} may be specified by the user via the
    arguments \code{phiLag} and \code{thetaLag}.

    There are two types of residuals which may be used in each
    recursion, Pearson residuals or score residuals, and in addition,
    for the binomial distribution, identity residuals may be used. The
    infinite moving average, \eqn{Z_t}{Z_t}, depends on the type of
    residuals used, as do the final parameters obtained from the
    filter. Standardisation of past observed counts is necessary to
    avoid instability, therefore the user should choose the appropriate
    type of residuals depending on the situation.

    The method of estimation for parameters implemented in the function
    aims to maximise the log likelihood by an iterative method commencing
    from suitably chosen initial values for the parameters. Starting from
    initial values \eqn{\hat\delta^{(0)}}{delta hat^(0)} for the vector of
    parameters updates are obtained using the iterations

\deqn{\hat{\delta}^{(k+1)}=\hat{\delta}^{(k)}+\Omega(\hat{\delta}^{(k)}%
)\frac{\partial l(\hat{\delta}^{(k)})}{\partial\delta}}{delta_(k+1) =%
delta_(k) + Omega(delta_k) * first derivative of log(delta_k)}

   where \eqn{\Omega (\hat\delta^{(k)})}{Omega(delta hat ^(k))} is some
   suitably chosen matrix.

   Iterations continue for \eqn{k\geq1}{k >= 1} until convergence is
   reached or the number of iterations \eqn{k} reaches a user specified
   upper limit on maximum iterations in which case they will stop. The
   convergence criterion used in our implementation is that based on
   \eqn{\eta}{eta}, the maximum of absolute values of the first
   derivatives.

   When \eqn{\eta}{eta} is less than a user specified value \code{grad}
   the iterations stop. There are two methods of optimization of the
   likelihood, Newton-Raphson and Fisher scoring. The method used is
   specified by the argument \code{method}. It should be noticed that if
   the initial value for parameters are not chosen well, the
   optimization of the likelihood might fail to converge. Care is needed
   when fitting mixed ARMA specifications because there is potential for
   the AR and MA parameters to be non-identifiable if the orders \eqn{p} and
   \eqn{q} are too large. Lack of identifiability manifests itself in the
   algorithm to optimize the likelihood failing to converge and/or the
   hessian being singular---check the warning messages and convergence
   error codes.
 }

\value{
  The function \code{\link{summary}} (i.e., \code{\link{summary.glarma}})
  can be used to obtain or print a summary of the results.

  The generic accessor functions \code{\link{coef}} (i.e.,
  \code{\link{coef.glarma}}), \code{\link{logLik}} (i.e.,
  \code{\link{logLik.glarma}}), \code{\link{fitted}} (i.e.,
  \code{\link{fitted.glarma}}), \code{\link{residuals}} (i.e.,
  \code{\link{residuals.glarma}}), \code{\link{nobs}} (i.e.,
  \code{\link{nobs.glarma}}), \code{\link{model.frame}} (i.e.,
  \code{\link{model.frame.glarma}}) and \code{\link{extractAIC}} (i.e.,
  \code{\link{extractAIC.glarma}}) can be used to extract various useful
  features of the value returned by \code{glarma}.

  \code{glarma} returns an object of class "glarma" with components:

  \item{delta}{a vector of coefficients for \code{beta}, \code{AR} and
    \code{MA}.}

  \item{logLik}{the loglikelihood of the specific distribution.}

  \item{logLikDeriv}{the derivative of the loglikelhood of the specified
    distribution.}

  \item{logLikDeriv2}{the second derivative of the loglikelihood of the
    specified distribution.}

  \item{eta}{the estimated linear predictor.}

  \item{mu}{the GLARMA estimated mean.}

  \item{fitted.values}{the GLARMA fitted values.}

  \item{residuals}{the residuals of the type specified.}

  \item{cov}{the estimated covariance matrix of the maximum
    likelihood estimators.}

  \item{phiLags}{vector of AR orders.}

  \item{thetaLags}{vector of MA orders.}

  \item{r}{the number of columns in the model matrix.}

  \item{pq}{the number of \code{phiLags} plus the number of
    \code{thetaLags}.}

  \item{null.deviance}{the deviance from the initial GLM fit.}

  \item{df.null}{the degrees of freedom from the initial GLM fit.}

  \item{y}{the \eqn{y} vector used in the GLARMA model.}

  \item{X}{the model matrix.}

  \item{offset}{the offset, \code{NULL} if there is no offset.}

  \item{type}{the distribution of the counts.}

  \item{method}{the method of iteration used.}

  \item{residType}{the type of the residuals returned.}

  \item{call}{the matched call.}

  \item{iter}{the number of iterations.}

  \item{errCode}{the error code; 0 indicating successful convergence of
    the iteration method, 1 indicating failure.}

  \item{WError}{error code for finiteness of \eqn{W}; 0 indicating all
    values of \eqn{W} are finite, 1 indicating at least one infinite
    value.}

  \item{min}{the minimum of the absolute value of the gradient.}

  \item{aic}{A version of Akaike's An Information Criterion, minus
   twice the maximized log-likelihood plus twice the number of
   parameters, computed by the aic component of the family. For
   binomial and Poisson families the dispersion is fixed at one and the
   number of parameters is the number of coefficients.}
}

\references{
  Dunsmuir, William T. M. and Scott, David J. (2015) The \pkg{glarma} Package
  for Observation-Driven Time Series Regression of Counts.
  \emph{Journal of Statistical Software}, \bold{67(7)}, 1--36.
  \doi{10.18637/jss.v067.i07}
}

\seealso{
  Additional examples may be found in \code{\link{Asthma}},
  \code{\link{OxBoatRace}}, \code{\link{RobberyConvict}}, and
  \code{\link{DriverDeaths}}.
}

\author{
  The original GLARMA routine for Poisson responses was developed in
  collaboration with Richard A. Davis and Ying Wang. The binomial
  response version was developed with the assistance of Haolan Lu. The
  extension to negative binomial response was carried out by Bo
  Wang. Daniel Drescher contributed to the initial structure of the
  software used as the basis of the package.

  The main author of the package is  "William T.M. Dunsmuir"
  <w.dunsmuir@unsw.edu.au>. Package development was carried out by
  Cenanning Li supervised by David J. Scott.
}

\examples{
### Example from Davis, Dunsmuir Wang (1999)
## MA(1,2,5), Pearson Residuals, Fisher Scoring
data(Polio)
y <- Polio[, 2]
X <- as.matrix(Polio[, 3:8])
glarmamod <- glarma(y, X, thetaLags = c(1,2,5), type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)
glarmamod
summary(glarmamod)

## Score Type (GAS) Residuals, Fisher Scoring
glarmamod <- glarma(y, X, thetaLags = c(1,2,5), type = "Poi", method = "FS",
                    residuals = "Score", maxit = 100, grad = 1e-6)
glarmamod
summary(glarmamod)

## Score Type (GAS)  Residuals, Newton Raphson
## Note: Newton Raphson fails to converge from GLM initial estimates.
## Setting up the initial estimates by ourselves
init.delta <- glarmamod$delta
beta <- init.delta[1:6]
thetaInit <- init.delta[7:9]

glarmamod <- glarma(y, X, beta = beta, thetaLags = c(1, 2, 5),
                    thetaInit = thetaInit, type ="Poi", method = "NR",
                    residuals = "Score", maxit = 100, grad = 1e-6)
glarmamod
summary(glarmamod)

## AR(1,5), Pearson Residuals, Fisher Scoring
glarmamod <- glarma(y, X, phiLags = c(1, 5), type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)
glarmamod
summary(glarmamod)



}


\keyword{GLARMA}
