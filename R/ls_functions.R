# Latent Stratification Functions
# Elea McDonnell Feit, eleafeit@gmail.com
# 27 January 2021


#' Simulate Data
#'
#' Simulates the data from the latent stratification model with four strata.
#'
#' @param n total sample size, the default value is 100000.
#' @param p treatment proportions, the default value is 0.5.
#' @param piA proportion of strata A, the default value is 0.05.
#' @param piB proportion of strata B, the default value is 0.10.
#' @param muA1 mean for strata A that received treatment, Z = 1, the default value is 5.
#' @param muA0 mean for strata A that did not received treatment, Z = 0, the default value is 4.
#' @param muB1 mean for strata B that received treatment, Z = 1, the default value is 5.
#' @param sigma variance for all strata, the default value is 1.
#'
#' @details
#' The four strata are defined as: \cr
#' A = positive under treatment and control, or always buyer. \cr
#' B = positive under treatment only, or influenced buyer. \cr
#' C = never positive, or never buyer. \cr
#' The model assumes that those who are positive under control, also known as defiers only doesn't exist.
#' The model also assumes that all strata share the same variance. \cr
#' The outcome y of strata A and B are generated through mixture models of normal distributions.
#' Strata A is generated using two normal distributions, one with mean muA1 and the other muA0.
#' Strata B is generated using a normal distribution with mean muB1 and 0, representing those in strata B won't be positive without treatment.
#' Strata C is 0 at all times. \cr
#' For the data frame in the output, column z is the dummy variable for treatment. If z = 1, then the observation has received treatment.
#' If z = 0, then the observation has not received treatment.
#'
#' @return
#' A list containing a data frame, a numeric value, and two vectors. \cr
#' The data frame \emph{data} contains the outcome variable y, treatment dummy z, and strata,
#' and the mean-centered effects-coded dummies for strata.\cr
#' The numeric value \emph{ATE} is the true average treatment effect. \cr
#' The vector \emph{par} consists of the true parameters of which the data is simulated. \cr
#'
#' @examples
#' sim = sim_latent_strat(n=10000, piA=0.2, piB=0.1, muA1=5, muA0=4.5, muB1=3, sigma=0.3)
#' sim$par
#' # a vector piA 0.2 piB 0.1 muA1 5 muA0 4.5 mUB 3 sigma 0.3
#' sim$data
#' # a data frame containing outcome variable y, treatment dummy z, and strata. The first 5000 rows have z=1.
#' sim$ATE
#' # 4.0
#'
#' @export
sim_latent_strat <- function(n=100000, # total sample size
                             p=0.5, # treatment proportion
                             piA=0.05, piB=0.10, # segment probabilities
                             muA1=5, muA0=4, muB1=5, # positive segment means
                             sigma=1) # positive segment variances
{
  piC <- 1 - piA - piB
  # treatment assignment is completely randomized
  n1 <- round(n*p)
  n0 <- n - n1
  z <- c(rep(1, n1), rep(0, n0))
  # compute size of each strata with fixed sizes
  nB <- round(piB*n)
  nC <- round(piC*n)
  nA <- n-nB-nC
  # create strata assignment vector and randomize the order
  s <- c(rep("A", nA), rep("B", nB), rep("C", nC))
  s <- sample(s)
  # create true/false vectors that are used as indices for y
  y <- rep(NA, n)
  xA <- s=="A"
  xB <- s=="B"
  xC <- s=="C"
  y[xA] <- z[xA]*rnorm(nA, muA1, sigma) + (1-z[xA])*rnorm(nA, muA0, sigma)
  y[xB] <- z[xB]*rnorm(nB, muB1, sigma) + (1-z[xB])*0
  y[xC] <- 0
  # Mean-centered effects-coded dummies for strata (used for regression estimates of ATE)
  sB <- 1*xB - 1*xA
  sB <- sB - mean(sB)
  sC <- 1*xC - 1*xA
  sC <- sC - mean(sC)
  # Compute true ATE
  ATE = piA*(muA1-muA0) + piB*muB1 + piC*0

  return(list(data=data.frame(y=y, z=z, s=s, sB=sB, sC=sC),
       ATE=ATE,
       par=c(piA=piA, piB=piB, muA1=muA1, muA0=muA0, muB1=muB1, sigma=sigma)))
}

#' Log-likelihood
#'
#' Computes the log-likelihood for each of the four observational groups under the assumption of three strata and common variance.
#'
#' @param par vector c(piA, piB, muA1, muA0, muB1, sigma), c(piA, piB/(1-piA), muA1, muA0, muB1, sigma) if trans=TRUE.
#' @param data data frame containing columns y (positive outcome with zeros) and z (treatment).
#' @param trans boolean signifying if piB has been transformed.
#'
#' @details
#' For the input data frame, column z is the dummy variable for treatment. If z = 1, then the observation has received treatment.
#' If z = 0, then the observation has not received treatment. \cr
#' Sometimes piB is transformed to relative proportions from absolute proportions. This transformation allows the reparameterization
#' of the piA and piB to allow constraint bounds between 0 and 1 in the optimization procedure. \cr
#' The log likelihoods are calculated based on equation 11-14 in the paper. Note that the equations presented in the paper are for
#' normal likelihood and for a single individual, thus corresponding adjustments have been made in the code to calculate the
#' log likelihood for the group. \cr
#' If you receive the warning \emph{Error in ll_ls(): Numeric overruns}, this means that the log likelihood has exceeded R's value storage
#' capacity, therefore storing it as infinity values. In this case, the value will be negative infinity, as likelihood does not exceed 1.
#'
#' @return Log-likelihood for the latent stratification model.
#'
#' @examples
#' sim = sim_latent_strat(n=10000, piA=0.2, piB=0.1, muA1=5, muA0=4.5, muB1=3, sigma=0.3)
#' ll_ls(sim$par, sim$data)
#'
#' # if the strata proportions are in relative sizes
#' ll_ls(sim$par, sim$data, trans=TRUE)
#'
#' @export
ll_ls <- function(par, data, trans=FALSE) {
  # assign values from the par vector
  piA <- par[1]
  piB <- par[2]
  # if tranformed, then par == c(piA, piB/(1-piA), muA1, muA0, muB1, sigma)
  if (trans) piB <- piB * (1 - piA)
  piC <- 1 - piA - piB
  muA1 <- par[3]
  muA0 <- par[4]
  muB1 <- par[5]
  sigma <- par[6]

  if (piA < 0 | piB < 0 | piC < 0 | sigma < 0)
    stop("Error in ll_ls(): piA, piB, piC or sigma < 0")

  # pre-compute indexing (for speed)
  # todo: maybe put this outside the likelihood function so you only do it once?
  nx <- data$y==0 # no purchase
  nz <- data$z==0 # not treated
  xz <- !nx & !nz
  xnz <- !nx & nz
  nxz <- nx & !nz
  nxnz <- nx & nz

  # Compute log-likelihood for each of four observational groups
  ll <- sum(nxz) * log(piC) +
        sum(nxnz) * log(piB + piC)
  ll_xnz <- log(piA) + dnorm(data$y[xnz], mean=muA0, sd=sigma, log=TRUE)
  ll_xz <- log(piA * dnorm(data$y[xz], mean=muA1, sd=sigma) +
                 piB * dnorm(data$y[xz], mean=muB1, sd=sigma))
  if (sum(ll_xz == -Inf) | sum(ll_xnz == -Inf)) {
    warning("Error in ll_ls(): Numeric overruns")
  }
  ll <- ll + sum(ll_xz) + sum(ll_xnz)
  return(unname(ll))
}

#' Gradient of the log-likelihood
#'
#' Computes the first order partial derivative of the log-likelihood for the latent stratification model, with respect to
#' every variable in the vector \emph{par}.
#'
#' @param par vector c(piA, piB, muA1, muA0, muB1, sigma), c(piA, piB/(1-piA), muA1, muA0, muB1, sigma) if trans=TRUE.
#' @param data data frame containing columns y (positive outcome with zeros) and z (treatment).
#' @param trans boolean signifying if piB has been transformed.
#'
#' @details
#' For the input data frame, column z is the dummy variable for treatment. If z = 1, then the observation has received treatment.
#' If z = 0, then the observation has not received treatment. \cr
#' Sometimes piB is transformed to relative proportions from absolute proportions. This transformation allows the reparameterization
#' of the piA and piB to allow constraint bounds between 0 and 1 in the optimization procedure. \cr
#' The output vector is named, each representing the gradient taken with respect to that variable in the parameter.
#'
#' @return Gradient of the log-likilihood for the latent stratification model as a named vector.
#'
#' @examples
#' sim = sim_latent_strat(n=10000, piA=0.2, piB=0.1, muA1=5, muA0=4.5, muB1=3, sigma=0.3)
#' gr_ll_ls(sim$par, sim$data)
#'
#' # if the strata proportions are in relative sizes
#' gr_ll_ls(sim$par, sim$data, trans=TRUE)
#'
#' @export
gr_ll_ls <- function(par, data, trans=FALSE) {
  piA <- par[1]
  piB <- par[2]
  if (trans) piB <- piB * (1 - piA)
  piC <- 1 - piA - piB
  muA1 <- par[3]
  muA0 <- par[4]
  muB1 <- par[5]
  sigma <- par[6]

  if (piA < 0 | piB < 0 | piC < 0 | sigma < 0)
    stop("Error in gr_ll_ls(): piA, piB, piC or sigma < 0")

  # pre-compute indexing (for speed)
  # todo: maybe put this outside the likelihood function so you only do it once?
  x <- 1*(data$y > 0) # purchased
  z <- data$z # treated
  y_xz <- data$y[x & z] # influenced

  # compute normal density for positive, treated (xz) observations
  fA1 <- dnorm(y_xz, muA1, sigma) # only compute for "treated, purchase" group
  fB1 <- dnorm(y_xz, muB1, sigma)
  if (sum(is.infinite(c(fA1, fB1))))
    warning("Numeric overrun in normal density calculation in gr_ll_ls()")

  # partial derivatives
  dpiA <- sum( fA1/(piA*fA1 + piB*fB1) ) - sum( (1-x)*z )/(1-piA-piB) +
    sum(x*(1-z))/piA - sum((1-x)*(1-z))/(1-piA)
  dpiB <- sum( fB1/(piA*fA1 + piB*fB1) ) - sum( (1-x)*z )/(1-piA-piB)
  dmuA1 <- sum( (piA*fA1/(piA*fA1 + piB*fB1))*((y_xz-muA1)/sigma^2) )
  dmuA0 <- sum( x*(1-z)*(data$y-muA0)/sigma^2 )
  dmuB1 <- sum( (piB*fB1/(piA*fA1 + piB*fB1))*((y_xz-muB1)/sigma^2) )
  dsigma <- sum( (piA*fA1*((y_xz-muA1)^2-sigma^2)/sigma^3 +
                    piB*fB1*((y_xz-muB1)^2-sigma^2)/sigma^3 )/(piA*fA1 + piB*fB1) ) +
    sum( x*(1-z)*((data$y-muA0)^2-sigma^2)/sigma^3 )
  out <- c(piA=dpiA, piB=dpiB, muA1=dmuA1, muA0=dmuA0, muB1=dmuB1, sigma=dsigma)
  return(out)
}

#' Hessian Matrix
#'
#' Computes the second order partial derivative with respect to each of the \emph{par} variables, resulting in a Hessian matrix.
#'
#' @param par vector c(piA, piB, muA1, muA0, muB1, sigma), c(piA, piB/(1-piA), muA1, muA0, muB1, sigma) if trans=TRUE.
#' @param data data frame containing columns y (positive outcome with zeros) and z (treatment).
#' @param trans boolean signifying if piB has been transformed.
#'
#' @details
#' For the input data frame, column z is the dummy variable for treatment. If z = 1, then the observation has received treatment.
#' If z = 0, then the observation has not received treatment. \cr
#' Sometimes piB is transformed to relative proportions from absolute proportions. This transformation allows the reparameterization
#' of the piA and piB to allow constraint bounds between 0 and 1 in the optimization procedure. \cr
#' The returned Hessian is the second order derivative with respect to \eqn{\theta} where
#' \eqn{\theta} is in the order of piA, piB, muA1, muA0, muB1, and sigma.
#'
#' @return Hessian matrix for the latent stratification model.
#'
#' @examples
#' sim = sim_latent_strat(n=10000, piA=0.2, piB=0.1, muA1=5, muA0=4.5, muB1=3, sigma=0.3)
#' hes_ll_ls(sim$par, sim$data)
#'
#' # if the strata proportions are in relative sizes
#' hes_ll_ls(sim$par, sim$data, trans=TRUE)
#'
#' @export
hes_ll_ls <- function(par, data, trans=FALSE) {
  piA <- par[1]
  piB <- par[2]
  if (trans) piB <- piB * (1 - piA)
  piC <- 1 - piA - piB
  muA1 <- par[3]
  muA0 <- par[4]
  muB1 <- par[5]
  s <- par[6]

  if (piA < 0 | piB < 0 | piC < 0 | s < 0)
    stop("Error in gr_ll_ls(): piA, piB, piC or sigma < 0")

  hes <- matrix(NA, nrow=6, ncol=6)
  z <- data$z
  y <- data$y
  x <- 1*(y > 0)
  Q1 <- (piA*exp((muB1 - y)^2/(2*s^2)) + piB*exp((muA1 - y)^2/(2*s^2)))^2
  Q2 <- exp( (muA1^2+muB1^2-2*y*(muA1+muB1)+2*y^2) / (2*s^2) )
  Q3 <- exp((muB1-y)^2/(2*s^2))
  Q4 <- exp((muA1-y)^2/(2*s^2))
  # Hessian Matrix is symmetric, so double assignment is used.
  hes[1,1] <- - sum((x-1)*(z-1)) / (piA-1)^2 + sum(x*(z-1))/piA^2 +
    sum((x-1)*z)/(piA + piB - 1)^2 +
    - sum(x*z*exp((muB1-y)^2/s^2) / Q1)
  hes[1,2] <- hes[2,1] <- sum((x-1)*z)/(piA+piB-1)^2 - sum(x*z*Q2/Q1)
  hes[1,3] <- hes[3,1] <- -sum(x*z*piB*(muA1-y)*Q2/(Q1*s^2))
  hes[1,4] <- hes[4,1] <- 0
  hes[1,5] <- hes[5,1] <- sum(x*z*piB*(muB1-y)*Q2/(Q1*s^2))
  hes[1,6] <- hes[6,1] <- sum(x*z*piB*(muA1-muB1)*(muA1+muB1-2*y)*Q2/(Q1*s^3))
  hes[2,2] <- sum(z*(x-1))/(piA+piB-1)^2 - sum(x*z*exp((muA1-y)^2/s^2)/Q1)
  hes[2,3] <- hes[3,2] <- sum(x*z*piA*(muA1-y)*Q2/(Q1*s^2))
  hes[2,4] <- hes[4,2] <- 0
  hes[2,5] <- hes[5,2] <- sum(x*z*piA*(y-muB1)*Q2/(Q1*s^2))
  hes[2,6] <- hes[6,2] <- -sum(x*z*piA*(muA1-muB1)*(muA1+muB1-2*y)*Q2/(Q1*s^3))
  hes[3,3] <- sum(-Q3*piA*x*(Q3*piA*s^2+Q4*piB*(-muA1^2+s^2+2*muA1*y-y^2))*z/(Q1*s^4))
  hes[3,4] <- hes[4,3] <- 0
  hes[3,5] <- hes[5,3] <- sum(x*z*piA*piB*(muA1-y)*(y-muB1)*Q2 / (Q1*s^4))
  hes[3,6] <- sum(-Q3*piA*x*z*(y-muA1)*
                    (2*Q3*piA*s^2+Q4*piB*(muB1^2-muA1^2+2*s^2+2*muA1*y-2*muB1*y)) /
                    (Q1*s^5))
  hes[6,3] <- hes[3,6]
  hes[4,4] <- sum(x*(z-1)/s^2)
  hes[4,5] <- hes[5,4] <- 0
  hes[4,6] <- hes[6,4] <- sum(-2*x*(z-1)*(muA0-y)/s^3)
  hes[5,5] <- sum(-x*z*piB*Q4*(Q4*piB*s^2 + Q3*piA*(s^2+2*muB1*y-muB1^2-y^2))/(Q1*s^4))
  hes[5,6] <- sum(-x*z*piB*(y-muB1)*Q4*(2*Q4*piB*s^2+Q3*piA*(muA1^2-muB1^2+2*s^2-2*muA1*y+2*muB1*y))/(Q1*s^5))
  hes[6,5] <- hes[5,6]
  hes[6,6] <- sum((1/(Q1*s^6))*exp(-(muA1+muB1)*y/s^2)*x*
                    ( exp((muB1^2-muB1*y+y*(muA1+y))/s^2)*piA^2*s^2*
                        ( s^2+3*muA0^2*(z-1)-6*muA0*y*(z-1)-3*(y^2+muA1^2*z-2*muA1*y*z) ) +
                        exp((muA1^2-muA1*y+y*(muB1+y))/s^2)*piB^2*s^2*
                        (s^2+3*muA0^2*(z-1)-6*muA0*y*(z-1)-3*(y^2+muB1^2*z-2*muB1*y*z)) +
                        exp((muA1^2+muB1^2+2*y^2)/(2*s^2))*piA*piB*
                        (2*s^4 + 6*muA0^2*s^2*(z-1) - 12*muA0*s^2*y*(z-1) +
                           (muA1-muB1)^2*(muA1+muB1-2*y)^2*z-3*s^2*
                           (2*y^2+(muA1^2+muB1^2)*z-2*(muA1+muB1)*y*z))))
  return(hes)
}

#' Standard Error of ATE
#'
#' Estimates the standard error and the exponential of the standard error of ATE via the delta method.
#'
#' @param par vector c(piA, piB, muA1, muA0, muB1, sigma).
#' @param vcv variance-covariance matrix of the parameters, can be calculated using ls_vcv().
#'
#' @details The delta method estimates the variance by expanding the function of a random variable through Taylor approximation,
#' which can be expanded to vectorized calculations.
#' For a more detailed explanation, see https://www.stata.com/support/faqs/statistics/delta-method/.
#'
#' @examples
#' sim = sim_latent_strat(n=10000, piA=0.2, piB=0.1, muA1=5, muA0=4.5, muB1=3, sigma=0.3)
#' varcovar = ls_vcv(sim$par, sim$data, "hessian")
#' varATEldelta(sim$par, varcovar)
#'
#' @return standard error of the average treatment effect.
#' @export
varATEldelta <-function(par, vcv) {
  piA <- par[1]
  piB <- par[2]
  piC <- 1 - piA - piB
  muA1 <- par[3]
  muA0 <- par[4]
  muB1 <- par[5]
  sigma <- par[6]
  # ATE <- piA*(muA1-muA0) + piB*muB1
  # expATE <- piA*(exp(muA1)-exp(muA0)) + piB*exp(muB1)
  g1 <- c(muA1-muA0, muB1, piA, -piA, piB, 0) # Jacobian for ATE
  g2 <- c(exp(muA1)-exp(muA0), exp(muB1)-1, piA*exp(muA1), -piA*exp(muA0),
          piB*exp(muB1), 0) # Jacobian for expATE  #todo: fix this
  vars <- c(ATE=t(g1) %*% vcv %*% g1, expATE=t(g2) %*% vcv %*% g2)
  return(vars)
}

#' Starting Values
#'
#' Computes the starting values for proportion, mean, and variance for optimization.
#' The proportions \emph{pi} are transformed to relative sizes.
#'
#' @param data data frame containing cols y (positive outcome with zeros) and z (treatment).
#' @param rand boolean controlling whether random adjusted starting mean values are produced.
#'
#' @details
#' For the input data frame, column z is the dummy variable for treatment. If z = 1, then the observation has received treatment.
#' If z = 0, then the observation has not received treatment. \cr
#' By allowing rand=TRUE, the optim() in mle_ls() can be started in multiple places, providing
#' random adj values while remaining within reasonable bounds.
#'
#' @return starting values used for maximum-likelihood optimization.
#'
#' @examples
#' sim = sim_latent_strat(n=10000, piA=0.2, piB=0.1, muA1=5, muA0=4.5, muB1=3, sigma=0.3)
#' start_ls(sim$data)
#'
#' # if wish to have multiple random starting values
#' start_ls(sim$data, rand=TRUE)
#'
#' @export
start_ls <- function(data, rand=FALSE) {
  # computes starting values for mle optimization (transformed)
  bounds <- bounds_ls(data)
  n <- nrow(data)
  xz <- data$z==1 & data$y>0
  xnz <- data$z==0 & data$y>0
  nxz <- data$z==1 & data$y==0
  piA <- mean(xnz)/mean(data$z==0)
  piC <- mean(nxz)/mean(data$z==1)
  piB <- 1 - piA - piC
  muA0 <- mean(data$y[data$z==0 & data$y > 0])
  sigma <- 0.8*max(sd(data$y[xz]), sd(data$y[xnz]))
  if (rand) {
    adj <- sample(c(-1,1), 1)*runif(1,0.1,1)*rand*sigma/2
  } else {
    adj <- 0.1*sigma/2 #todo: modify mle_ls so that we start with muA1 > muB1 and muA1 < muB1
  }
  muA1 <- mean(data$y[xz]) + adj
  muB1 <- mean(data$y[xz]) - adj
  return(c(piA=piA, piB=piB/(piB + piC), muA1=muA1, muA0=muA0, muB1=muB1, sigma=sigma))
}

#' Bounds
#'
#' Computes parameter bounds to be used in L-BFGS-B mle optimization.
#'
#' @param data data frame containing cols y (positive outcome with zeros) and z (treatment).
#'
#' @details
#' For the input data frame, column z is the dummy variable for treatment. If z = 1, then the observation has received treatment.
#' If z = 0, then the observation has not received treatment. \cr
#' Each strata must have at least three observations.\cr
#' According to RDocumentation, the L-BFGS-B method for optimization is taken from Byrd et. al. (1995), allowing box constraints with
#' upper and lower bound.
#'
#' @return The maximum and lowest possible values for segment proportions, mean, and variance.
#'
#' @examples
#' sim = sim_latent_strat(n=10000, piA=0.2, piB=0.1, muA1=5, muA0=4.5, muB1=3, sigma=0.3)
#' bounds_ls(sim$data)
#'
#' @export
bounds_ls <- function(data) {
  # Computes parameter bounds to be used in L-BFGS-B mle optimization
  y_xz <- data$y[data$z==1 & data$y>0]
  y_xnz <- data$y[data$z==0 & data$y>0]
  # each strata must be at least 3 observations
  min_p <- 0.001
  max_p <- 0.999
  # bounds for mu_A1 and mu_B1
  min_mu <- min(y_xz)
  max_mu <- max(y_xz)
  # sigma
  min_sigma <- 0.8*min(sd(y_xnz), sd(y_xz))
  max_sigma <- max(sd(y_xnz), sd(y_xz))
  return(list(lower=c(min_p, min_p, min_mu, -Inf, min_mu, min_sigma),
       upper=c(max_p, max_p, max_mu, Inf, max_mu, max_sigma)))
}

#' Variance-Covariance Matrix
#'
#' @param par vector c(piA, piB, muA1, muA0, muB1, sigma).
#' @param data data frame containing cols y (positive outcome with zeros) and z (treatment).
#' @param method method used for computation: score, hessian, robust, and bootstrap.
#'
#' @details
#' Computes the variance-covariance matrix: \cr
#' if method = "hessian" then the standard errors are computed by the numeric hessian \cr
#' if method = "score" then standard errors are computed from the gradient \cr
#' if method = "robust" then white robust standard errors are computed \cr
#' if method = "bootstrap" then the standard errors are computed by bootstrap \cr
#' For the input data frame, column z is the dummy variable for treatment. If z = 1, then the observation has received treatment.
#' If z = 0, then the observation has not received treatment. \cr
#'
#' @return the variance-covariance matrix based on the specificed method.
#'
#' @examples
#' sim = sim_latent_strat(n=10000, piA=0.2, piB=0.1, muA1=5, muA0=4.5, muB1=3, sigma=0.3)
#' ls_vcv(sim$par, sim$data, "hessian")
#'
#' @export
ls_vcv <- function(par, data, method) {
  if (method=="score" | method=="robust") {
    B <- 0 # as defined in White 1982
    for (i in 1:nrow(data)) { # todo: vectorize this loop for speed?
      gr <- gr_ll_ls(par, data[i, ])
      B <- B +  gr %*% t(gr)
    }
  }
  if (method=="hessian" | method=="robust")
    A <- hes_ll_ls(par, data) # aka A in White 1982
  if (method=="hessian") {
    if (is.finite(det(-A))) {
      out <- solve(-A)
    } else {
      out <- NULL
    }
  }
  if (method=="score") {
    if (is.finite(det(B))) {
      out <- solve(B)
    } else {
      out <- NULL
    }
  }
  if (method=="robust") {
    if (is.finite(det(A))) {
      out <- solve(A) %*% B %*% solve(A)
    } else {
      out <- NULL
    }
  }
  if (method=="bootstrap") out <- NULL # todo
  return(out)
}

#' Maximum Likeihood Estimate
#'
#' Computes the maximum likelihood estimate for the latent stratification model,
#' optionally takes starting values in original parameter space.
#'
#' @param data data frame containing cols y (positive outcome with zeros) and z (treatment).
#' @param start vector starting values for parameters c(piA, piB, muA1, muA0, muB1, sigma).
#' @param starts number of starting values.
#' @param vcv the variance-covariance matrix of the data, can be calculated using ls_vcv().
#' @param quiet boolean controlling if the computation time should be printed after execution.
#'
#' @details
#' If starts=1, then the optimization is run once from the starting values.
#' If starts>1, then mle optimization is done with multiple starting values. \cr
#' The output is an object, which can be called using summary() to show the ATE parameters, their standard errors, and the maximum likelihood.
#' To access the variance-covariance matrix, use $vcv. \cr
#' For the input data frame, column z is the dummy variable for treatment. If z = 1, then the observation has received treatment.
#' If z = 0, then the observation has not received treatment. \cr
#' For the \empth{vcv} parameter: \cr
#' if vcv = "hessian" then the standard errors are computed by the numeric hessian \cr
#' if vcv = "score" then standard errors are computed from the gradient \cr
#' if vcv = "robust" then white robust standard errors are computed \cr
#' if vcv = "bootstrap" then the standard errors are computed by bootstrap \cr
#' The function uses L-BFGS-B method for optimization. According to RDocumentation, it is taken from Byrd et. al. (1995),
#' allowing box constraints with upper and lower bound. \cr
#' If quiet=TRUE, then the time to optimize and calculate the variance covariance matrix will be displayed along with the results.
#'
#' @return object containing the MLE results.
#'
#' @examples
#' sim = sim_latent_strat(n=10000, piA=0.2, piB=0.1, muA1=5, muA0=4.5, muB1=3, sigma=0.3)
#' mle_ls(sim$data)
#'
#' # if you wish to start the optimization in 3 places
#' mle_ls(sim$data, starts=3)
#'
#' # if you wish to identify the starting values yourself
#' startv = c(0.2, 0.1, 5, 4.5, 3, 0.3)
#' mle_ls(sim$data, start=startv)
#'
#' @export
mle_ls <- function(data, start=NULL, starts=1, vcv="hessian", quiet=FALSE) {
  start.time <- Sys.time()
  starts <- max(starts, 1) # at least two optimization runs
  out <- NULL
  s <- NULL
  ll <- rep(NA, starts)
  # if no starting value input, then calculate using start_ls()
  if (is.null(start)) {
    s[[1]] <- start_ls(data)
  } else {
    start[2] <-  start[2] / (1 - start[1]) # transform second parameter
    s[[1]] <- start
  }
  if (starts >= 2) {
    for (i in 2:starts) s[[i]] <- start_ls(data, rand=TRUE)
  }
  bounds <- bounds_ls(data)
  for (i in 1:starts) {
    out[[i]] <- optim(par=s[[i]], fn=ll_ls, gr=gr_ll_ls, trans=TRUE,
                      data=data, method="L-BFGS-B",
                      lower=bounds$lower, upper=bounds$upper,
                      control=list(fnscale=-1, parscale=c(0.1, 0.1, 1, 1, 1, 1),
                                   maxit=2000, trace=0))
    ll[i] <- out[[i]]$value
  }
  if (!quiet) {
    cat("Optimization runs took",
        round(Sys.time() - start.time, 1), "secs", fill=TRUE)
    cat("Greatest ll was achieved in run", which.max(ll), "at", max(ll),
        "(worst was", min(ll),")", fill=TRUE)
  }
  out <- out[[which.max(ll)]] # save the run with the minimum ll
  est <- out$par
  est[2] <- est[2] * (1 - est[1]) # untransform second parameter
  # computes ATE using estimated values
  ATE <- unname(est[1]*(est[3] - est[4]) + est[2]*est[5])
  expATE <- unname(est[1]*(exp(est[3]) - exp(est[4])) + est[2]*(exp(est[5])-1))
  pars <- data.frame(par=c("ATE", "expATE", names(est)),
                     est=c(ATE, expATE, unname(est)))
  if (out$convergence!=0)
    cat("Error in MLE convergence. See $convergence in output.", fill=TRUE)
  # compute standard errors
  # hes <- optimHess(par=est, fn=ll_ls, trans=FALSE, data=data, # not transformed parameters
  #                 control=list(fnscale=-1, parscale=c(10, 10, rep(1, 4))))
  vcv <- ls_vcv(est, data, method=vcv)
  if (is.null(vcv)) {
    pars$se <- rep(NA, 2 + 6)
  } else {
    pars$se <- sqrt(c(varATEldelta(est, vcv), diag(vcv)))
  }
  if (!quiet) {
    cat("Computing variance-covariance matrix took",
        round(Sys.time() - start.time, 1), " secs", fill=TRUE)
  }
  return(list(pars=pars, vcv=vcv, optim.out=out))  # todo: return as a ls object
  # read on how to define an object and write a summary function
  # look at t.test, summary just to display

  # summary to display par only for now, maximum ll, number of observation, soft warning for multistart*
}

#' Difference in means estimate of the ATE
#'
#' Computes the ATE by taking the difference in the means.
#'
#' @param data data frame containing cols y (positive outcome with zeros) and z (treatment).
#'
#' @details
#' For the input data frame, column z is the dummy variable for treatment. If z = 1, then the observation has received treatment.
#' If z = 0, then the observation has not received treatment. \cr
#' t.test() from base R can be used as an alternative.
#'
#' @return difference in means estimate of the ATE
#'
#' @examples
#' sim = sim_latent_strat(n=10000, piA=0.2, piB=0.1, muA1=5, muA0=4.5, muB1=3, sigma=0.3)
#' ATEd(sim$data)
#'
#' @export
ATEd <- function(data) {
  # difference in means estimate of the ATE
  ATE <- mean(data$y[data$z==1]) - mean(data$y[data$z==0])
  return(ATE)
}

#' Oracle model estimate of the ATE
#'
#' Computes the ATE under the assumption that all strata information are known and true.
#'
#' @param data data frame containing cols y (positive outcome with zeros) and z (treatment).
#'
#' @details
#' For the input data frame, column z is the dummy variable for treatment. If z = 1, then the observation has received treatment.
#' If z = 0, then the observation has not received treatment. \cr
#' This serves as a benchmark result under the ideal scenario where all strata are known, and can used to compare with other analysis results.
#'
#' @return oracle estimate of ATE
#'
#' @examples
#' sim = sim_latent_strat(n=10000, piA=0.2, piB=0.1, muA1=5, muA0=4.5, muB1=3, sigma=0.3)
#' ATEo(data)
#'
#' @export
ATEo <- function(data) {
  # oracle model (known strata) estimate of the ATE
  xA <- data$s=="A"
  xB <- data$s=="B"
  xC <- data$s=="C"
  nz <- data$z==0
  piA <- sum(xA)/nrow(data)
  piB <- sum(xB)/nrow(data)
  piC <- sum(xC)/nrow(data)
  #piD <- sum(data$s=="D")/nrow(data)
  ATE <- piA * (mean(data$y[xA & !nz]) - mean(data$y[xA & nz])) +
    piB * (mean(data$y[xB & !nz]) - mean(data$y[xB & nz])) +
    piC * (mean(data$y[xC & !nz]) - mean(data$y[xC & nz]))
  #piD * (mean(data$y[data$s=="D" & data$z==1]) - mean(data$y[data$s=="D" & data$z==0]))
  return(ATE)
}


