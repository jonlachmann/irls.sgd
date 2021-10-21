# Title     : IRLS with subsampling
# Objective : IRLS algorithm with stochastic subsampling
# Created by: jonlachmann
# Created on: 2021-04-13

#' Subsampling IRLS
#'
#' @param X A matrix containing the covariates (including an intercept if one wants to use one)
#' @param y The dependent variable
#' @param family A glm family for the distribution to use, i.e. "binomial()"
#' @param ctrl A list of control parameters
#' @details The control parameters that can be set are
#' \itemize{
#'  \item{"subs" }{The subsample proportion to use for each iteration}
#'  \item{"maxit" }{The maximum number of iterations to run for}
#'  \item{"tol" }{The relative tolerance to consider when measuring convergence (not really useful for subsampling)}
#'  \item{"cooling" }{The cooling schedule to use, 3 numericals. First is number of iterations to stay constant for, the second is the temperature to start the first cooling iteration at and third is the exponential decay.}
#'  \item{"expl" }{The explosion detection settings, 3 numericals. First is number of iterations to not detect at, the second is the relative change in deviance to consider an explosion and the third is by which factor to change the subsample size at an explosion.}
#' }
#'
#' @export irls
irls <- function (X, y, family, ctrl=list(subs=1, maxit=100, tol=1e-7, cooling = c(3,0.9,0.95), expl = c(3,1.5,1))) {
  if (length(ctrl$expl) != 3) stop("expl needs to be a numeric of length 3")
  if (length(ctrl$expl) != 3) stop("cooling needs to be a numeric of length 3")

  temp <- ctrl$cooling[2]
  nobs <- nrow(X)
  nvars <- ncol(X)
  w <- matrix(0.43, nobs, 1)

  # Get the initial subsample if we are doing subsampling
  sub_size <- nobs*ctrl$subs
  if (ctrl$subs != 1) {
    subsi <- sample.int(nobs, sub_size, replace=T)
  }
  else subsi <- 1:nobs

  # Get initial eta, mu and deviance
  eta <- family$linkfun((y[subsi]+0.5)/2)
  mu <- family$linkinv(eta)
  dev <- sum(family$dev.resids(y[subsi], family$linkinv(rowSums(X[subsi,,drop=F])), 1))/ctrl$subs

  # Initialize matrices for algorithm history
  devhist <- matrix(NA,ctrl$maxit,1)
  rankhist <- matrix(NA,ctrl$maxit,1)
  betahist <- matrix(NA,ctrl$maxit,nvars)
  devhist[1] <- dev
  betahist[1,] <- 0
  explosions <- 0

  # Start algorithm loop
  iter <- 1
  conv <- F
  while (!conv & iter < ctrl$maxit) {
    # Get var_mu
    var_mu <- family$variance(mu)
    # Get mu eta
    mu_eta <- family$mu.eta(eta)
    # Get z
    z <- eta + (y[subsi] - mu) / (mu_eta)
    # Get w
    w[subsi] <- sqrt((mu_eta^2)/var_mu)

    if (iter > 1) betaold <- beta
    # Do WLS
    #fit <- .Call(stats:::C_Cdqrls, X[subsi,,drop=F] * as.numeric(w[subsi,drop=F]), z * w[subsi,drop=F], 1e-7, check=FALSE)
    beta <- qrls(X[subsi,,drop=F], w[subsi,drop=F], z)
    # Extract betas and rank
    #rankhist[iter] <- fit$rank
    rankhist[iter] <- length(beta)
    #beta <- fit$coefficients

    # Do new subsampling
    if (ctrl$subs != 1) subsi <- sample_int_expj(nobs, sub_size, prob=(w+0.1))

    # Do cooling schedule
    if (iter > ctrl$cooling[1]+2) {
      beta <- temp * beta + (1-temp) * betaold
      temp <- temp * ctrl$cooling[3]
    }

    # Get eta
    eta <- X[subsi,,drop=F] %*% beta

    # Get mu
    mu <- family$linkinv(eta)

    # Calculate deviance
    devold <- dev
    dev <- sum(family$dev.resids(y[subsi], mu, 1))/ctrl$subs

    # Check convergence
    if ((abs(dev - devold)) / (0.1 + abs(dev)) < ctrl$tol) {
      conv <- T
    }
    # Check for exploding deviance in subsampling
    if (iter > ctrl$expl[1] && (dev-devold) / abs(devold) > ctrl$expl[2]) {
      # Reset beta, lower temperature and start counting explosions
      beta <- betahist[max(iter-1,1),]
      temp <- temp*0.5
      explosions <- explosions + 1
      if (explosions > 5) {
        print(paste0("5 exploding deviances at iteration ", iter))
        break
      }

      # Increase subsample size to avoid further explosions
      if (ctrl$subs != 1) {
        # Increase subsample size by increase factor expl[3]
        sub_size <- min(sub_size*ctrl$expl[3], nobs)
        subsi <- sample_int_expj(nobs, sub_size, prob=(w+0.1))
      }
      # Get eta
      eta <- X[subsi,,drop=F] %*% beta

      # Get mu
      mu <- family$linkinv(eta)

    } else explosions <- 0

    iter <- iter + 1
    betahist[iter,] <- beta
    devhist[iter] <- dev
  }
  return(list(coefficients=beta,
              iter=iter,
              rank=rankhist[iter-1],
              devhist=devhist[1:iter],
              betahist=betahist[1:iter, ,drop=F],
              rankhist=rankhist,
              finaltemp=temp,
              finalsub=sub_size
  ))
}