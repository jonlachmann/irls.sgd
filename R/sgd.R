# Title     : TODO
# Objective : TODO
# Created by: jonlachmann
# Created on: 2021-04-21

#' Gradient for GLM
#'
#' @param x A matrix containing the covariates (including an intercept if one wants to use one).
#' @param y The dependent variable.
#' @param theta The coefficient vector to use.
#' @param family A glm family for the distribution to use, i.e. "binomial()".
#'
#' @export glm.grad
glm.grad <- function (x, y, theta, family) {
  # TODO: Does not work for the Gamma family yet
  m <- length(y)
  eta <- family$linkinv(x%*%theta)
  return(t(x)%*%(eta-y)/m)
}

#' (Batch) (Stochastic) Gradient Descent for GLM
#'
#' @param x A matrix containing the covariates (including an intercept if one wants to use one).
#' @param y The dependent variable.
#' @param family A glm family for the distribution to use, i.e. "binomial()".
#' @param sgd.ctrl control parameters to pass to the sgd function (see ?sgd for details).
#'
#' @export glm.sgd
glm.sgd <- function (x, y, family, sgd.ctrl=NULL) {
  grad.fun <- function (theta, data) {
    glm.grad(data[,-1,drop=F], data[,1], theta, family)
  }
  sgd_res <- sgd(grad.fun, data=cbind(y,x), sgd.ctrl)
  sgd_res$coefficients <- sgd_res$x
  return(sgd_res)
}

#' (Batch) (Stochastic) Gradient Descent
#'
#' @param gradient The gradient of the objective function.
#' @param data The dataset to use for the gradient, if applicable, if NULL, no data is used.
#' @param ctrl A list of control parameters.
#' @details The control parameters that can be set are
#' \itemize{
#'  \item{"start" }{The starting point, NULL means random around 0.}
#'  \item{"alpha" }{The initial step size, defaults to 0.0005.}
#'  \item{"decay" }{The exponential decay of the step size per iteration, defaults to 0.999.}
#'  \item{"subs" }{The subsample percentage to use at each iteration, defaults to 1 (100\%).}
#'  \item{"maxit" }{The maximum number of iterations to run for, defaults to 1000.}
#'  \item{"histfreq" }{The number of iterations between stored history frames, defaults to 50.}
#' }
#'
#' @export sgd
sgd <- function (gradient, data=NULL, ctrl=list(start, alpha, decay, subs, maxit, histfreq)) {
  # Initialize coefficients from start or randomly
  if (is.null(ctrl$start)) x <- rnorm(ncol(data)-1)
  else x <- ctrl$start
  nvars <- length(x)

  # Set up default values if not set
  if (is.null(ctrl$alpha)) ctrl$alpha <- 0.0005
  if (is.null(ctrl$decay)) ctrl$decay <- 0.999
  if (is.null(ctrl$subs)) ctrl$subs <- 1
  if (is.null(ctrl$maxit)) ctrl$maxit <- 1000
  if (is.null(ctrl$histfreq)) ctrl$histfreq <- 50

  # Get observation count and set subsample size
  if (!is.null(data)) {
    nobs <- nrow(data)
    sub_size <- nobs*ctrl$subs
  }

  # Set decay
  if (is.null(ctrl$decay)) ctrl$decay <- (maxit-2)/maxit

  # Set up matrix for history
  xhist <- matrix(NA, ctrl$maxit/ctrl$histfreq+1, nvars)
  xhist[1,] <- x

  # Run (S)GD
  for (i in 1:ctrl$maxit) {
    if (ctrl$subs != 1) sub_idx <- sample(1:nobs, sub_size, replace=T)
    else sub_idx <- 1:nobs
    ctrl$alpha <- ctrl$alpha * ctrl$decay
    if (!is.null(data)) x <- x - ctrl$alpha * gradient(x, data[sub_idx,,drop=F])
    else x <- x - ctrl$alpha * gradient(x)
    if (i %% ctrl$histfreq == 0) xhist[i/ctrl$histfreq+1,] <- x
  }
  return(list(x=x, xhist=xhist))
}
