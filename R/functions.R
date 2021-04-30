# Title     : General functions
# Objective : General functions used by the package irls.sgd
# Created by: jonlachmann
# Created on: 2021-04-20

#' Get the deviance for a model
#'
#' @param beta The vector of coefficients
#' @param X A matrix containing the covariates (including an intercept if one wants to use one)
#' @param y The dependent variable
#' @param family A glm family for the distribution to use, i.e. "binomial()"
#'
#' @export get_deviance
get_deviance <- function(beta, x, y, family) {
  mu <- family$linkinv(x %*% beta)
  sum(family$dev.resids(y, mu, rep(1,nrow(x))))
}

#' Create a plot to examine the convergence of the algorithm.
#'
#' @param result The results from irls.sgd.
#' @param true_coeffs A vector of true coefficients if they are available.
#' @param x A matrix of independent variables to estimate the true coefficients from.
#' @param y A vector with the dependent variable to estimate the true coefficients from.
#' @param family A glm family for the distribution to use, i.e. "binomial()" when estimating the true coefficients.
#'
#' @export conv.plot
conv.plot <- function (result, true_coeffs=NULL, x=NULL, y=NULL, family=binomial()) {
  if (is.null(true_coeffs) & !is.null(x) & !is.null(y)) {
    cat("Estimating full model...\n")
    true_coeffs <- glm.fit(x, y, family=gaussian())$coefficients
  }
  for (i in 1:ncol(result$irls_hist)) {
    if (!is.null(true_coeffs)) {
      multiplot(cbind(c(result$irls_hist[,i], result$sgd_hist[,i]), true_coeffs[i]),
                frame.plot=F, ylab=bquote(beta[.(i-1)]), xlab="Iteration")
    }
    else {
      multiplot(c(result$irls_hist[,i], result$sgd_hist[,i]),
                frame.plot=F, ylab=bquote(beta[.(i-1)]), xlab="Iteration")
    }
    abline(v=length(result$irls_hist[,i]))
  }
}