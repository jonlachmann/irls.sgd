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