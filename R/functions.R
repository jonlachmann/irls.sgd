# Title     : General functions
# Objective : General functions used by the package irls.sgd
# Created by: jonlachmann
# Created on: 2021-04-20

# Calculate the full deviance for a model
get_deviance <- function(beta, X, y, family) {
  mu <- family$linkinv(X %*% beta)
  sum(family$dev.resids(y, mu, rep(1,nrow(X))))
}