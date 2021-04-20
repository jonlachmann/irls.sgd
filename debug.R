nvars <- 15
nobs <- 10^6
full_model_count <- 2^nvars

#install.packages("../irls.sgd_1.0.tar.gz")

library(irls.sgd)
library(mvtnorm)

# Generate the data to use
{
  set.seed(1911)
  covmat <- matrix(rnorm(nvars^2, runif(nvars^2, -5,5), runif(nvars^2, 0, 5)), nvars)
  covmat <- covmat %*% t(covmat)

  million_x <- cbind(1, matrix(rmvnorm(nobs, runif(nvars, -5,5), covmat), nobs))
  covars <- sample.int(nvars, 8)
  betas <- runif(9, -10, 10)
  million_y_g <- million_x[,c(1,covars)] %*% betas + rnorm(nobs, 0, 3)
  million_y_l <- rbinom(nobs, 1, (1/(1+exp(-million_y_g))))
}

logistic.loglik.aic <- function (y, x, model, complex, params) {
  suppressWarnings({mod <- glm.fit(as.matrix(x[,model]), y, family=binomial())})
  ret <- -(mod$deviance/2) - mod$rank
  return(ret)
}

logistic.loglik.aic.irlssgd <- function (y, x, model, complex, params) {
  mod <- irls.sgd(as.matrix(x[,model]), y, binomial(),
            irls.control=list(subs=params$subs, maxit=75, tol=1e-7, cooling = c(3,0.9,0.95), expl = c(3,1.5)))
  return(-(mod$deviance/2) - mod$rank)
}

# Calculate the full model set using 3% at each iteration
full_10K_sub_3 <- vector("list", full_model_count)
progress <- 0
for (i in 1:full_model_count) {
  print(i)
  modelvector <- as.logical(c(T,intToBits(i)[1:15]))
  loglik <- logistic.loglik.aic.irlssgd(million_y_l[1:10000], million_x[1:10000,], modelvector, NULL, list(subs = 0.03))
  full_10K_sub_3[[i]] <- list(prob=NA, model=modelvector[-1], crit=loglik, alpha=NA)
  if (i %% floor(full_model_count/40) == 0) progress <- print.progressbar(progress, 40)
}