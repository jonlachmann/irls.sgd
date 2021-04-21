#' Subsampling IRLS with SGD
#'
#' @param X A matrix containing the covariates (including an intercept if one wants to use one)
#' @param y The dependent variable
#' @param family A glm family for the distribution to use, i.e. "binomial()"
#' @param irls.control A list of control parameters for the Subsampling IRLS
#' @param sgd.control A list of control parameters for the SGD
#'
#' @export irls.sgd
irls.sgd <- function (X, y, family,
                      irls.control=list(subs=1, maxit=100, tol=1e-7, cooling = c(3,0.9,0.95), expl = c(3,1.5)),
                      sgd.control=list()) {
  # Calculate the IRLS-S model
  irls_res <- irls(X, y, family, irls.control)
  # Set the results from the IRLS-S as the start for SGD
  sgd.control$start <- irls_res$betahist[nrow(irls_res$betahist),]
  # Make sure that we shuffle the data for SGD to avoid problems where there is an inherent order in the data
  sgd.control$shuffle <- TRUE
  # Calculate the SGD model
  sgd_res <- sgd(X, y, model="glm", model.control=list(family=family$family), sgd.control)
  # Get the final deviance measure
  deviance <- get_deviance(sgd_res$coefficients, X, y, family)
  # Format results and return
  return(list(coefficients=sgd_res$coefficients,
              deviance=deviance,
              rank=irls_res$rank,
              final.irls.temp=irls_res$finaltemp,
              irls.iters=irls_res$iter,
              final.irls.subsize=irls_res$finalsub
  ))
}