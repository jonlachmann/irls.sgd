#' @useDynLib irls.sgd, .registration=TRUE
NULL

#' Subsampling IRLS with SGD
#'
#' @param x A matrix containing the covariates (including an intercept if one wants to use one)
#' @param y The dependent variable
#' @param family A glm family for the distribution to use, i.e. "binomial()"
#' @param irls.control A list of control parameters for the Subsampling IRLS (see ?irls for details).
#' @param sgd.control A list of control parameters for the SGD (see ?sgd for details, note that sgd.control$start will be overwritten).
#'
#' @export irls.sgd
irls.sgd <- function (x, y, family,
                      irls.control=list(subs=1, maxit=100, tol=1e-7, cooling = c(3,0.9,0.95), expl = c(3,1.5)),
                      sgd.control=list(subs=1, alpha=0.0005), save_hist=F) {

  # Calculate the S-IRLS model
  irls_res <- irls(x, y, family, irls.control)

  # Set the results from the IRLS-S as the start for SGD
  sgd.control$start <- irls_res$coefficients

  # Calculate the SGD model
  sgd_res <- glm.sgd(x, y, family, sgd.control)

  # Get the final deviance measure
  deviance <- get_deviance(sgd_res$coefficients, x, y, family)

  # Format results and return
  results <- list(coefficients=sgd_res$coefficients,
              deviance=deviance,
              rank=irls_res$rank,
              final.irls.temp=irls_res$finaltemp,
              irls.iters=irls_res$iter,
              final.irls.subsize=irls_res$finalsub
  )
  if (save_hist) {
    results$irls_hist <- irls_res$betahist
    results$sgd_hist <- sgd_res$xhist
  }
  return(results)
}