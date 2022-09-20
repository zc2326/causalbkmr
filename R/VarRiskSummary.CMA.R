#' VarRiskSummary.CMA
#'
#' A combination function of VarRiskSummary and riskSummary.approx for MI BKMR fits, compare estimated \code{h} function when a single variable (or a set of variables) is at the 75th versus 25th percentile, when all of the other variables are fixed at a particular percentile
#'
#' @param whichz vector indicating which variables (columns of Z) for which the summary should be computed
#' @param BKMRfits An object contatinint the results return by the kmbayes function
#' @param e.y effect modifier for the outcome variable
#' @param e.y.names column name of the effect modifier for the outcome variable
#' @param qs.diff vector indicating the two quantiles q_1 and q_2 at which to compute \code{h(z_{q2}) -h(z_{q1})}
#' @param q.fixed vector of quantiles at which to fix the remaining predictors in Z
#' @param alpha 1-confidence interval
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#'
#' @return a data frame containing the (posterior mean) estimate and posterior standard deviation of the  predictor risk measures
#'
#' @export
VarRiskSummary.CMA <-function(whichz = 1,
                              BKMRfits,
                              e.y = NULL, e.y.names = NULL,
                              qs.diff = c(0.25, 0.75), q.fixed = 0.5,
                              alpha, sel, seed) {

  fit <- BKMRfits
  y <- fit$y
  Z <- fit$Z
  X <- fit$X
  if (!is.null(e.y.names)){
    Z = Z[,-which(e.y.names == colnames(Z))]
  }
  X.predict <- matrix(colMeans(X),nrow=1)

  a <- astar <- apply(Z, 2, quantile, q.fixed)
  astar[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile, qs.diff[1])# fix z-th variable to the quantile
  a[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile, qs.diff[2])

  preds = TE.bkmr(a=a, astar=astar, e.y=e.y, fit.y.TE=fit,
                  X.predict.Y = X.predict, alpha = alpha, sel=sel, seed=seed)

  toreturn = c(preds$est[,"mean"], preds$est[,"sd"])
  names(toreturn) = c("est", "sd")
  return(toreturn)
}
