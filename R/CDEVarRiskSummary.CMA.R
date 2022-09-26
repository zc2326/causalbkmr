
#' Single variable plot for CDE.
#'
#' A combination function of VarRiskSummary and riskSummary.approx for MI BKMR fits.
#' Compare estimated \code{h} function when a single variable (or a set of variables) is at the 75th versus 25th percentile, when all of the other variables are fixed at a particular percentile

#' @param whichz vector indicating which variables (columns of Z) for which the summary should be computed
#' @param BKMRfits An object contatinint the results return by the kmbayes function
#' @param e.y effect modifier for the outcome variable
#' @param e.y.names column name of the effect modifier for the outcome variable
#' @param m.value values that the mediator is set to
#' @param m.quant values of the quantile that the mediator is set to
#' @param m.name column name of the mediator
#' @param qs.diff vector indicating the two quantiles q_1 and q_2 at which to compute h(z_{q2})-h(z_{q1})
#' @param q.fixed a second quantile at which to compare the estimated h function
#' @param alpha 1-confidence interval
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#'
#' @return a data frame containing the (posterior mean) estimate and posterior standard deviation of the single predictor CDE risk measures
#' @details
#' For guided examples, go to https://zc2326.github.io/causalbkmr/articles/BKMRCMA_QuickStart.html
#' @export
CDEVarRiskSummary.CMA <-function(whichz = 1,
                                 BKMRfits,
                                 e.y = NULL, e.y.names = NULL,
                                 m.value = NULL, m.quant = c(0.1, 0.5, 0.75), m.name,
                                 qs.diff = c(0.25, 0.75), q.fixed = 0.5,
                                 alpha, sel, seed) {

  # if (!is.null(m.value) & !is.null(m.quant)){
  #   m.quant = NULL
  # }
  fit <- BKMRfits
  y <- fit$y
  Z <- fit$Z
  # X <- fit$X
  m = Z[,m.name]
  if (!is.null(e.y.names)){
    Z = Z[,-which(e.y.names == colnames(Z))]
  }
  Z = Z[,-which(m.name == colnames(Z))]
  # X.predict <- matrix(colMeans(X),nrow=1)

  a <- astar <- apply(Z, 2, quantile, q.fixed)
  astar[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile, qs.diff[1]) # fix z-th variable to the quantile
  a[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile, qs.diff[2])

  preds = CDE.bkmr(a=a, astar=astar, e.y=e.y, m.value = m.value, m.quant = m.quant,
                   fit.y=fit, alpha = alpha, sel=sel, seed=seed)
  toreturn = preds$est[,c("mean","sd")]
  colnames(toreturn) = c("est", "sd")
  return(toreturn)
}
