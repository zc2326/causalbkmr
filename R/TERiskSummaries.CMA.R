#' TERiskSummaries.CMA
#'
#'
#' Compare estimated Total Effect when all predictors are at a particular quantile to when all are at a second fixed quantile
#'
#' @param fit.TE total effect model fit regressing outcome on exposures, effect modifiers and confounders on outcome
#' @param e.y effect modifier for the outcome variable
#' @param e.y.names column name of the effect modifier for the outcome variable
#' @param qs vector of quantiles at which to calculate the overall risk summary
#' @param q.fixed a second quantile at which to compare the estimated h function
#' @param alpha 1-confidence interval
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#'
#'
#' @return a data frame containing the (posterior mean) estimate and posterior standard deviation of the total effect risk measures
#'
#'
#' @details
#' For guided examples, go to
#' https://zc2326.github.io/causalbkmr/articles/BKMRCMA_QuickStart.html
#'
#' @examples
#' \dontrun{
#' library(causalbkmr)
#' riskSummary10 = TERiskSummaries.CMA(fit.TE = fit.y.TE, e.y=e.y10, e.y.name = "E.Y", sel=sel)
#'
#' ggplot(riskSummary10,
#'        aes(quantile,
#'            est,
#'            ymin = est - 1.96 * sd,
#'            ymax = est + 1.96 * sd)) +
#'   geom_pointrange()
#' }
#'
#' @export
TERiskSummaries.CMA <- function(fit.TE,
                                e.y=NULL, e.y.names=NULL,
                                qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5,
                                alpha = 0.05, sel, seed = 122) {

  toreturn <- data.frame(quantile=qs,
                         est=rep(NA,times=length(qs)),
                         sd=rep(NA,times=length(qs)))
  fit <- fit.TE
  Z <- fit$Z
  if (!is.null(e.y.names)){
    Z = Z[,-which(e.y.names == colnames(Z))]
  }
  X <- fit$X
  X.predict <- matrix(colMeans(X),nrow=1)
  for(i in 1:length(qs)){
    quant <- qs[i]
    astar <- apply(Z, 2, quantile, q.fixed)
    a <- apply(Z, 2, quantile, quant)
    preds = TE.bkmr(a=a, astar=astar, e.y=e.y, fit.y.TE=fit.TE,
                    X.predict.Y=X.predict, alpha = alpha, sel=sel, seed=seed)
    toreturn[i, c(2,3)] = c(preds$est[,"mean"], preds$est[,"sd"])
  }
  return(toreturn)
}
