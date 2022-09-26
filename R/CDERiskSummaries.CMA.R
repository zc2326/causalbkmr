
#' CDERiskSummaries.CMA
#'
#' Compare estimated Controlled Direct Effect when all predictors are at a particular quantile to when all are at a second fixed quantile
#'
#' @param fit.y An object contatinint the results return by the kmbayes function, a model fit regressing outcome on exposures, effect modifiers, mediator and confounders on outcome
#' @param e.y effect modifier for the outcome variable
#' @param e.y.names column name of the effect modifier for the outcome variable
#' @param m.value values that the mediator is set to
#' @param m.quant values of the quantile that the mediator is set to
#' @param m.name column name of the mediator
#' @param qs vector of quantiles at which to calculate the overall risk summary
#' @param q.fixed a second quantile at which to compare the estimated h function
#' @param alpha 1-confidence interval
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#' @return a data frame containing the (posterior mean) estimate and posterior standard deviation of the CDE risk measures
#' @details
#' For guided examples, go to
#' https://zc2326.github.io/causalbkmr/articles/BKMRCMA_QuickStart.html
#'
#' @examples
#'
#' \dontrun{
#' CDEriskSummary10 = CDERiskSummaries.CMA(fit.y = fit.y, e.y = e.y10, e.y.name = "E.Y", m.name = "m", sel = sel)
#' }
#' @export
CDERiskSummaries.CMA <- function(fit.y,
                                 e.y=NULL, e.y.names=NULL,
                                 m.value = NULL, m.quant = c(0.1, 0.5, 0.75), m.name,
                                 qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5,
                                 alpha = 0.05, sel, seed = 122) {
  if (!is.null(m.value) & !is.null(m.quant)){
    m.quant = NULL # if both m.value and m.quant are specified, default set to m.value
  }
  df <- dplyr::tibble()
  fit <- fit.y
  Z <- fit$Z
  m = Z[,m.name]
  if (!is.null(e.y.names)){
    Z = Z[,-which(e.y.names == colnames(Z))]
  }
  Z = Z[,-which(m.name == colnames(Z))]
  X <- fit$X
  X.predict <- matrix(colMeans(X),nrow=1)
  for(i in seq_along(qs)) {
    print(paste("Running", i, "out of", length(qs), "quantile values:"))
    quant <- qs[i]
    astar <- apply(Z, 2, quantile, q.fixed)
    a <- apply(Z, 2, quantile, quant)
    preds = CDE.bkmr(a, astar, e.y, m.value, m.quant, fit.y, alpha, sel, seed)
    if (is.null(m.quant)){newm = m.value} else(newm = m.quant)
    df0 <-dplyr::data_frame(quantile = quant,
                            m = newm,
                            est = preds$est[,"mean"],
                            sd = preds$est[,"sd"])
    df <- dplyr::bind_rows(df, df0)
  }
  df$m = as.factor(df$m)
  df$quantile = as.factor(df$quantile)
  return(df)
}
