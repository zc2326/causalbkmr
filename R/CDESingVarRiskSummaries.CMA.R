#' Single Variable Risk Summaries
#
#' Compute summaries of the risks associated with a change in a single variable in \code{Z} from a single level (quantile) to a second level (quantile), for the other variables in \code{Z} fixed to a specific level (quantile)
#'
#' @param BKMRfits An object contatinint the results return by the kmbayes function
#' @param e.y effect modifier for the outcome variable
#' @param e.y.names column name of the effect modifier for the outcome variable
#' @param which.z vector indicating which variables (columns of Z) for which the summary should be computed
#' @param z.names optional vector of names for the columns of z
#' @param m.value values that the mediator is set to
#' @param m.quant values of the quantile that the mediator is set to
#' @param m.name column name of the mediator
#' @param qs.diff vector indicating the two quantiles q_1 and q_2 at which to compute h(z_{q2})-h(z_{q1})
#' @param q.fixed a second quantile at which to compare the estimated h function
#' @param alpha 1-confidence interval
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#'
#' @return a data frame containing the (posterior mean) estimate and posterior standard deviation of the singlepredictor CDE risk measures
#'
#' @details
#' For guided examples, go to
#' https://zc2326.github.io/causalbkmr/articles/BKMRCMA_QuickStart.html
#'
#'
#' @export
CDESingVarRiskSummaries.CMA <-function(BKMRfits,
                                       e.y = NULL, e.y.names = NULL,
                                       which.z = 1:ncol(BKMRfits$Z),
                                       z.names = NULL,
                                       m.value = NULL, m.quant = c(0.1, 0.5, 0.75), m.name,
                                       qs.diff = c(0.25, 0.75),
                                       q.fixed = c(0.25, 0.50, 0.75),
                                       alpha = 0.05, sel, seed = 122) {
  if (!is.null(m.value) & !is.null(m.quant)){
    m.quant = NULL
  }
  Z = BKMRfits$Z
  if (!is.null(e.y.names)){
    Z = Z[,-which(e.y.names == colnames(Z))]
  }
  Z = Z[,-which(m.name == colnames(Z))]
  if(is.null(z.names)) z.names <- paste0("z", 1:ncol(Z))

  which.z = 1:ncol(Z)
 # z.names = colnames(Z)
  df <- tibble::tibble()
  for(i in seq_along(q.fixed)) {
    for(j in seq_along(which.z)) {
      risk = CDEVarRiskSummary.CMA(whichz = which.z[j],
                                   BKMRfits = BKMRfits,
                                   e.y = e.y, e.y.names = e.y.names,
                                   m.value = m.value, m.quant = m.quant, m.name = m.name,
                                   qs.diff = qs.diff, q.fixed = q.fixed[i],
                                   alpha = alpha, sel = sel, seed = seed)
      df0 <- tibble::tibble(q.fixed = q.fixed[i],
                            m.fixed = rownames(risk),
                            variable = z.names[which.z[j]],
                            est = risk[,"est"],
                            sd = risk[,"sd"])
      df <- dplyr::bind_rows(df, df0)
    }
  }
  df$m.fixed = gsub("CDE","", df$m.fixed)
  df$m.fixed = as.factor(df$m.fixed)
  df$variable = as.factor(df$variable)
  # df$variable <- factor(df$variable, levels = z.names[which.z])
  df$q.fixed = as.factor(df$q.fixed)
  attr(df, "qs.diff") <- qs.diff
  return(df)
}
