
#' Single Variable Risk Summaries for CMA
#
#' Compute summaries of the risks associated with a change in a single variable in \code{Z} from a single level (quantile) to a second level (quantile), for the other variables in \code{Z} fixed to a specific level (quantile)
#'
#' @param BKMRfits An object containing the results returned by a the kmbayes function
#' @param e.y effect modifier for the outcome variable
#' @param e.y.names column name of the effect modifier for the outcome variable
#' @param which.z vector indicating which variables (columns of Z) for which the summary should be computed
#' @param z.names optional vector of names for the columns of z
#' @param qs.diff vector indicating the two quantiles q_1 and q_2 at which to compute h(z_{q2}) -h(z_{q1})
#' @param q.fixed  vector of quantiles at which to fix the remaining predictors in Z
#' @param alpha 1-confidence interval
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#'
#'
#' @return a data frame containing the (posterior mean) estimate and posterior standard deviation of the single predictor risk measures
#' @export
#'
#' @examples
SingVarRiskSummaries.CMA <-function(BKMRfits,
                                    e.y = NULL, e.y.names = NULL,
                                    which.z = 1:ncol(BKMRfits$Z),
                                    z.names = colnames(BKMRfits$Z),
                                    qs.diff = c(0.25, 0.75),
                                    q.fixed = c(0.25, 0.50, 0.75),
                                    alpha = 0.05, sel, seed = 122) {
  Z = BKMRfits$Z
  if(is.null(z.names)) z.names <- paste0("z", 1:ncol(Z))
  if (!is.null(e.y.names)){
    Z = Z[,-which(e.y.names == colnames(Z))]
    which.z = 1:ncol(Z)
    z.names = colnames(Z)
  }
  df <- tibble::tibble()
  for(i in seq_along(q.fixed)) {
    for(j in seq_along(which.z)) {
      risk = VarRiskSummary.CMA(whichz = which.z[j],
                                BKMRfits = BKMRfits,
                                e.y = e.y, e.y.names = e.y.names,
                                qs.diff = qs.diff, q.fixed = q.fixed[i],
                                alpha = alpha, sel = sel, seed = seed)
      df0 <- tibble::tibble(q.fixed = q.fixed[i],
                            variable = z.names[which.z[j]],
                            est = risk["est"],
                            sd = risk["sd"])
      df <- dplyr::bind_rows(df, df0)
    }
  }
  df$variable <- factor(df$variable, levels = z.names[which.z])
  df$q.fixed = as.factor(df$q.fixed)
  attr(df, "qs.diff") <- qs.diff
  return(df)
}