#' Calculate Single Variable Risk Summaries when fixing multiple effect modifiers at certain levels
#' Compute summaries of the risks associated with a change in a single variable in \code{Z} from a single level (quantile) to a second level (quantile), for a set of effect modifiers (in \code{Z}) fixed to a specific level (quantile)
#'
#' @param list.fit.y.TE
#' @param which.z  vector indicating which variables (columns of Z) for which the summary should be computed, effect modifiers are not included
#' @param qs.diff vector indicating the two quantiles q_1 and q_2 at which to compute \code{h(z_{q2}) -h(z_{q1})}
#' @param q.fixed vector of quantiles at which to fix the remaining predictors in Z
#' @param q.alwaysfixed  the quantile values in the point which we want to keep fixed for all comparisons
#' @param EY.alwaysfixed.name names of all the effect modifiers that we want to fixed for all comparisons
#' @param sel selects which iterations of the MCMC sampler to use for inference
#' @param z.names column names of the selected columns of Z in which.z
#' @param method method for obtaining posterior summaries at a vector of new points. Options are"approx" and "exact"; defaults to "approx", which is faster particularly for large datasets
#' @param ...
#' @return a list of data frames containing the (posterior mean) estimate and posterior standard deviation of the  predictor risk measures, for each of the comparisons specified
#'
#' @export
#' @details
#' For guided examples, go to https://zc2326.github.io/causalbkmr/articles/BKMRCMA_Effectof_singleZ.html
#'
#'
SingVarRiskSummaries.fixEY <- function(list.fit.y.TE, which.z = 1:length(z.names),
                                       qs.diff = c(0.25, 0.75), q.fixed = c(0.25, 0.50, 0.75),
                                       q.alwaysfixed = NULL, EY.alwaysfixed.name = NULL,
                                       sel = NULL, z.names = colnames(BKMRfits[[1]]$Z), method="approx",...) {

  index.alwaysfixed <- which(colnames(list.fit.y.TE[[1]]$Z)   %in% EY.alwaysfixed.name )
  toreturn <- vector("list",  length(q.alwaysfixed))

  for(r in 1:length(q.alwaysfixed)){
    q.alwaysfixed_r <- q.alwaysfixed[r]
    toreturn[[r]] <- SingVarRiskSummaries.MI(BKMRfits=list.fit.y.TE,which.z =which.z,  qs.diff = qs.diff,
                                             q.fixed = q.fixed, q.alwaysfixed = q.alwaysfixed_r,
                                             index.alwaysfixed = index.alwaysfixed, sel = sel, method="approx")
  }
  toreturn
}
