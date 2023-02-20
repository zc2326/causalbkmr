#' Calculate overall risk summaries when fixing multiple effect modifiers at certain levels
#'
#' @param BKMRfits  The BKMR model fit as a 'List' form.
#' @param qs  vector of quantiles at which to calculate the overall risk summary
#' @param q.fixed  a second quantile at which to compare the estimated {h} function
#' @param q.alwaysfixed  the quantile values in the point which we want to keep fixed for all comparisons
#' @param EY.alwaysfixed.name names of all the effect modifiers that we want to fixed for all comparisons
#' @param sel selects which iterations of the MCMC sampler to use for inference
#' @param method method for obtaining posterior summaries at a vector of new points. Options are"approx" and "exact"; defaults to "approx", which is faster particularly for large datasets
#'
#' @importFrom stats var
#' @return a list of data frame containing the (posterior mean) estimate and posterior standard deviation of the overall risk measures, for each of the comparisons specified
#'
#' @details
#' For guided examples, go to https://zc2326.github.io/causalbkmr/articles/BKMRCMA_Effectof_singleZ.html
#' @export
#'
OverallRiskSummaries.fixEY <- function(BKMRfits, qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5,
                                       q.alwaysfixed = NULL,  EY.alwaysfixed.name = NULL,
                                       sel = NULL, method = "approx") {
  index.alwaysfixed <- which(colnames(list.fit.y.TE[[1]]$Z)   %in% EY.alwaysfixed.name )
  toreturn <- vector("list",  length(q.alwaysfixed))

  for(r in 1:length(q.alwaysfixed)){
    q.alwaysfixed_r <- q.alwaysfixed[r]

    toreturn[[r]] <- OverallRiskSummaries.MI(BKMRfits=list.fit.y.TE, qs = qs,
                                             q.fixed = q.fixed, q.alwaysfixed = q.alwaysfixed_r,
                                             index.alwaysfixed = index.alwaysfixed, sel = sel, method="approx")
  }
  toreturn
}
