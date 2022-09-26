
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
#' @details
#' For guided examples, go to
#' https://zc2326.github.io/causalbkmr/articles/BKMRCMA_QuickStart.html
#'
#' @examples
#' \dontrun{
#' library(causalbkmr)
#' dat <-  cma_sampledata(N=300, L=3, P=3, scenario=1, seed=7)
#'
#' A <- cbind(dat$z1, dat$z2, dat$z3)
#' X <- cbind(dat$x3)
#' y  <- dat$y
#' m  <- dat$M
#'
#' E.M <- NULL
#' E.Y <- dat$x2
#'
#' Z.M <- cbind(A,E.M)
#' Z.Y <- cbind(A, E.Y)
#' Zm.Y <- cbind(Z.Y, m)
#'
#' set.seed(1)
#' fit.y <- kmbayes(y=y, Z=Zm.Y, X=X, iter=5000, verbose=TRUE, varsel=FALSE)
#' #save(fit.y,file="bkmr_y.RData")
#'
#' set.seed(2)
#' fit.y.TE <- kmbayes(y=y, Z=Z.Y, X=X, iter=5000, verbose=TRUE, varsel=FALSE)
#' #save(fit.y.TE,file="bkmr_y_TE.RData")
#'
#' set.seed(3)
#' fit.m <- kmbayes(y=m, Z=Z.M, X=X, iter=5000, verbose=TRUE, varsel=FALSE)
#' #save(fit.m,file="bkmr_m.RData")
#'
#' X.predict <- matrix(colMeans(X),nrow=1)
#' astar <- c(apply(A, 2, quantile, probs=0.25))
#' a <- c(apply(A, 2, quantile, probs=0.75))
#'
#' e.y10 = quantile(E.Y, probs=0.1)
#' e.y90 = quantile(E.Y, probs=0.9)
#'
#' sel<-seq(2500,5000,by=5)
#'
#'
#' risks.singvar10 = SingVarRiskSummaries.CMA(BKMRfits = fit.y.TE, which.z = 1:3,
#' e.y=e.y10, e.y.names="E.Y",
#' sel=sel)
#' ggplot(risks.singvar10, aes(variable, est, ymin = est - 1.96*sd,
#'                             ymax = est + 1.96*sd, col = q.fixed)) +
#'   geom_pointrange(position = position_dodge(width = 0.75)) +
#'   coord_flip()
#'
#' }
SingVarRiskSummaries.CMA <-function(BKMRfits,
                                    e.y = NULL, e.y.names = NULL,
                                    which.z = 1:ncol(BKMRfits$Z),
                                    z.names = NULL,
                                    qs.diff = c(0.25, 0.75),
                                    q.fixed = c(0.25, 0.50, 0.75),
                                    alpha = 0.05, sel, seed = 122) {
  Z = BKMRfits$Z
  if (!is.null(e.y.names)){
    Z = Z[, -which(e.y.names == colnames(Z))]
    which.z = 1:ncol(Z)
    #z.names = colnames(Z)
  }

  if(is.null(z.names)) {
    z.names <- paste0("z", 1:ncol(Z))
  }


  df <- dplyr::data_frame()
  for(i in seq_along(q.fixed)) {
    for(j in seq_along(which.z)) {
      risk = VarRiskSummary.CMA(whichz = which.z[j],
                                BKMRfits = BKMRfits,
                                e.y = e.y, e.y.names = e.y.names,
                                qs.diff = qs.diff, q.fixed = q.fixed[i],
                                alpha = alpha, sel = sel, seed = seed)
      df0 <- dplyr::data_frame(q.fixed = q.fixed[i],
                               variable = z.names[which.z[j]],
                               est = risk["est"],
                               sd = risk["sd"])
      df <- dplyr::bind_rows(df, df0)


    }
  }

  df <- dplyr::mutate_(df, variable = ~factor(variable,levels = z.names[which.z]),
                       q.fixed = ~as.factor(q.fixed))
  attr(df, "qs.diff") <- qs.diff
  return(df)
}
