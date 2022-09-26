#' Single Variable Risk Summaries for MI BKMR fits
#' Compute summaries of the risks associated with a change in a single variable in \code{Z} from a single level (quantile) to a second level (quantile), for the other variables in \code{Z} fixed to a specific level (quantile) for MI BKMR fits
#'
#' @param BKMRfits  A list of multiple BKMR fits and that each of these fits were ran for the same number of MCMC iterations.
#' @param qs.diff vector indicating the two quantiles q_1 and q_2 at which to compute \code{h(z_{q2}) -h(z_{q1})}
#' @param q.fixed vector of quantiles at which to fix the remaining predictors in Z
#' @param q.alwaysfixed  the quantile values in the point which we want to keep fixed for all comparisons
#' @param index.alwaysfixed the index values in the point which we want to keep fixed for all comparisons
#' @param sel  selects which iterations of the MCMC sampler to use for inference
#' @param which.z  vector indicating which variables (columns of Z) for which the summary should be computed
#' @param z.names column names of the selected columns of Z in which.z
#' @param ... other arguments to pass on to the prediction function
#' @param method method for obtaining posterior summaries at a vector of new points. Options are"approx" and "exact"; defaults to "approx", which is faster particularly for large datasets
#'
#' @return a data frame containing the (posterior mean) estimate and posterior standard deviation of the  predictor risk measures
#'
#' @export
#' @details
#' For guided examples, go to https://zc2326.github.io/causalbkmr/articles/MI_BKMR.html
#' @examples
#' \dontrun{
#' library(causalbkmr)
#' data(BKMRfits10)
#' singvarrisk.MI.fixed <- SingVarRiskSummaries.MI(BKMRfits = BKMRfits10, which.z=c(1,3,4),
#' qs.diff = c(0.25, 0.75),  q.fixed = c(0.25, 0.50, 0.75),
#' q.alwaysfixed = 0.25, index.alwaysfixed = 2,
#' sel=sel.MI, method = "approx")
#'
#' ## plot the single variable dataframe for the MI fits
#' ggplot(singvarrisk.MI, aes(variable, est, ymin = est - 1.96*sd,
#'                            ymax = est + 1.96*sd, col = q.fixed)) +
#'   geom_hline(aes(yintercept=0), linetype="dashed", color="gray") +
#'   geom_pointrange(position = position_dodge(width = 0.75)) +
#'   coord_flip() + ggtitle("")+
#'   scale_x_discrete(name="Variable")+ scale_y_continuous(name="estimate")
#' }
#'
SingVarRiskSummaries.MI <- function(BKMRfits, which.z = 1:ncol(BKMRfits[[1]]$Z), qs.diff = c(0.25, 0.75), q.fixed = c(0.25, 0.50, 0.75), q.alwaysfixed = NULL, index.alwaysfixed = NULL, sel = NULL, z.names = colnames(BKMRfits[[1]]$Z), method="approx",...) {

  start.time <- proc.time()["elapsed"]
  if(is.null(z.names)) z.names <- paste0("z", 1:ncol(BKMRfits[[1]]$Z))
  Z.MI <- Z.complete.MI(BKMRfits)

  df <- dplyr::data_frame()
  for(i in seq_along(q.fixed)) {
    for(j in seq_along(which.z)) {
      risk <- VarRiskSummary.MI(whichz = which.z[j], BKMRfits=BKMRfits, Z.MI=Z.MI, qs.diff = qs.diff, q.fixed = q.fixed[i], q.alwaysfixed = q.alwaysfixed, index.alwaysfixed = index.alwaysfixed, sel = sel, method = method, ...)
      df0 <- dplyr::data_frame(q.fixed = q.fixed[i], variable = z.names[which.z[j]], est = risk["est"], sd = risk["sd"])
      df <- dplyr::bind_rows(df, df0)
    }

    end.time <- proc.time()["elapsed"]
    print(paste(i,"out of", length(q.fixed), "complete: ", round((end.time - start.time)/60, digits=2), "min run time" ))
  }
  df <- dplyr::mutate_(df, variable = ~factor(variable, levels = z.names[which.z]), q.fixed = ~as.factor(q.fixed))
  attr(df, "qs.diff") <- qs.diff
  df
}



VarRiskSummary.MI <- function(whichz = 1, BKMRfits, Z.MI, qs.diff = c(0.25, 0.75), q.fixed = 0.5, q.alwaysfixed = NULL, index.alwaysfixed = NULL, sel = NULL, method = "approx") {

  cc <- c(-1, 1)
  K <- length(BKMRfits)

  if(method=="exact") {
    preds.fun <- function(znew) ComputePostmeanHnew.exact.MI(fit = fit, y = y, Z = Z, X = X, Znew = znew, sel = sel)

    postmean.temp <- matrix(NA, nrow=length(sel)*K, ncol=2)
    postvar.temp  <- array(NA, dim = c(2,2,length(sel)*K))## 2 = nrow(newz)

    for(k in 1:K){
      fit <- BKMRfits[[k]]
      y <- fit$y
      Z <- fit$Z
      X <- fit$X

      point1 <- point2 <- apply(Z.MI, 2, quantile, q.fixed)
      point2[whichz] <- apply(Z.MI[, whichz, drop = FALSE], 2, quantile, qs.diff[2])
      point1[whichz] <- apply(Z.MI[, whichz, drop = FALSE], 2, quantile, qs.diff[1])

      ## if both q.alwaysfixed and index.alwaysfixed are specified,
      ## change the values in the point which we want to keep fixed for all comparisons
      if(!is.null(q.alwaysfixed) & !is.null(index.alwaysfixed)){
        point1[index.alwaysfixed] <- point2[index.alwaysfixed] <- apply(Z.MI[,index.alwaysfixed, drop=FALSE],2,quantile, q.alwaysfixed)
      }

      newz <- rbind(point1, point2)

      preds <- preds.fun(newz)

      postmean.temp[((k-1)*length(sel)+1):(length(sel)*k),] <- preds$postmean_mat
      postvar.temp[,,((k-1)*length(sel)+1):(length(sel)*k)] <- preds$postvar_arr

    }
    m  <- colMeans(postmean.temp)
    ve <- var(postmean.temp)
    ev <- apply(postvar.temp, c(1, 2), mean)
    v  <- ve + ev

    diff     <- drop(cc %*% m)
    diff.sd  <- drop(sqrt(cc %*% v %*% cc))
    toreturn <- c(est = diff, sd = diff.sd)

  } else if(method=="approx") {
    preds.fun <- function(znew) ComputePostmeanHnew.approx(fit = fit, y = y, Z = Z, X = X, Znew = znew, sel = sel)

    est.vec <- rep(NA, times=K)
    var.vec <- rep(NA, times=K)

    for(k in 1:K){
      fit <- BKMRfits[[k]]
      y <- fit$y
      Z <- fit$Z
      X <- fit$X

      point1 <- point2 <- apply(Z.MI, 2, quantile, q.fixed)
      point2[whichz] <- apply(Z.MI[, whichz, drop = FALSE], 2, quantile, qs.diff[2])
      point1[whichz] <- apply(Z.MI[, whichz, drop = FALSE], 2, quantile, qs.diff[1])

      ## if both q.alwaysfixed and index.alwaysfixed are specified,
      ## change the values in the point which we want to keep fixed for all comparisons
      if(!is.null(q.alwaysfixed) & !is.null(index.alwaysfixed)){
        point1[index.alwaysfixed] <- point2[index.alwaysfixed] <- apply(Z.MI[,index.alwaysfixed, drop=FALSE],2,quantile, q.alwaysfixed)
      }

      newz <- rbind(point1, point2)

      preds <- preds.fun(newz)

      est.vec[k] <- drop(cc %*% preds$postmean)
      var.vec[k] <- drop(cc %*% preds$postvar %*% cc)
    }
    if(K==1){toreturn <- c(est=est.vec, sd=sqrt(var.vec))}else{toreturn <- Rubin.MI(mean.vec = est.vec, variance.vec = var.vec)}
  } else stop("method must be one of c('approx', 'exact')")
  toreturn
}
