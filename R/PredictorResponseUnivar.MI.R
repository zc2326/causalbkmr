
#' Plot univariate predictor-response function on a new grid of point for MI BKMR
#'
#' @param BKMRfits  A list of multiple BKMR fits and that each of these fits were ran for the same number of MCMC iterations.
#' @param which.z vector identifying which predictors (columns of \code{Z}) should be plotted
#' @param ngrid number of grid points to cover the range of each predictor (column in \code{Z})
#' @param q.fixed vector of quantiles at which to fix the remaining predictors in \code{Z}
#' @param sel  logical expression indicating samples to keep; defaults to keeping the second half of all samples
#' @param min.plot.dist specifies a minimum distance that a new grid point needs to be from an observed data point in order to compute the prediction; points further than this will not be computed
#' @param center  flag for whether to scale the exposure-response function to have mean zero
#' @param method method for obtaining posterior summaries at a vector of new points. Options are "approx" and "exact"; defaults to "approx", which is faster particularly for large datasets
#' @param ... other arguments to pass on to the prediction function
#'
#' @importFrom dplyr %>%
#' @importFrom fields rdist
#' @return a long data frame with the predictor name, predictor value, posterior mean estimate, and posterior standard deviation
#' @export
#'
#' @examples
PredictorResponseUnivar.MI <- function(BKMRfits, which.z = 1:ncol(BKMRfits[[1]]$Z), ngrid = 50, q.fixed = 0.5, sel = NULL, min.plot.dist = Inf, center = TRUE, method = "approx", ...) {

  start.time <- proc.time()["elapsed"]
  Z.MI <- Z.complete.MI(BKMRfits)
  z.names <- colnames(Z.MI)

  df <- dplyr::data_frame()
  for(i in which.z) {
    res <- PredictorResponseUnivarVar.MI(whichz = i, BKMRfits = BKMRfits, Z.MI = Z.MI, ngrid = ngrid, q.fixed = q.fixed, sel = sel, min.plot.dist = min.plot.dist, center = center, method = method, ...)
    df0 <- dplyr::mutate(res, variable = z.names[i]) %>%
      dplyr::select_(~variable, ~z, ~est, ~se)
    df <- dplyr::bind_rows(df, df0)

    end.time <- proc.time()["elapsed"]
    print(paste(i,"out of", length(which.z), "complete: ", round((end.time - start.time)/60, digits=2), "min run time" ))
  }
  df$variable <- factor(df$variable, levels = z.names[which.z])
  df
}


PredictorResponseUnivarVar.MI <- function(whichz = 1, BKMRfits, Z.MI, ngrid = 50, q.fixed = 0.5, sel = NULL, min.plot.dist = Inf, center = TRUE, method = "approx",...) {

  K <- length(BKMRfits)
  ord <- c(whichz, setdiff(1:ncol(Z.MI), whichz))
  z1 <- seq(min(Z.MI[,ord[1]]), max(Z.MI[,ord[1]]), length = ngrid)
  z.others <- lapply(2:ncol(Z.MI), function(x) quantile(Z.MI[,ord[x]], q.fixed))
  z.all <- c(list(z1), z.others)
  newz.grid <- expand.grid(z.all)
  colnames(newz.grid) <- colnames(Z.MI)[ord]
  newz.grid <- newz.grid[,colnames(Z.MI)]

  if (!is.null(min.plot.dist)) {
    mindists <- rep(NA,nrow(newz.grid))
    for (i in seq_along(mindists)) {
      pt <- as.numeric(newz.grid[i, colnames(Z.MI)[ord[1]]])
      dists <- fields::rdist(matrix(pt, nrow = 1), Z.MI[, colnames(Z.MI)[ord[1]]])
      mindists[i] <- min(dists)
    }
  }



  if(method=="exact") {
    postmean.temp <- matrix(NA, nrow=length(sel)*K, ncol=ngrid)
    postvar.temp  <- array(NA, dim = c(ngrid,ngrid,length(sel)*K))

    for(k in 1:K){
      fit <- BKMRfits[[k]]
      y <- fit$y
      Z <- fit$Z
      X <- fit$X

      preds <- ComputePostmeanHnew.exact.MI(fit = fit, y = y, Z = Z, X = X, Znew = newz.grid, sel = sel)

      postmean.temp[((k-1)*length(sel)+1):(length(sel)*k),] <- preds$postmean_mat
      postvar.temp[,,((k-1)*length(sel)+1):(length(sel)*k)] <- preds$postvar_arr

    }

    ve <- var(postmean.temp)
    ev <- apply(postvar.temp, c(1, 2), mean)
    v  <- ve + ev

    preds.plot <- colMeans(postmean.temp)
    se.plot <- sqrt(diag(v))

  } else if(method=="approx") {

    postmean.temp <- matrix(NA, nrow=K, ncol=ngrid)
    postvar.temp  <- matrix(NA, nrow=K, ncol=ngrid)

    for(k in 1:K){
      fit <- BKMRfits[[k]]
      y <- fit$y
      Z <- fit$Z
      X <- fit$X

      preds <- ComputePostmeanHnew.approx(fit = fit, y = y, Z = Z, X = X, Znew = newz.grid, sel = sel)

      postmean.temp[k,] <- preds$postmean
      postvar.temp[k,]  <- diag(preds$postvar)

    }

    temp <- sapply(1:ngrid, function(x){Rubin.MI(mean.vec = postmean.temp[,x], variance.vec = postvar.temp[,x])})

    preds.plot <- temp["est",]
    se.plot    <- temp["sd",]

  } else stop("method must be one of c('approx', 'exact')")


  if(center) preds.plot <- preds.plot - mean(preds.plot)
  if(!is.null(min.plot.dist)) {
    preds.plot[mindists > min.plot.dist] <- NA
    se.plot[mindists > min.plot.dist] <- NA
  }

  res <- dplyr::data_frame(z = z1, est = preds.plot, se = se.plot)
}

