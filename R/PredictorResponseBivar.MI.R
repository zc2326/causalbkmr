#' Predict the exposure-response function at a new grid of points for MI BKMR
#'
#' @param BKMRfits  A list of multiple BKMR fits and that each of these fits were ran for the same number of MCMC iterations.
#' @param z.pairs data frame showing which pairs of predictors to plot
#' @param method  method for obtaining posterior summaries at a vector of new points. Options are "approx" and "exact"; defaults to "approx", which is faster particularly for large datasets
#' @param ngrid  number of grid points in each dimension
#' @param q.fixed  vector of quantiles at which to fix the remaining predictors in Z
#' @param sel logical expression indicating samples to keep; defaults to keeping the second half of all samples
#' @param min.plot.dist specifies a minimum distance that a new grid point needs to be from an observed data point in order to compute the prediction; points further than this will not be computed
#' @param center flag for whether to scale the exposure-response function to have mean zero
#' @param z.names optional vector of names for the columns of \code{z}
#' @param verbose TRUE or FALSE: flag of whether to print intermediate output to the screen
#' @param ...  other arguments to pass on to the prediction function
#'
#' @return a long data frame with the name of the first predictor, the name of the second predictor, the value of the first predictor, the value of the second predictor, the posterior mean estimate, and the posterior standard deviation of the estimated exposure response function
#' @export
#'
#' @importFrom dplyr select
#' @importFrom magrittr %<>%
PredictorResponseBivar.MI <- function(BKMRfits, z.pairs = NULL, method = "approx", ngrid = 50, q.fixed = 0.5, sel = NULL, min.plot.dist = 0.5, center = TRUE, z.names = colnames(BKMRfits[[1]]$Z), verbose = TRUE, ...) {

  start.time <- proc.time()["elapsed"]
  Z.MI <- Z.complete.MI(BKMRfits)
  l <- ncol(Z.MI)
  z.names <- colnames(Z.MI)

  K <- length(BKMRfits)

  npairs <- length(z.pairs)
  if(is.null(z.pairs)) npairs <- factorial(l)/factorial(l-2)/factorial(2)
  est.matrix <- matrix(NA,nrow=ngrid*ngrid*npairs, ncol=K)
  for(k in 1:K){
    fit <- BKMRfits[[k]]
    y <- fit$y
    Z <- fit$Z
    X <- fit$X

    temp <- PredictorResponseBivar.singfit.MI(fit=fit, z.pairs = z.pairs, method = method, ngrid = ngrid, q.fixed = q.fixed, sel = sel, min.plot.dist = min.plot.dist, center = center, z.names = z.names, verbose = TRUE, Z.MI=Z.MI,k=k,K=K, ...)

    est.matrix[,k] <- temp %>% select(est) %>% unlist(use.names=FALSE)

    end.time <- proc.time()["elapsed"]
    message(paste("MI fit", k, "out of", K, "complete: ", round((end.time - start.time)/60, digits=2), "min run time" ))
  }

  est.MI <- apply(est.matrix,1,mean)


  data.toreturn <- temp ### okay since the grid numbers and variable order are the same
  data.toreturn[,"est"] <- est.MI

  data.toreturn
}


PredictorResponseBivar.singfit.MI <- function(fit, y = NULL, Z = NULL, X = NULL, z.pairs = NULL, method = "approx", ngrid = 50, q.fixed = 0.5, sel = NULL, min.plot.dist = 0.5, center = TRUE, z.names = colnames(Z), verbose = TRUE, Z.MI,k=1,K=1, ...) {

  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
  }

  if (is.null(z.names)) {
    z.names <- colnames(Z.MI)
    if (is.null(z.names)) {
      z.names <- paste0("z", 1:ncol(Z))
    }
  }

  if (is.null(z.pairs)) {
    z.pairs <- expand.grid(z1 = 1:ncol(Z), z2 = 1:ncol(Z))
    z.pairs <- z.pairs[z.pairs$z1 < z.pairs$z2, ]
  }

  df <- dplyr::data_frame()
  for(i in 1:nrow(z.pairs)) {
    compute <- TRUE
    whichz1 <- z.pairs[i, 1] %>% unlist %>% unname
    whichz2 <- z.pairs[i, 2] %>% unlist %>% unname
    if(whichz1 == whichz2) compute <- FALSE
    z.name1 <- z.names[whichz1]
    z.name2 <- z.names[whichz2]
    names.pair <- c(z.name1, z.name2)
    if(nrow(df) > 0) { ## determine whether the current pair of variables has already been done
      completed.pairs <- df %>%
        dplyr::select_('variable1', 'variable2') %>%
        dplyr::distinct() %>%
        dplyr::transmute(z.pair = paste('variable1', 'variable2', sep = ":")) %>%
        unlist %>% unname
      if(paste(names.pair, collapse = ":") %in% completed.pairs | paste(rev(names.pair), collapse = ":") %in% completed.pairs) compute <- FALSE
    }
    if(compute) {
      if(verbose) message("MI fit ", k," out of ",K, ":  Pair ", i, " out of ", nrow(z.pairs))
      res <- PredictorResponseBivarPair.MI(fit = fit, y = y, Z = Z, X = X, whichz1 = whichz1, whichz2 = whichz2, method = method, ngrid = ngrid, q.fixed = q.fixed, sel = sel, min.plot.dist = min.plot.dist, center = center, z.names = z.names, Z.MI=Z.MI, ...)
      df0 <- res
      df0$variable1 <- z.name1
      df0$variable2 <- z.name2
      df0 %<>%
        dplyr::select_(~variable1, ~variable2, ~z1, ~z2, ~est, ~se)
      df <- dplyr::bind_rows(df, df0)
    }
  }
  df$variable1 <- factor(df$variable1, levels = z.names)
  df$variable2 <- factor(df$variable2, levels = z.names)
  df
}
