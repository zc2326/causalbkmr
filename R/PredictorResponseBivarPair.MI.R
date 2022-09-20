#' Plot bivariate predictor-response function on a new grid of points for MI BKMR
#'
#' @param fit An object containing the results returned by a the kmbayes function
#' @param y  a vector of outcome data of length \code{n}.
#' @param Z an \code{n-by-M} matrix of predictor variables to be included in the \code{h} function. Each row represents an observation and each column represents an predictor.
#' @param X an \code{n-by-K} matrix of covariate data where each row represents an observation and each column represents a covariate. Should not contain an intercept column.
#' @param whichz1 vector identifying the first predictor that (column of \code{Z}) should be plotted
#' @param whichz2 vector identifying the second predictor that (column of \code{Z}) should be plotted
#' @param whichz3 vector identifying the third predictor that will be set to a pre-specified fixed quantile (determined by \code{prob})
#' @param method method for obtaining posterior summaries at a vector of new points. Options are "approx" and "exact"; defaults to "approx", which is faster particularly for large datasets
#' @param prob pre-specified quantile to set the third predictor (determined by \code{whichz3}); defaults to 0.5 (50th percentile)
#' @param q.fixed vector of quantiles at which to fix the remaining predictors in \code{Z}
#' @param sel logical expression indicating samples to keep; defaults to keeping the second half of all samples
#' @param ngrid number of grid points to cover the range of each predictor (column in \code{Z})
#' @param min.plot.dist specifies a minimum distance that a new grid point needs to be from an observed data point in order to compute the prediction; points further than this will not be computed
#' @param center flag for whether to scale the exposure-response function to have mean zero
#' @param Z.MI Multiple Imputed  \code{Z}
#' @param ... other arguments to pass on to the prediction function
#'
#' @importFrom bkmr ComputePostmeanHnew
#' @importFrom fields rdist
#' @return a data frame with value of the first predictor, the value of the second predictor, the posterior mean estimate, and the posterior standard deviation
#' @export
#'
#' @examples
PredictorResponseBivarPair.MI <- function(fit, y, Z, X, whichz1 = 1, whichz2 = 2, whichz3 = NULL, method = "approx", prob = 0.5, q.fixed = 0.5, sel = NULL, ngrid = 50, min.plot.dist = 0.5, center = TRUE, Z.MI, ...) {
  if(ncol(Z) < 3) stop("requires there to be at least 3 Z variables")

  if(is.null(colnames(Z))) colnames(Z) <- paste0("z", 1:ncol(Z))

  if(is.null(whichz3)) {
    ord <- c(whichz1, whichz2, setdiff(1:ncol(Z), c(whichz1, whichz2)))
  } else {
    ord <- c(whichz1, whichz2, whichz3, setdiff(1:ncol(Z), c(whichz1, whichz2, whichz3)))
  }
  z1 <- seq(min(Z.MI[,ord[1]]), max(Z.MI[,ord[1]]), length=ngrid)
  z2 <- seq(min(Z.MI[,ord[2]]), max(Z.MI[,ord[2]]), length=ngrid)
  z3 <- quantile(Z.MI[, ord[3]], probs = prob)
  z.all <- c(list(z1), list(z2), list(z3))
  if(ncol(Z) > 3) {
    z.others <- lapply(4:ncol(Z), function(x) quantile(Z.MI[,ord[x]], q.fixed))
    z.all <- c(z.all, z.others)
  }
  newz.grid <- expand.grid(z.all)
  z1save <- newz.grid[, 1]
  z2save <- newz.grid[, 2]
  colnames(newz.grid) <- colnames(Z)[ord]
  newz.grid <- newz.grid[,colnames(Z)]

  if(!is.null(min.plot.dist)) {
    mindists <- rep(NA, nrow(newz.grid))
    for(k in seq_along(mindists)) {
      pt <- as.numeric(newz.grid[k,c(colnames(Z)[ord[1]],colnames(Z)[ord[2]])])
      dists <- fields::rdist(matrix(pt, nrow = 1), Z[, c(colnames(Z)[ord[1]],colnames(Z)[ord[2]])])
      mindists[k] <- min(dists)
    }
  }

  if (method %in% c("approx", "exact")) {
    preds <- ComputePostmeanHnew(fit = fit, y = y, Z = Z, X = X, Znew = newz.grid, sel = sel, method = method)
    preds.plot <- preds$postmean
    se.plot <- sqrt(diag(preds$postvar))
  } else {
    stop("method must be one of c('approx', 'exact')")
  }
  if(center) preds.plot <- preds.plot - mean(preds.plot)
  if(!is.null(min.plot.dist)) {
    preds.plot[mindists > min.plot.dist] <- NA
    se.plot[mindists > min.plot.dist] <- NA
  }
  #     hgrid <- matrix(preds.plot, ngrid, ngrid, dimnames=list(z1=round(z1,2), z2=round(z2,2)))
  #     se.grid <- matrix(se.plot, ngrid, ngrid, dimnames=list(z1=round(z1,2), z2=round(z2,2)))

  res <- dplyr::data_frame(z1 = z1save, z2 = z2save, est = preds.plot, se = se.plot)
}
