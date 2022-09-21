

#' Compute the posterior mean and variance of \code{h} at a new predictor values using approx method for MI BKMR fits
#'
#'
#' Compute the posterior mean and variance of \code{h} at a new predictor values
#'
#' @param fit An object contatin the results return by the kmbayes function
#' @param y  a vector of outcome data of length \code{n}.
#' @param Z an \code{n-by-M} matrix of predictor variables to be included in the \code{h} function. Each row represents an observation and each column represents an predictor.
#' @param X an \code{n-by-K} matrix of covariate data where each row represents an observation and each column represents a covariate. Should not contain an intercept column.
#' @param Znew matrix of new predictor values at which to predict new h, where each row represents a new observation. If set to NULL then will default to using the observed exposures Z
#' @param sel selects which iterations of the MCMC sampler to use for inference
#'
#' @return A list of mean, variance, entire mean matrix and variance array
#' @export
#' @importFrom stats var
#' @importFrom bkmr ExtractEsts
#' @examples
ComputePostmeanHnew.approx <- function (fit, y = NULL, Z = NULL, X = NULL, Znew = NULL, sel = NULL) {
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y))
      y <- fit$y
    if (is.null(Z))
      Z <- fit$Z
    if (is.null(X))
      X <- fit$X
  }
  if (!is.null(Znew)) {
    if (is.null(dim(Znew))){
      Znew <- matrix(Znew, nrow = 1)
    }
    else{
      Znew <- data.matrix(Znew)
    }
  }
  if (is.null(dim(X)))
    X <- matrix(X, ncol = 1)
  ests <- ExtractEsts(fit, sel = sel)
  sigsq.eps <- ests$sigsq.eps[, "mean"]
  r <- ests$r[, "mean"]
  beta <- ests$beta[, "mean"]
  lambda <- ests$lambda[, "mean"]
  if (fit$family == "gaussian") {
    ycont <- y
  }
  else if (fit$family == "binomial") {
    ycont <- ests$ystar[, "mean"]
  }
  Kpart <- makeKpart(r, Z)
  K <- exp(-Kpart)
  V <- diag(1, nrow(Z), nrow(Z)) + lambda[1] * K
  cholV <- chol(V)
  Vinv <- chol2inv(cholV)
  if (!is.null(Znew)) {
    n0 <- nrow(Z)
    n1 <- nrow(Znew)
    nall <- n0 + n1
    Kpartall <- makeKpart(r, rbind(Z, Znew))
    Kmat <- exp(-Kpartall)
    Kmat0 <- Kmat[1:n0, 1:n0, drop = FALSE]
    Kmat1 <- Kmat[(n0 + 1):nall, (n0 + 1):nall, drop = FALSE]
    Kmat10 <- Kmat[(n0 + 1):nall, 1:n0, drop = FALSE]
    lamK10Vinv <- lambda[1] * Kmat10 %*% Vinv
    postvar <- lambda[1] * sigsq.eps * (Kmat1 - lamK10Vinv %*%
                                          t(Kmat10))
    postmean <- lamK10Vinv %*% (ycont - X %*% beta)
  }
  else {
    lamKVinv <- lambda[1] * K %*% Vinv
    postvar <- lambda[1] * sigsq.eps * (K - lamKVinv %*%
                                          K)
    postmean <- lamKVinv %*% (ycont - X %*% beta)
  }
  ret <- list(postmean = drop(postmean), postvar = postvar)
  ret
}

