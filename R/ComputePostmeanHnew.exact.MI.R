#' Compute the posterior mean and variance of \code{h} at a new predictor values using exact method for MI BKMR fits
#'
#'
#' #' Compute the posterior mean and variance of \code{h} at a new predictor values
#' Function to estimate the posterior mean and variance by obtaining the posterior mean and
#' variance at #' particular iterations and then using the iterated mean and variance formulas.
#' This function returns the entire mean matrix and variance array
#' needed for each of the future functions to obtain an unbiased estimate
#' of the SE used to create CI in plots.
#'
#' @param fit An object contatinint the results return by the kmbayes function
#' @param y  a vector of outcome data of length \code{n}.
#' @param Z an \code{n-by-M} matrix of predictor variables to be included in the \code{h} function. Each row represents an observation and each column represents an predictor.
#' @param X an \code{n-by-K} matrix of covariate data where each row represents an observation and each column represents a covariate. Should not contain an intercept column.
#' @param Znew matrix of new predictor values at which to predict new h, where each row represents a new observation. If set to NULL then will default to using the observed exposures Z
#' @param sel selects which iterations of the MCMC sampler to use for inference
#'
#' @return A list of mean, variance, entire mean matrix and variance array
#' @export
#' @importFrom stats var
ComputePostmeanHnew.exact.MI <- function(fit, y = NULL, Z = NULL, X = NULL, Znew = NULL, sel = NULL) {

  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
  }

  if (!is.null(Znew)) {
    if (is.null(dim(Znew))) Znew <- matrix(Znew, nrow = 1)
    if (class(Znew) == "data.frame") Znew <- data.matrix(Znew)
    if (ncol(Z) != ncol(Znew)) {
      stop("Znew must have the same number of columns as Z")
    }
  }

  if (is.null(dim(X))) X <- matrix(X, ncol=1)

  # if (!is.null(fit$Vinv)) {
  #   sel <- attr(fit$Vinv, "sel")
  # }

  if (is.null(sel)) {
    sel <- with(fit, seq(floor(iter/2) + 1, iter, 10))
    if (length(sel) < 100) {
      sel <- with(fit, seq(floor(iter/2) + 1, iter, length.out = 100))
    }
    sel <- unique(floor(sel))
  }

  family <- fit$family
  data.comps <- fit$data.comps
  post.comps.store <- list(postmean = vector("list", length(sel)),
                           postvar = vector("list", length(sel))
  )

  for (i in seq_along(sel)) {
    s <- sel[i]
    beta <- fit$beta[s, ]
    lambda <- fit$lambda[s, ]
    sigsq.eps <- fit$sigsq.eps[s]
    r <- fit$r[s, ]

    if (family == "gaussian") {
      ycont <- y
    } else if (family == "binomial") {
      ycont <- fit$ystar[s, ]
    }

    Kpart <- makeKpart(r, Z)
    K <- exp(-Kpart)
    Vcomps <- makeVcomps(r = r, lambda = lambda, Z = Z, data.comps = data.comps)
    Vinv <- Vcomps$Vinv
    # if (is.null(fit$Vinv)) {
    # V <- diag(1, nrow(Z), nrow(Z)) + lambda[1]*K
    # cholV <- chol(V)
    # Vinv <- chol2inv(cholV)
    # } else {
    #   Vinv <- fit$Vinv[[i]]
    # }

    if (!is.null(Znew)) {
      # if(is.null(data.comps$knots)) {
      n0 <- nrow(Z)
      n1 <- nrow(Znew)
      nall <- n0 + n1
      Kpartall <- makeKpart(r, rbind(Z, Znew))
      Kmat <- exp(-Kpartall)
      Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
      Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
      Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]

      lamK10Vinv <- lambda[1]*Kmat10 %*% Vinv
      postvar <- lambda[1]*sigsq.eps*(Kmat1 - lamK10Vinv %*% t(Kmat10))
      postmean <- lamK10Vinv %*% (ycont - X%*%beta)
      # } else {
      # stop("GPP not yet implemented")
      # }
    } else {
      lamKVinv <- lambda[1]*K%*%Vinv
      postvar <- lambda[1]*sigsq.eps*(K - lamKVinv%*%K)
      postmean <- lamKVinv %*% (ycont - X%*%beta)
    }

    post.comps.store$postmean[[i]] <- postmean
    post.comps.store$postvar[[i]] <- postvar

  }

  postmean_mat <- t(do.call("cbind", post.comps.store$postmean))
  m <- colMeans(postmean_mat)
  postvar_arr <- with(post.comps.store,
                      array(unlist(postvar),
                            dim = c(nrow(postvar[[1]]), ncol(postvar[[1]]), length(postvar)))
  )
  ve <- var(postmean_mat)
  ev <- apply(postvar_arr, c(1, 2), mean)
  v <- ve + ev
  ret <- list(postmean = m, postvar = v, postmean_mat = postmean_mat, postvar_arr = postvar_arr)

  ret
}
