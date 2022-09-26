

#' makeKpart
#'
#' @param Z1 need to update
#' @param Z2 need to update
#' @param r need to update
#' @importFrom fields rdist
#' @return need to update
#' @export
#'
#'
#'
makeKpart <- function(r, Z1, Z2 = NULL) {
Z1r <- sweep(Z1, 2, sqrt(r), "*")
if (is.null(Z2)) {
  Z2r <- Z1r
} else {
  Z2r <- sweep(Z2, 2, sqrt(r), "*")
}
Kpart <- fields::rdist(Z1r, Z2r)^2
Kpart
}
#' Title
#'
#' @param r need to update
#' @param lambda need to update
#' @param Z need to update
#' @param data.comps need to update
#'
#' @return need to update
#' @export
#'
#' @examples
makeVcomps <- function(r, lambda, Z, data.comps) {
  if (is.null(data.comps$knots)) {
    Kpart <- makeKpart(r, Z)
    V <- diag(1, nrow(Z), nrow(Z)) + lambda[1]*exp(-Kpart)
    if (data.comps$nlambda == 2) {
      V <- V + lambda[2]*data.comps$crossTT
    }
    cholV <- chol(V)
    Vinv <- chol2inv(cholV)
    logdetVinv <- -2*sum(log(diag(cholV)))
    Vcomps <- list(Vinv = Vinv, logdetVinv = logdetVinv)
  } else {## predictive process approach
    ## note: currently does not work with random intercept model
    nugget <- 0.001
    n0 <- nrow(Z)
    n1 <- nrow(data.comps$knots)
    nall <- n0 + n1
    # Kpartall <- makeKpart(r, rbind(Z, data.comps$knots))
    # Kall <- exp(-Kpartall)
    # K0 <- Kall[1:n0, 1:n0 ,drop=FALSE]
    # K1 <- Kall[(n0+1):nall, (n0+1):nall ,drop=FALSE]
    # K10 <- Kall[(n0+1):nall, 1:n0 ,drop=FALSE]
    K1 <- exp(-makeKpart(r, data.comps$knots))
    K10 <- exp(-makeKpart(r, data.comps$knots, Z))
    Q <- K1 + diag(nugget, n1, n1)
    R <- Q + lambda[1]*tcrossprod(K10)
    cholQ <- chol(Q)
    cholR <- chol(R)
    Qinv <- chol2inv(cholQ)
    Rinv <- chol2inv(cholR)
    Vinv <- diag(1, n0, n0) - lambda[1]*t(K10) %*% Rinv %*% K10
    logdetVinv <- 2*sum(log(diag(cholQ))) - 2*sum(log(diag(cholR)))
    Vcomps <- list(Vinv = Vinv, logdetVinv = logdetVinv, cholR = cholR, Q = Q, K10 = K10, Qinv = Qinv, Rinv = Rinv)
  }
  Vcomps
}
