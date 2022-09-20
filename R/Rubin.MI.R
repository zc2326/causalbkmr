

#' Rubins 1987 method
#'
#' Function to implement Rubins 1987 method to estimate point estimates and std errors when combining information across MI fits (for approx method)
#'
#' @param mean.vec A vector consist of mean values
#' @param variance.vec A vector consist of variance values
#'
#' @return Point estimates and std errors when combining information across MI fits (for approx method)
#' @export
#'
#' @examples
Rubin.MI <- function(mean.vec, variance.vec){
  K    <- length(mean.vec)
  qbar <- mean(mean.vec)
  wbar <- mean(variance.vec)
  B    <- var(mean.vec)
  var  <- wbar+(1+1/K)*B
  c(est=qbar, sd=sqrt(var))
}
