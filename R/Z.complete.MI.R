
#' Function to compile the exposures from all MI datasets
#'
#' @param BKMRfits @param BKMRfits  A list of multiple BKMR fits and that each of these fits were ran for the same number of MCMC iterations.
#'
#' @return a dataframe
#'
#' @details
#' This step is needed to so that the contrast used when estimating the effects across the MI BKMR fits is consistent. (i.e. the same 25th, 50th, and 75th percentile of the metals is used when estimating the effect of a change in all metals at their 25th to all their 75th on the outcome for each MI BKMR fit).
#' If the data are not imputed in the exposures, the original Z matrix can used as Z.complete.MI

#' For guided examples, go to https://zc2326.github.io/causalbkmr/articles/MI_BKMR.html
#' @examples
#' Z.MI <- Z.complete.MI(BKMRfits10)
#'
#'
#'
#' @export
Z.complete.MI <- function(BKMRfits){
  n <- nrow(BKMRfits[[1]]$Z)
  l <- ncol(BKMRfits[[1]]$Z)
  K <- length(BKMRfits)

  Z.full <- matrix(NA, nrow=n*K, ncol=l)
  ifelse(is.null(colnames(BKMRfits[[1]]$Z)), colnames(Z.full) <- paste0("z", 1:l), colnames(Z.full) <- colnames(BKMRfits[[1]]$Z))

  for(k in 1:K){
    fit <- BKMRfits[[k]]
    Z.full[(n*(k-1)+1):(n*k),] <- fit$Z
  }
  Z.full
}
