#' Posterior/Bootstrap Summary Function
#'
#' Calculate mean/median/CI/sd of a vector of posterior/bootstrap samples
#'
#' @param posteriorsamp sample prediction
#' @param alpha 1-confidence interval
#' @return Mean. median , CI and sd of the sample prediction
#' @export

postresults <- function(posteriorsamp, alpha){
  toreturn <- vector()
  toreturn["mean"]   <- mean(posteriorsamp, na.rm=TRUE)
  toreturn["sd"]     <- stats::sd(posteriorsamp, na.rm=TRUE)
  toreturn["lower"]  <- stats::quantile(posteriorsamp, probs=alpha/2, na.rm=TRUE)
  toreturn["median"]  <- stats::quantile(posteriorsamp, probs=0.5, na.rm=TRUE)
  toreturn["upper"] <- stats::quantile(posteriorsamp, probs=1-alpha/2, na.rm=TRUE)
  return(toreturn)
}
