
#' Estimate TE for BKMR
#'
#' Get Ya and Yastar used in effects caculation
#'
#' @param a exposure variables at current level
#' @param astar exposure variables at counterfactual level
#' @param e.y effect modifier for the outcome variable
#' @param fit.y.TE total effect model fit regressing outcome on exposures, effect modifiers and confounders on outcome
#' @param X.predict.Y counfounders for outcome
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#' @return A list containing the sample prediction for Ya and Yastar
#' @export

YaYastar.SamplePred <- function(a, astar, e.y, fit.y.TE, X.predict.Y, sel, seed){
  set.seed(seed)
  z.y = c(a, e.y)
  zstar.y = c(astar, e.y)
  newz <- rbind(z.y, zstar.y)

  # give prediction Y for both a and astar
  TE.mat <- bkmr::SamplePred(fit.y.TE, Znew = newz, Xnew = X.predict.Y, sel = sel)
  Ya <- TE.mat[,"znew1"]
  Yastar <- TE.mat[,"znew2"]

  return(list(Ya = Ya, Yastar = Yastar))
}


#' Estimate total effect for BKMR
#'
#' @param a exposure variables at current level
#' @param astar exposure variables at counterfactual level
#' @param e.y effect modifier for the outcome variable
#' @param fit.y.TE total effect model fit regressing outcome on exposures, effect modifiers and confounders on outcome
#' @param X.predict.Y counfounders for outcome
#' @param alpha 1-confidence interval
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#' @return Totak effect for BKMR
#' @export
TE.bkmr <- function(a, astar, e.y, fit.y.TE, X.predict.Y, alpha=0.05, sel, seed){

  toreturn <- list()
  toreturn$est <- matrix(NA, nrow=1, ncol=5, dimnames=list("TE", c("mean","median","lower","upper","sd")))


  YaYastar <- YaYastar.SamplePred(a, astar, e.y, fit.y.TE, X.predict.Y, sel, seed)
  Ya     <- YaYastar$Ya
  Yastar <- YaYastar$Yastar

  toreturn$TE.samp     <- as.vector(Ya - Yastar)
  toreturn$Ya.samp     <- as.vector(Ya)
  toreturn$Yastar.samp <- as.vector(Yastar)

 # TE.sum <- postresults(toreturn$TE, alpha=alpha)
  toreturn$est[1, ]  <- postresults(toreturn$TE, alpha=alpha)[c("mean","median","lower","upper","sd")]
  return(toreturn)
}
