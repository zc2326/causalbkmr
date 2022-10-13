
#' Estimate TE for BKMR-CMA-MI
#'
#' Get Ya and Yastar used in effects caculation
#'
#' @param a exposure variables at current level
#' @param astar exposure variables at counterfactual level
#' @param e.y effect modifier for the outcome variable
#' @param BKMRfits.y.TE A list of BKMR models for total effect model fit regressing outcome on exposures, effect modifiers and confounders on outcome, using multiple imputed data.
#' @param X.predict.Y counfounders for outcome
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#' @return A list containing the sample prediction for Ya and Yastar
#'
#' @details
#' For guided examples, go to
#' https://zc2326.github.io/causalbkmr/articles/BKMRCMA_QuickStart.html
#'
#' @examples
#' \dontrun{
#' library(causalbkmr)
#' dat <-  cma_sampledata(N=300, L=3, P=3, scenario=1, seed=7)
#'
#' A <- cbind(dat$z1, dat$z2, dat$z3)
#' X <- cbind(dat$x3)
#' y  <- dat$y
#' m  <- dat$M
#'
#' E.M <- NULL
#' E.Y <- dat$x2
#'
#' Z.M <- cbind(A,E.M)
#' Z.Y <- cbind(A, E.Y)
#' Zm.Y <- cbind(Z.Y, m)
#'
#' set.seed(1)
#' fit.y <- kmbayes(y=y, Z=Zm.Y, X=X, iter=5000, verbose=TRUE, varsel=FALSE)
#' #save(fit.y,file="bkmr_y.RData")
#'
#' set.seed(2)
#' fit.y.TE <- kmbayes(y=y, Z=Z.Y, X=X, iter=5000, verbose=TRUE, varsel=FALSE)
#' #save(fit.y.TE,file="bkmr_y_TE.RData")
#'
#' set.seed(3)
#' fit.m <- kmbayes(y=m, Z=Z.M, X=X, iter=5000, verbose=TRUE, varsel=FALSE)
#' #save(fit.m,file="bkmr_m.RData")
#'
#' X.predict <- matrix(colMeans(X),nrow=1)
#' astar <- c(apply(A, 2, quantile, probs=0.25))
#' a <- c(apply(A, 2, quantile, probs=0.75))
#'
#' e.y10 = quantile(E.Y, probs=0.1)
#' e.y90 = quantile(E.Y, probs=0.9)
#'
#' sel<-seq(2500,5000,by=5)
#'
#' TE.ey10 <- TE.bkmr(a=a, astar=astar, e.y = e.y10, fit.y.TE=fit.y.TE, X.predict=X.predict, alpha=0.05, sel=sel, seed=122)
#' TE.ey90 <- TE.bkmr(a=a, astar=astar, e.y = e.y90, fit.y.TE=fit.y.TE, X.predict=X.predict, alpha=0.05, sel=sel, seed=122)
#' TE.ey10$est
#' TE.ey90$est
#' }
#'
#' @export

YaYastar.SamplePred <- function(a, astar, e.y, BKMRfits.y.TE, X.predict.Y, sel, seed){
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
TE.bkmr.MI <- function(a, astar, e.y, fit.y.TE, X.predict.Y, alpha=0.05, sel, seed){

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
