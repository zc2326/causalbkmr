#' Estimate controlled direct effect for BKMR-CMA-MI
#'
#' @param a exposure variables at current level
#' @param astar exposure variables at counterfactual level
#' @param e.y effect modifier for the outcome variable
#' @param m.value values that the mediator is set to
#' @param m.quant values of the quantile that the mediator is set to
#' @param BKMRfits.Y A list of `K` model fit regressing outcome on exposures, effect modifiers, mediator and confounders on outcome
#' @param alpha 1-confidence interval
#' @param sel a vector selecting which iterations of thfe fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#'
#' @return Controlled direct effect for BKMR
#'
#' @details
#' For guided examples, go to
#' https://zc2326.github.io/causalbkmr/articles/BKMRCMA_QuickStart.html
#'
#' @export

CDE.bkmr.MI <- function(a, astar, e.y, m.value=NULL, m.quant=c(0.1,0.5,0.75), BKMRfits.Y, alpha=0.05, sel, seed){
  if (!is.null(m.value) & !is.null(m.quant)){
    m.quant = NULL # if both m.value and m.quant are specified, default set to m.value
  }
  toreturn <- list()
  m <- fit.y$Z[,ncol(fit.y$Z)]  ### okay as long as m is the LAST variable in Zm birthlength
  Z <- fit.y$Z[,-ncol(fit.y$Z)] # exposure + effect modifier
  X.predict<-rep(0,ncol(fit.y$X)) # in the calculation for CDE this value doesn't matter since it will get cancelled out
  if (is.null(m.value)){
    toreturn$est <- matrix(NA, nrow=length(m.quant), ncol=5, dimnames=list(paste0("CDE",m.quant*100, "%"), c("mean","median","lower","upper","sd")))
    print(paste("Running", length(m.quant), "mediator values for CDE:"))
    for(i in seq_along(m.quant)){
      print(paste(i, "out of", length(m.quant)))
      mnew <- stats::quantile(m, probs=m.quant[i]) # mediator set to certain quantile

      set.seed(seed)
      z.y = c(a, e.y)
      zstar.y = c(astar, e.y)
      newZm  <- rbind(c(z.y,mnew),c(zstar.y,mnew)) # exposure(3), 10% age, 10% mediator
      CDE.mat <- bkmr::SamplePred(fit.y, Znew = newZm, Xnew = X.predict, sel=sel) # X.predict is mean of confounders

      Yam     <- CDE.mat[,"znew1"]
      Yastarm <- CDE.mat[,"znew2"]

      CDE <- as.vector(Yam - Yastarm)
      toreturn[[paste0("CDE",m.quant[i]*100,"%.samp")]] <- CDE
      toreturn$est[paste0("CDE",m.quant[i]*100,"%"),] <- postresults(CDE, alpha=alpha)[c("mean","median","lower","upper","sd")]
    }
  }
  else if (is.null(m.quant)){
    toreturn$est <- matrix(NA, nrow=length(m.value), ncol=5, dimnames=list(paste0("CDE",m.value), c("mean","median","lower","upper","sd")))
    print(paste("Running", length(m.value), "mediator values for CDE:"))
    for(i in seq_along(m.value)){
      print(paste(i, "out of", length(m.value)))
      mnew = m.value[i]

      set.seed(seed)
      z.y = c(a, e.y)
      zstar.y = c(astar, e.y)
      newZm  <- rbind(c(z.y,mnew),c(zstar.y,mnew)) # exposure(3), 10% age, 10% mediator
      CDE.mat <- bkmr::SamplePred(fit.y, Znew = newZm, Xnew = X.predict, sel=sel) # X.predict is mean of confounders

      Yam     <- CDE.mat[,"znew1"]
      Yastarm <- CDE.mat[,"znew2"]

      CDE <- as.vector(Yam - Yastarm)
      toreturn[[paste0("CDE",mnew,".samp")]] <- CDE
      toreturn$est[paste0("CDE",mnew),] <- postresults(CDE, alpha=alpha)[c("mean","median","lower","upper","sd")]
    }
  }

  return(toreturn)
}
