
#' Estimate controlled direct effect for BKMR
#'
#' @param a exposure variables at current level
#' @param astar exposure variables at counterfactual level
#' @param e.y effect modifier for the outcome variable
#' @param m.value values that the mediator is set to
#' @param m.quant values of the quantile that the mediator is set to
#' @param fit.y model fit regressing outcome on exposures, effect modifiers, mediator and confounders on outcome
#' @param alpha 1-confidence interval
#' @param sel a vector selecting which iterations of thfe fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#' @return Controlled direct effect for BKMR
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
#'
#' CDE.ey10 <- CDE.bkmr(a=a, astar=astar, e.y = e.y10, m.quant=c(0.1,0.5,0.75),
#'                      fit.y=fit.y, alpha=0.05, sel=sel, seed=777)
#' CDE.ey10$est
#'
#' plotdf <- as.data.frame(CDE.ey10$est)
#' plotdf["Effect"] <- rownames(plotdf)
#' ggplot(plotdf, aes(Effect, mean, ymin = lower, ymax = upper ))  +
#'   geom_pointrange(position = position_dodge(width = 0.75))  +  coord_flip()
#'
#'
#' }
#'
#'
#' @export

CDE.bkmr <- function(a, astar, e.y, m.value=NULL, m.quant=c(0.1,0.5,0.75), fit.y, alpha=0.05, sel, seed){
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
