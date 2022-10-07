
#' YaMastar counterfactual
#'
#' Get YaMastar used in NDE and NIE caculation
#'
#' @param a exposure variables at current level
#' @param astar exposure variables at counterfactual level
#' @param e.m effect modifier for the mediator variable
#' @param e.y effect modifier for the outcome variable
#' @param fit.m model fit regressing mediator on exposures and confounders on mediator
#' @param fit.y model fit regressing outcome on exposures, effect modifiers, mediator and confounders on outcome
#' @param X.predict.M counfounders for mediator
#' @param X.predict.Y counfounders for outcome
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#' @param K number of samples to generate for each MCMC iteration
#' @return A vector containing the sample prediction for YaMastar
#' @export
YaMastar.SamplePred <- function(a, astar, e.m, e.y, fit.m, fit.y, X.predict.M, X.predict.Y, sel, seed, K){
  start.time <- proc.time()

  set.seed(seed)
  z.m <- c(a, e.m)
  EM.samp <- bkmr::SamplePred(fit.m, Znew = z.m, Xnew = X.predict.M, sel=sel)
  Mastar     <- as.vector(EM.samp)

  sigma.samp  <- sqrt(fit.m$sigsq.eps[sel])
  random.samp <- matrix(stats::rnorm(length(sel)*K),nrow=length(sel),ncol=K)

  Mastar.samp  <- Mastar + sigma.samp*random.samp

  YaMastar.samp.mat     <- matrix(NA,nrow=length(sel),ncol=K)
  z.y = c(a, e.y)
  for(j in 1:length(sel)){
    Mastar.j <-  Mastar.samp[j,]
    aMastar.j <- cbind(matrix(z.y, nrow=K, ncol=length(z.y), byrow=TRUE), Mastar.j)
    YaMastar.j <- bkmr::SamplePred(fit.y, Znew = aMastar.j, Xnew = X.predict.Y, sel=sel[j])
    YaMastar.samp.mat[j,] <- as.vector(YaMastar.j)

    end.time.temp <- proc.time()
    if(j%%50==0) print(paste("iter", j, "time: ", round((end.time.temp - start.time)["elapsed"]/60,2),"min"))
  }
  toreturn <- apply(YaMastar.samp.mat,1,mean)

  return(toreturn)
}
