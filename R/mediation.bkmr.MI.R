#' Estimate NDE/NIE for BKMR-CMA-MI
#'
#' @param a exposure variables at current level
#' @param astar exposure variables at counterfactual level
#' @param e.m effect modifier for the mediator variable
#' @param e.y effect modifier for the outcome variable
#' @param BKMRfits.m   A list of model fits using multiple imputed data regressing mediator on exposures and confounders on mediator
#' @param BKMRfits.y A list of model fits using multiple imputed data regressing outcome on exposures, effect modifiers, mediator and confounders on outcome
#' @param BKMRfits.y.TE A list of total effect model fits using multiple imputed data  regressing outcome on exposures, effect modifiers and confounders on outcome
#' @param X.predict.M counfounders for mediator
#' @param X.predict.Y counfounders for outcome
#' @param effects type(s) of effects that users want to output
#' @param m.value values that the mediator is set to
#' @param m.quant values of the quantile that the mediator is set to
#' @param alpha 1-confidence interval
#' @param sel a vector selecting which iterations of the fit should be retained or inference
#' @param seed the random seed to use to evaluate the code
#' @param K number of samples to generate for each MCMC iteration in YaMastar calculation
#'
#' @return A list contaning the sample prediction for TE, NDE, NIE and their summary statistics
#'
#' @details
#' For guided examples, go to
#' https://zc2326.github.io/causalbkmr/articles/BKMRCMA_QuickStart.html
#'
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
#' #' mediationeffects.ey10  <- mediation.bkmr(a=a, astar=astar, e.y = e.y10, fit.m=fit.m, fit.y=fit.y, fit.y.TE=fit.y.TE, X.predict.M=X.predict, X.predict.Y=X.predict, alpha=0.05, sel=sel, seed=22, K=10)
#' mediationeffects.ey10$est
#' plotdf <- as.data.frame(mediationeffects.ey10$est)
#' plotdf["Effect"] <- rownames(plotdf)
#' ggplot(plotdf, aes(Effect, mean, ymin = lower, ymax = upper ))  +
#'   geom_pointrange(position = position_dodge(width = 0.75))  +  coord_flip()
#'
#' }
#'
#' @export
mediation.bkmr.MI <- function(a, astar, e.m = NULL, e.y, BKMRfits.m=NULL, BKMRfits.y=NULL,
                              BKMRfits.y.TE=NULL,
                           X.predict.M=NULL, X.predict.Y=NULL,
                           effects = "all",  # c("all", "TE", "CDE", (meidation: "PNDE", "TNIE"), ("TNDE", "PNIE"))
                           m.quant=c(0.1,0.5,0.75),
                           m.value=NULL,
                           alpha = 0.05, sel, seed, K){

  if (sum(!effects %in% c("all", "TE", "CDE", "NDE", "NIE"))) {
    stop("effects must be in c('all', 'TE', 'CDE', 'NDE', 'NIE')")
  }
  else{
    if ("all" %in% effects){
      effects = c("TE", "CDE", "NDE", "NIE")
    }
    if ("TE" %in% effects){
      if (is.null(fit.y.TE)){
        stop("Must specify 'fit.y.TE'")
      }
      if (is.null(X.predict.Y)){
        stop("Must specify 'X.predict.Y'")
      }
    }
    if (sum(c("NIE", "NDE") %in% effects)){
      if (is.null(fit.y.TE) | is.null(fit.m) | is.null(fit.y)){
        stop("Must specify all three fits: 'fit.y.TE', 'fit.y', 'fit.m'")
      }
      if (is.null(X.predict.Y)){
        stop("Must specify 'X.predict.Y'")
      }
      if (is.null(X.predict.M)){
        stop("Must specify 'X.predict.M'")
      }
    }
    if ("CDE" %in% effects){
      if (is.null(fit.y)){
        stop("Must specify 'fit.y'")
      }
      # if (is.null(m.value) & is.null(m.quant)){
      #   stop("Must specify either 'm.value' or 'm.quant'")
      # }
      if (!is.null(m.value) & !is.null(m.quant)){
        m.quant = NULL # if both m.value and m.quant are specified, default set to m.value
      }
    }

    # start code
    TE = NULL; CDE = NULL; NDE = NULL; NIE = NULL;
    toreturn <- list()
    effects.temp = effects
    if ("CDE" %in% effects){
      effects.temp = effects[effects!="CDE"]
    }
    toreturn$est <- matrix(NA, nrow=length(effects.temp), ncol=5, dimnames=list(effects.temp, c("mean","median","lower","upper","sd")))

    if ("TE" %in% effects){
      TE <- TE.bkmr(a=a, astar=astar, e.y=e.y, fit.y.TE=fit.y.TE, X.predict.Y=X.predict.Y, alpha=alpha, sel=sel, seed=seed)
      toreturn$TE.samp <- TE$TE.samp
      toreturn$est["TE",] <- postresults(TE$TE.samp, alpha=alpha) [c("mean","median","lower","upper","sd")]
    }
    if (sum(c("NIE", "NDE") %in% effects)){
      if (is.null(TE)){
        TE <- TE.bkmr(a=a, astar=astar, e.y=e.y, fit.y.TE=fit.y.TE, X.predict.Y=X.predict.Y, alpha=alpha, sel=sel, seed=seed)
      }

      Ya     <- TE$Ya.samp
      Yastar <- TE$Yastar.samp

      YaMastar <- YaMastar.SamplePred(a=a, astar=astar, e.m = e.m, e.y = e.y, fit.m=fit.m, fit.y=fit.y,
                                      X.predict.M=X.predict.M, X.predict.Y=X.predict.Y, sel=sel, seed=seed, K=K)
      if ("NDE" %in% effects){
        NDE <- YaMastar - Yastar
        toreturn$NDE.samp <- NDE
        toreturn$est["NDE",] <- postresults(NDE, alpha=alpha) [c("mean","median","lower","upper","sd")]
      }
      if ("NIE" %in% effects){
        NIE <- Ya - YaMastar
        toreturn$NIE.samp <- NIE
        toreturn$est["NIE",] <- postresults(NIE, alpha=alpha) [c("mean","median","lower","upper","sd")]
      }
    }
    if ("CDE" %in% effects){
      CDE <- CDE.bkmr(a, astar, e.y, m.value=m.value, m.quant=m.quant, fit.y, alpha, sel, seed)
      toreturn$est = rbind(toreturn$est, CDE$est)
      toreturn = append(toreturn, CDE[2:length(CDE)])
    }
    return (toreturn)
  }
}
