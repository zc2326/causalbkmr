
#' CMA Sample Data Simulation
#'
#' Simulate simple sample dataset for BKMRCMA illustration, function to generate simulated data with multiple exposures and true mediated effect
#'
#' @param N sample size of dataset
#' @param L number of mixture components
#' @param P number of covariates
#' @param scenario Scenarios  that are similar to the scenarios presented in the paper
#' @param seed the random seed to use to evaluate the code
#'
#'
#' @return a data frame
#'
#' @example
#' dat <-  cma_sampledata(N=300, L=3, P=3, scenario=1, seed=7)
#' head(dat$data, n = 3L)
#' @export
#'
#' @importFrom stats quantile plogis rnorm
#' @importFrom SimDesign rmvnorm
cma_sampledata <- function(N, L, P, scenario, seed) {
  set.seed(seed)

  ## assume covariance of 0.3 for all L elements
  covmat <- matrix(0.3, nrow=L, ncol=L)
  diag(covmat) <- 1

  ## first generate the exposures from a MVN
  dat <- list(Z = rmvnorm(N, mean=rep(0,L), sigma = covmat),
              X = rmvnorm(N, mean = rep(0, P)), sigma = diag(1, P, P))
  colnames(dat$Z) <- paste("z",1:L,sep="")
  colnames(dat$X) <- paste("x",1:L,sep="")

  ## scale signal-to-noise ratio for the given sample size
  signoiseM <- 0.36 * sqrt(380/N)
  signoiseY <- 0.42 * sqrt(380/N)

  ## assign the same dose-response surface to generate M and Y for the four scenarios in paper
  if(scenario==1){
    hfunM <- function(z, ind1 = 1, ind2 = 2, ind3 = 4) 1/4*z[ind1] + 1/4*z[ind3] + 1/4*z[ind1]*z[ind3]
    hfunY <- function(z, ind1 = 1, ind2 = 4, ind3 = 5) 1/12*((z[ind1]+3) + (z[ind2]+3) + 1/3*(z[ind1]+3)*(z[ind2]+3)-12)+ 1/4*z[ind1]*z[ind3]
    ## scale random noise based on h functions
    sig.trueM <- 0.05/signoiseM
    sig.trueY <- 0.044/signoiseY
  }else if(scenario==2){
    hfunM <- function(z, ind1 = 1, ind2 = 2, ind3 = 4) (plogis(z[ind1], 0, 0.2)-0.5)+ 1/4*z[ind3]
    hfunY <- function(z, ind1 = 1, ind2 = 2, ind3 = 5) (plogis(1/6*((z[ind1]+3) + (z[ind2]+3) + 1/3*(z[ind1]+3)*(z[ind2]+3)), 2, 0.2)-0.5)+ 1/4*z[ind3]
    ## scale random noise based on h functions
    sig.trueM <- 0.4229/signoiseM
    sig.trueY <- 0.2405/signoiseY
  }else if(scenario==3){
    hfunM <- function(z, ind1 = 1, ind2 = 2, ind3 = 4) -(1/8)*(z[ind1])^2+0.5 + 1/4*z[ind3]
    hfunY <- function(z, ind1 = 1, ind2 = 2, ind3 = 5) (plogis(1/6*((z[ind1]+3) + (z[ind2]+3) + 1/3*(z[ind1]+3)*(z[ind2]+3)), 2, 0.2)-0.5)+ 1/4*z[ind3]
    ## scale random noise based on h functions
    sig.trueM <- 0.1925/signoiseM
    sig.trueY <- 0.2405/signoiseY
  }else if(scenario==4){
    hfunM <- hfunY <- function(z, ind1 = 1, ind2 = 2, ind3 = 4) plogis(-(1/4)*(z[ind1]-1)^2-1/4*(z[ind2]-1)^2+1/3*(z[ind1]+3)*(z[ind2]+3), 2, 0.5)-0.5+ 1/4*z[ind3]
    ## scale random noise based on h functions
    sig.trueM <- 0.4095/signoiseM
    sig.trueY <- 0.3702/signoiseY
  }else{
    stop("scenario must be one of the following: 1,2,3,4")
  }

  dat$epsM <- rnorm(N, sd=sig.trueM)
  dat$epsY <- rnorm(N, sd=sig.trueY)
  dat$hM   <- apply(cbind(dat$Z, dat$X), 1, hfunM)
  dat$M    <- with(dat, hM + epsM)
  dat$Zm   <- cbind(dat$Z,dat$M, dat$X)
  dat$hY   <- apply(dat$Zm, 1, hfunY)
  dat$y    <- with(dat, hY + epsY)

  dat$data <- data.frame(with(dat,cbind(Z, M, X, y)))
  return(dat)
}
