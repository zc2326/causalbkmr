#' Calculate overall risk summaries for BKMR MI
#'
#' Compare estimated h function when all predictors are at a particular quantile to when all are at a second fixed quantile
#'
#' @param BKMRfits  A list of multiple BKMR fits and that each of these fits were ran for the same number of MCMC iterations.
#' @param qs  vector of quantiles at which to calculate the overall risk summary
#' @param q.fixed  a second quantile at which to compare the estimated {h} function
#' @param q.alwaysfixed  the quantile values in the point which we want to keep fixed for all comparisons
#' @param index.alwaysfixed the index values in the point which we want to keep fixed for all comparisons
#' @param sel  selects which iterations of the MCMC sampler to use for inference
#' @param method method for obtaining posterior summaries at a vector of new points. Options are"approx" and "exact"; defaults to "approx", which is faster particularly for large datasets
#'
#' @importFrom stats var
#' @return  a data frame containing the (posterior mean) estimate and posterior standard deviation of the overall risk measures
#' @export
OverallRiskSummaries.MI <- function(BKMRfits, qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, q.alwaysfixed = NULL, index.alwaysfixed = NULL, sel = NULL, method = "approx") {

  start.time <- proc.time()["elapsed"]
  cc <- c(-1, 1)
  K <- length(BKMRfits)

  Z.MI <- Z.complete.MI(BKMRfits)

  toreturn <- data.frame(quantile=qs, est=rep(NA,times=length(qs)), sd=rep(NA,times=length(qs)))

  if(method=="exact") {
    print("exact method")
    preds.fun <- function(znew) ComputePostmeanHnew.exact.MI(fit = fit, y = y, Z = Z, X = X, Znew = znew, sel = sel)
    for(i in 1:length(qs)){
      quant <- qs[i]

      ## 2 = nrow(newz)
      postmean.temp <- matrix(NA, nrow=length(sel)*K, ncol=2)
      postvar.temp  <- array(NA, dim = c(2,2,length(sel)*K))
      for(k in 1:K){
        fit <- BKMRfits[[k]]
        y <- fit$y
        Z <- fit$Z
        X <- fit$X

        point1 <- apply(Z.MI, 2, quantile, q.fixed)
        point2 <- apply(Z.MI, 2, quantile, quant)

        ## if both q.alwaysfixed and index.alwaysfixed are specified,
        ## change the values in the point which we want to keep fixed for all comparisons
        if(!is.null(q.alwaysfixed) & !is.null(index.alwaysfixed)){
          point1[index.alwaysfixed] <- point2[index.alwaysfixed] <- apply(Z.MI[,index.alwaysfixed, drop=FALSE],2,quantile, q.alwaysfixed)
        }

        newz <- rbind(point1, point2)

        preds <- preds.fun(newz)
        postmean.temp[((k-1)*length(sel)+1):(length(sel)*k),] <- preds$postmean_mat
        postvar.temp[,,((k-1)*length(sel)+1):(length(sel)*k)] <- preds$postvar_arr
      }
      m  <- colMeans(postmean.temp)
      ve <- var(postmean.temp)
      ev <- apply(postvar.temp, c(1, 2), mean)
      v  <- ve + ev

      toreturn[i,"est"] <- drop(cc %*% m)
      toreturn[i,"sd"]  <- drop(sqrt(cc %*% v %*% cc))

      end.time <- proc.time()["elapsed"]
      print(paste(i,"out of", length(qs), "complete: ", round((end.time - start.time)/60, digits=2), "min run time" ))
    }
  } else if(method=="approx") {
    print("approx method")
    preds.fun <- function(znew) causalbkmr::ComputePostmeanHnew.approx(fit = fit, y = y, Z = Z, X = X, Znew = znew, sel = sel)

    for(i in 1:length(qs)){
      quant <- qs[i]
      ## 2 = nrow(newz)

      est.vec <- rep(NA, times=K)
      var.vec <- rep(NA, times=K)
      for(k in 1:K){
        fit <- BKMRfits[[k]]
        y <- fit$y
        Z <- fit$Z
        X <- fit$X

        point1 <- apply(Z.MI, 2, quantile, q.fixed)
        point2 <- apply(Z.MI, 2, quantile, quant)

        ## if both q.alwaysfixed and index.alwaysfixed are specified,
        ## change the values in the point which we want to keep fixed for all comparisons
        if(!is.null(q.alwaysfixed) & !is.null(index.alwaysfixed)){
          point1[index.alwaysfixed] <- point2[index.alwaysfixed] <- apply(as.matrix(Z.MI[,index.alwaysfixed]),2,quantile, q.alwaysfixed)
        }

        newz <- rbind(point1, point2)

        preds <- preds.fun(newz)
        est.vec[k] <- drop(cc %*% preds$postmean)
        var.vec[k] <- drop(cc %*% preds$postvar %*% cc)
      }

      if(K==1){MIest <- c(est=est.vec, sd=sqrt(var.vec))}else{MIest <- Rubin.MI(mean.vec = est.vec, variance.vec = var.vec)}
      toreturn[i,"est"] <- MIest["est"]
      toreturn[i,"sd"]  <- MIest["sd"]

      end.time <- proc.time()["elapsed"]
      print(paste(i,"out of", length(qs), "complete: ", round((end.time - start.time)/60, digits=2), "min run time" ))
    }
  } else stop("method must be one of c('approx', 'exact')")
  toreturn
}



