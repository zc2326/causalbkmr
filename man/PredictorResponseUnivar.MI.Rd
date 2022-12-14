% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PredictorResponseUnivar.MI.R
\name{PredictorResponseUnivar.MI}
\alias{PredictorResponseUnivar.MI}
\title{Plot univariate predictor-response function on a new grid of point for MI BKMR}
\usage{
PredictorResponseUnivar.MI(
  BKMRfits,
  which.z = 1:ncol(BKMRfits[[1]]$Z),
  ngrid = 50,
  q.fixed = 0.5,
  sel = NULL,
  min.plot.dist = Inf,
  center = TRUE,
  method = "approx",
  ...
)
}
\arguments{
\item{BKMRfits}{A list of multiple BKMR fits and that each of these fits were ran for the same number of MCMC iterations.}

\item{which.z}{vector identifying which predictors (columns of \code{Z}) should be plotted}

\item{ngrid}{number of grid points to cover the range of each predictor (column in \code{Z})}

\item{q.fixed}{vector of quantiles at which to fix the remaining predictors in \code{Z}}

\item{sel}{logical expression indicating samples to keep; defaults to keeping the second half of all samples}

\item{min.plot.dist}{specifies a minimum distance that a new grid point needs to be from an observed data point in order to compute the prediction; points further than this will not be computed}

\item{center}{flag for whether to scale the exposure-response function to have mean zero}

\item{method}{method for obtaining posterior summaries at a vector of new points. Options are "approx" and "exact"; defaults to "approx", which is faster particularly for large datasets}

\item{...}{other arguments to pass on to the prediction function}
}
\value{
a long data frame with the predictor name, predictor value, posterior mean estimate, and posterior standard deviation
}
\description{
Plot univariate predictor-response function on a new grid of point for MI BKMR
}
\details{
For guided examples, go to https://zc2326.github.io/causalbkmr/articles/MI_BKMR.html
}
\examples{
\dontrun{
library(causalbkmr)
data(BKMRfits10)
univar.MI <- PredictorResponseUnivar.MI(BKMRfits10, ngrid = 50, q.fixed = 0.5, sel = sel.MI, method="approx")
ggplot(univar.MI, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) +
geom_smooth(stat = "identity") + ylab("h(z)") + facet_wrap(~ variable)+ggtitle("")
}


}
