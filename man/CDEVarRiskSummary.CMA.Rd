% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CDEVarRiskSummary.CMA.R
\name{CDEVarRiskSummary.CMA}
\alias{CDEVarRiskSummary.CMA}
\title{Single variable plot for CDE.}
\usage{
CDEVarRiskSummary.CMA(
  whichz = 1,
  BKMRfits,
  e.y = NULL,
  e.y.names = NULL,
  m.value = NULL,
  m.quant = c(0.1, 0.5, 0.75),
  m.name,
  qs.diff = c(0.25, 0.75),
  q.fixed = 0.5,
  alpha,
  sel,
  seed
)
}
\arguments{
\item{whichz}{vector indicating which variables (columns of Z) for which the summary should be computed}

\item{BKMRfits}{An object contatinint the results return by the kmbayes function}

\item{e.y}{effect modifier for the outcome variable}

\item{e.y.names}{column name of the effect modifier for the outcome variable}

\item{m.value}{values that the mediator is set to}

\item{m.quant}{values of the quantile that the mediator is set to}

\item{m.name}{column name of the mediator}

\item{qs.diff}{vector indicating the two quantiles q_1 and q_2 at which to compute h(z_{q2})-h(z_{q1})}

\item{q.fixed}{a second quantile at which to compare the estimated h function}

\item{alpha}{1-confidence interval}

\item{sel}{a vector selecting which iterations of the fit should be retained or inference}

\item{seed}{the random seed to use to evaluate the code}
}
\value{
a data frame containing the (posterior mean) estimate and posterior standard deviation of the single predictor CDE risk measures
}
\description{
A combination function of VarRiskSummary and riskSummary.approx for MI BKMR fits.
Compare estimated \code{h} function when a single variable (or a set of variables) is at the 75th versus 25th percentile, when all of the other variables are fixed at a particular percentile
}
\details{
For guided examples, go to https://zc2326.github.io/causalbkmr/articles/BKMRCMA_QuickStart.html
}
