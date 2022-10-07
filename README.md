
<!-- README.md is generated from README.Rmd. Please edit that file -->

# causalbkmr: a suite of functions for causal Bayesian Kernel Machine Regression (causal-BKMR)

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/zc2326/causalbkmr.svg?branch=master)](https://travis-ci.com/zc2326/causalbkmr)
<!-- badges: end -->

## About the Package

The R package `causalbkmr` consists of three parts: `BKMR-CMA` for
Bayesian Kernel Machine Regression-Causal Mediation Analysis, `MI-BKMR`
for Multiple Imputation BKMR, and `g-BKMR` for Bayesian Kernel Machine
Regression for time-varying exposures and time-varying confounders.

We welcome your feedback and questions (email <zc2326@columbia.edu>)!

## Installation

You can install the latest version of `causalbkmr` via:

``` r
devtools::install_github("zc2326/causalbkmr")
```

Load `causalbkmr`:

``` r
library(causalbkmr) 
```

### BKMR-CMA

A command that implements the BKMR-CMA method, wich allows the
estimation of direct and indirect health effects of multiple
environmental exposures through a single mediator.

We use BKMR for the mediator and outcome regression models since BKMR
allows for all possible nonlinearities and interactions among the
elements included in the kernel with the specified outcome, without a
prior specification, and credible intervals a BKMR fit inherently
control for multiple testing due to the Bayesian nature of the model and
the prior specification.

We predict counterfactuals using the posterior predictive distributions
of the mediator and the outcome and present an algorithm for estimation
of mediation effects. We also conduct a simulation study to compare how
our approach performs relative to current mediation methods that assume
a restrictive linear relationship between the exposure, mediator, and
outcome.

Cite the paper: [Bayesian kernel machine regression-causal mediation
analysis](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.9255?casa_token=lf0hlEtqtjgAAAAA%3AyPLEllmuJGIeEN9ZMIr7kT33RjXQmjiLbUq4JjqfI5dxlSvkdeVjzkEiOqG9Rbh70Frxe3ONzA2aql8)

See the [Quick Start
guide](https://zc2326.github.io/causalbkmr/articles/BKMRCMA_QuickStart.html)
for `BKMR-CMA`. See the package
[website](https://zc2326.github.io/causalbkmr/articles/BKMRCMA_method_overview.html)
for an overview of statistical modeling approaches.

### BKMR-MI

A command that is used for valid estimation of environmental mixture
effects and evaluation of uncertainty in the presence of missing data,
which are imputed using multiple imputation techniques. The commands
combine information from multiple Bayesian kernel machine regression
(BKMR) models fit using the bkmr R package (Bobb et al. 2015, Valeri et
al. 2017, Bobb et al. 2018, Anglen Bauer et al. 2019). The commands also
produce effective visualizations of the estimated causal effects and
dose-response relationships. The package contains functions to be used
with BKMR MI fits to create a data frame for plotting with `ggplot`.

The data `BKMRfits10` is simulated data with 10 BKMRfits, each is a BKMR
fit using multiple imputed data with size *n* = 500.

See the [Quick Start
guide](https://zc2326.github.io/causalbkmr/articles/MI_BKMR.html) for
`BKMR-CMA`. See the package
[website](https://zc2326.github.io/causalbkmr/articles/BKMRMI_method_overview.html)
for an overview of statistical modeling approaches.

### g-BKMR

Exposure to environmental chemicals has been shown to rewire development
affecting later health status. Quantifying the joint effect of
environmental mixtures over time is crucial to determine intervention
timing. However, causal interpretation of longitudinal environmental
mixture studies encounters challenges. There is no statistical approach
that allows simultaneously for time-varying confounding, flexible
modeling, and variable selection when examining the effect of multiple,
correlated, and time-varying exposures. To address these gaps, we
develop a causal inference method, g-BKMR, which enables us to estimate
nonlinear, non-additive effects of time-varying exposures and
time-varying confounders, while also allowing for variable selection.
