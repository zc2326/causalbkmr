
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

A command that allows the estimation of direct and indirect health
effects of multiple environmental exposures through a single mediator.

Cite the paper: [Bayesian kernel machine regression-causal mediation
analysis](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.9255?casa_token=lf0hlEtqtjgAAAAA%3AyPLEllmuJGIeEN9ZMIr7kT33RjXQmjiLbUq4JjqfI5dxlSvkdeVjzkEiOqG9Rbh70Frxe3ONzA2aql8)

See the [Quick Start
guide](https://zc2326.github.io/causalbkmr/articles/BKMRCMA_QuickStart.html)
for `BKMR-CMA`. See the package
[website](https://zc2326.github.io/causalbkmr/articles/BKMRCMA_method_overview.html)
for an overview of statistical modeling approaches.

### MI-BKMR

A command that is used for valid estimation of environmental mixture
effects and evaluation of uncertainty in the presence of missing data,
which are imputed using multiple imputation techniques. The commands
also produce effective visualizations of the estimated causal effects
and dose-response relationships.

-   This code assumes you have K BKMR fits and that each of these fits
    were ran for the same number if MCMC iterations.

-   All effects are calculated as the average change in the outcome for
    a change in the exposure elements from a particular quantile to
    another quantile calculated across ALL imputed datasets. If there
    are no missing values in the Z matrix (in the mixture exposure),
    this is same contrast considered when only using the observed Z.  

-   All functions have the option to choose between an “approx” or
    “exact” method. The “exact” method combines the posterior samples
    from all MI fits and uses this posterior chain of length
    \#iterations times \#MI datasets for inference. The “approx” method
    uses the bkmr approx estimates and std errors from each MI fit and
    calculates an overall estimate and sd using Rubin’s 1987 method.
    (When using the “exact” method, the functions take a while to run,
    so make sure you save the data frames to be used for plotting).

-   This code can also be used with only 1 BKMR fit. The additional
    flexibility with this code over what is included in the bkmr package
    is that you have the option to fix specified variables to a given
    quantile when creating overall risk and single variable risk plots.

### g-BKMR

------------------------------------------------------------------------

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
