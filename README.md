
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GDINA Package for Cognitively Diagnostic Analyses

[![Project Status: Active ? The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build
Status](https://travis-ci.org/Wenchao-Ma/GDINA.svg?branch=master)](https://travis-ci.org/Wenchao-Ma/GDINA)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/GDINA)](https://cran.r-project.org/package=GDINA)
[![](https://cranlogs.r-pkg.org/badges/GDINA?color=brightgreen)](https://cran.r-project.org/package=GDINA)
[![](http://cranlogs.r-pkg.org/badges/grand-total/GDINA?color=green)](https://cran.r-project.org/package=GDINA)

## How to cite the package

Ma, W. & de la Torre, J. (2019). GDINA: The generalized DINA model
framework. R package version 2.7.3. Retrived from
<https://CRAN.R-project.org/package=GDINA>

Visit the package website <https://wenchao-ma.github.io/GDINA> for
examples, tutorials and more information.

## Learning resources

  - Watch [Ma and de la Torre’s (2019) NCME digital
    module 5](https://ncme.elevate.commpartners.com/) on a gentle
    introduction to the G-DINA model framework and the use of graphical
    user interface for CDM analyses

  - Check [de la Torre and Akbay’s (2019)
    article](https://doi.org/10.14689/ejer.2019.80.9) on how to conduct
    various CDM analyses using the graphical user interface

## Features of the package

  - Estimating G-DINA model and a variety of widely-used models subsumed
    by the G-DINA model, including the DINA model, DINO model,
    additive-CDM (A-CDM), linear logistic model (LLM), reduced
    reparametrized unified model (RRUM), multiple-strategy DINA model
    for dichotomous responses
  - Estimating models within the G-DINA model framework using
    user-specified design matrix and link functions
  - Estimating Bugs-DINA, DINO and G-DINA models for dichotomous
    responses
  - Estimating sequential G-DINA model for ordinal and nominal responses
  - Estimating the generalized multiple-strategy cognitive diagnosis
    models (experimental)
  - Estimating the diagnostic tree model (experimental)
  - Modelling independent, saturated, higher-order, loglinear smoothed,
    and structured joint attribute distribution
  - Accommodating multiple-group model analysis
  - Imposing monotonic constrained success probabilities
  - Accommodating binary and polytomous attributes
  - Validating Q-matrix under the general model framework
  - Evaluating absolute and relative item and model fit
  - Comparing models at the test and item levels
  - Detecting differential item functioning using Wald and likelihood
    ratio test
  - Simulating data based on all aforementioned CDMs
  - Providing graphical user interface for users less familiar with R

## Installation

This package depends on R version 3.5.0 or newer

The stable version of GDINA should be installed from R CRAN at
[here](https://CRAN.R-project.org/package=GDINA)

To install this package from source:

1)  Windows users may need to install the
    [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/) and include
    the checkbox option of installing Rtools to their path for easier
    command line usage. Mac users will have to download the necessary
    tools from the
    [Xcode](https://itunes.apple.com/ca/app/xcode/id497799835?mt=12) and
    its related command line tools (found within Xcode’s Preference Pane
    under Downloads/Components); most Linux distributions should already
    have up to date compilers (or if not they can be updated easily).

2)  Install the `devtools` package (if necessary), and install the
    package from the Github source code.

<!-- end list -->

``` r
# install.packages("devtools")
devtools::install_github("Wenchao-Ma/GDINA")
```
