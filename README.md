
<!-- README.md is generated from README.Rmd. Please edit that file -->
GDINA Package for Cognitive Diagnosis Modelling
===============================================

[![Project Status: Active ? The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html) [![Build Status](https://travis-ci.org/Wenchao-Ma/GDINA.svg?branch=master)](https://travis-ci.org/Wenchao-Ma/GDINA)

------------------------------------------------------------------------

[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.1.0-6666ff.svg)](https://cran.r-project.org/) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/GDINA)](https://cran.r-project.org/package=GDINA) [![packageversion](https://img.shields.io/badge/Package%20version-2.1.7-orange.svg?style=flat-square)](commits/master)

------------------------------------------------------------------------

[![Last-changedate](https://img.shields.io/badge/last%20change-2018--05--03-yellowgreen.svg)](/commits/master)

<!-- [![](http://www.r-pkg.org/badges/version/GDINA)](https://cran.r-project.org/package=GDINA) -->
<!-- [![Travis CI status](https://travis-ci.org/Wenchao-Ma/GDINA.svg?branch=master)](https://travis-ci.org/Wenchao-Ma/GDINA) -->
[![](https://cranlogs.r-pkg.org/badges/GDINA?color=brightgreen)](https://cran.r-project.org/package=GDINA)

[![](http://cranlogs.r-pkg.org/badges/grand-total/GDINA?color=green)](https://cran.r-project.org/package=GDINA)

Features of the package
-----------------------

-   Estimating G-DINA model and a variety of widely-used models subsumed by the G-DINA model, including the DINA model, DINO model, additive-CDM (A-CDM), linear logistic model (LLM), reduced reparametrized unified model (RRUM), multiple-strategy DINA model for dichotomous responses
-   Estimating models within the G-DINA model framework using user-specified design matrix and link functions
-   Estimating Bugs-DINA, DINO and G-DINA models for dichotomous responses
-   Estimating sequential G-DINA model for ordinal and nominal responses
-   Modelling independent, saturated, higher-order, loglinear smoothed, and structured joint attribute distribution
-   Accommodating multiple-group model analysis
-   Imposing monotonic constrained success probabilities
-   Accommodating binary and polytomous attributes
-   Validating Q-matrix under the general model framework
-   Evaluating absolute and relative item and model fit
-   Comparing models at the test and item levels
-   Detecting differential item functioning using Wald and likelihood ratio test
-   Simulating data based on all aforementioned CDMs
-   Providing graphical user interface for users less familiar with R

Installation
------------

You can install the development version of GDINA from github with:

``` r
# install.packages("devtools")
devtools::install_github("Wenchao-Ma/GDINA")
```

The stable version of GDINA should be installed from R CRAN at [here](https://CRAN.R-project.org/package=GDINA)

Upcoming Training Sessions
--------------------------

-   Frontiers in Educational Measurement pre-conference workshop on September 11, 2018 at the University of Oslo - see [here](http://www.uv.uio.no/cemo/english/conferences/fremo2018/pre-conference/)

Past Training Sessions
----------------------

-   NCME training session on April 12, 2018 - see [here](https://www.ncme.org/NCME)
-   NCME training session on April 26, 2017
-   Workshop at the Brazil Congress of Item Response Theory November 30-December 1, 2016
-   Pre-conference short course in the Fourth Conference on the Statistical Methods in Psychometrics from August 30 to September 1, 2016 at Columbia University [here](https://sites.google.com/site/2016psychometrics/)
