# extract elements from objects of various classes

A generic function to extract elements from objects of class `GDINA`,
`itemfit`, `modelcomp`, `Qval` or `simGDINA`. This page gives the
elements that can be extracted from the class `GDINA`. To see what can
be extracted from
[`itemfit`](https://wenchao-ma.github.io/GDINA/reference/itemfit.md),
[`modelcomp`](https://wenchao-ma.github.io/GDINA/reference/modelcomp.md),
and [`Qval`](https://wenchao-ma.github.io/GDINA/reference/Qval.md), go
to the corresponding function help page.

Objects which can be extracted from `GDINA` objects include:

- AIC:

  AIC

- att.prior:

  attribute prior weights for calculating marginalized likelihood in the
  last EM iteration

- attributepattern:

  all attribute patterns involved in the current calibration

- BIC:

  BIC

- CAIC:

  Consistent AIC

- catprob.cov:

  covariance matrix of item probability parameter estimates; Need to
  specify `SE.type`

- catprob.parm:

  item parameter estimates

- catprob.se:

  standard error of item probability parameter estimates; Need to
  specify `SE.type`

- convergence:

  `TRUE` if the calibration is converged.

- dat:

  raw data

- del.ind:

  deleted observation number

- delta.cov:

  covariance matrix of delta parameter estimates; Need to specify
  `SE.type`

- delta.parm:

  delta parameter estimates

- delta.se:

  standard error of delta parameter estimates; Need to specify `SE.type`

- designmatrix:

  A list of design matrices for each item/category

- deviance:

  deviance, or negative two times observed marginal log likelihood

- discrim:

  GDINA discrimination index

- expectedCorrect:

  expected \# of examinees in each latent group answering item correctly

- expectedTotal:

  expected \# of examinees in each latent group

- higher.order:

  higher-order model specifications

- LCprob.parm:

  success probabilities for all latent classes

- logLik:

  observed marginal log likelihood

- linkfunc:

  link functions for each item

- initial.catprob:

  initial item category probability parameters

- natt:

  number of attributes

- ncat:

  number of categories

- ngroup:

  number of groups

- nitem:

  number of items

- nitr:

  number of EM iterations

- nobs:

  number of observations, or sample size

- nLC:

  number of latent classes

- prevalence:

  prevalence of each attribute

- posterior.prob:

  posterior weights for each latent class

- reduced.LG:

  Reduced latent group for each item

- SABIC:

  Sample size Adusted BIC

- sequential:

  is a sequential model fitted?

## Usage

``` r
extract(object, what, ...)
```

## Arguments

- object:

  objects from class `GDINA`,`itemfit`, `modelcomp`, `Qval` or
  `simGDINA`

- what:

  what to extract

- ...:

  additional arguments

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
extract(fit,"discrim")
extract(fit,"designmatrix")
} # }
```
