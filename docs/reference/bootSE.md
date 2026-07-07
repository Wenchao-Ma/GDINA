# Calculating standard errors and variance-covariance matrix using bootstrap methods

This function conducts nonparametric and parametric bootstrap to
calculate standard errors of model parameters. Parametric bootstrap is
only applicable to single group models.

## Usage

``` r
bootSE(
  GDINA.obj,
  bootsample = 50,
  type = "nonparametric",
  cores = 1,
  seed = 12345
)
```

## Arguments

- GDINA.obj:

  an object of class GDINA

- bootsample:

  the number of bootstrap samples

- type:

  type of bootstrap method. Can be `parametric` or `nonparametric`

- cores:

  number of cores to be used for calculation; Default = 1; if `cores>1`,
  bootstrap will run in parallel with `foreach` package.

- seed:

  random seed for resampling

## Value

itemparm.se standard errors for item probability of success in list
format

delta.se standard errors for delta parameters in list format

lambda.se standard errors for structural parameters of joint attribute
distribution

boot.est resample estimates

## References

Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive
Diagnosis Modeling. *Journal of Statistical Software, 93(14)*, 1-26.

## Author

Wenchao Ma, The University of Minnesota, <wma@umn.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
# For illustration, only 5 resamples are run
# results are definitely not reliable

dat <- sim30GDINA$simdat
Q <- sim30GDINA$simQ
fit <- GDINA(dat = dat, Q = Q, model = "GDINA",att.dist = "higher.order")
boot.fit <- bootSE(fit,bootsample = 5,seed=123)
boot.fit$delta.se
boot.fit$lambda.se
} # }
```
