# Item fit statistics from the power-divergence family

Calculate item fit statistics from the power-divergence family

## Usage

``` r
itemfitPD(
  GDINA.obj,
  lambda = 2/3,
  bootstrap = FALSE,
  R = 1000,
  Stone = FALSE,
  init.parm = FALSE,
  p.adjust.method = "holm",
  person.sim = "post",
  cores = 2,
  digits = 4,
  bound = 1e-10,
  seed = 123456
)
```

## Arguments

- GDINA.obj:

  Object containing a model fitted with the
  [`GDINA::GDINA`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md)
  function.

- lambda:

  Numeric; parameter for the power-divergence fit statistic.

- bootstrap:

  Logical; whether parametric bootstrap should be used.

- R:

  Integer; number of replicates in the bootstrap procedure if used.

- Stone:

  Logical; whether Stone indices should be computed (only available if
  `bootstrap = TRUE`).

- init.parm:

  Logical; whether the estimated item parameters are used in the
  estimation of the bootstrap replications.

- p.adjust.method:

  p-values can be adjusted for multiple comparisons at item level. This
  is conducted using `p.adjust` function in stats, and therefore all
  adjustment methods supported by `p.adjust` can be used, including
  `"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BH"` and `"BY"`.
  See `p.adjust` for more details. `"holm"` is the default.

- person.sim:

  Character; how to simulate attribute profiles in the bootstrap
  replications.

- cores:

  Integer; number of cores for parallelization during bootstrap.

- digits:

  Integer; number of decimal digits to report.

- bound:

  Numeric; minimum possible value for probabilities.

- seed:

  random seed.

## Value

an object of class `itemfitPD` consisting of several elements including:

- X2:

  Chi square statistics, adjusted and unadjusted p values for each item

- G2:

  G square statistics, adjusted and unadjusted p values for each item

- PD:

  PD statistics, adjusted and unadjusted p values for each item

- time:

  time used for the computation

## References

Najera, P., Ma, W., Sorrel, M. A. and Abad, F. J. (2025). Assessing
Item-Level Fit for the Sequential G-DINA Model.*Behaviormetrika*.

## Author

Pablo Najera Universidad Pontificia Comillas <pnajera@comillas.edu>

Wenchao Ma University of Minnesota <wma@umn.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ

mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
mod1
PDfit <- itemfitPD(mod1)
PDfit

dat <- sim21seqDINA$simdat
Q <- sim21seqDINA$simQ
sDINA <- GDINA(dat,Q,model="DINA",sequential = TRUE)
PDfit <- itemfitPD(sDINA)
PDfit
PDfit <- itemfitPD(sDINA, bootstrap = TRUE, Stone = TRUE, cores = 10)
PDfit
} # }
```
