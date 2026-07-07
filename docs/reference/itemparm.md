# extract item parameters (deprecated)

This function has been deprecated; use `coef` instead.

## Usage

``` r
itemparm(
  object,
  what = c("catprob", "gs", "delta", "rrum", "itemprob", "LCprob"),
  withSE = FALSE,
  SE.type = 2,
  digits = 4,
  ...
)

# S3 method for class 'GDINA'
itemparm(
  object,
  what = c("catprob", "gs", "delta", "rrum", "itemprob", "LCprob"),
  withSE = FALSE,
  SE.type = 2,
  digits = 4,
  ...
)
```

## Arguments

- object:

  estimated GDINA object returned from
  [`GDINA`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md)

- what:

  what to show.

- withSE:

  show standard errors or not?

- SE.type:

  Type of standard errors.

- digits:

  how many decimal places for the output?

- ...:

  additional arguments

## References

Philipp, M., Strobl, C., de la Torre, J., & Zeileis, A.(2017). On the
estimation of standard errors in cognitive diagnosis models. *Journal of
Educational and Behavioral Statistics, 43*, 88-115.

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
# deprecated
itemparm(fit)
coef(fit)
} # }
```
