# Extract log-likelihood for each individual

Extract individual log-likelihood.

## Usage

``` r
indlogLik(object, ...)
```

## Arguments

- object:

  GDINA object

- ...:

  additional arguments

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ

fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
iL <- indlogLik(fit)
iL[1:6,]
} # }
```
