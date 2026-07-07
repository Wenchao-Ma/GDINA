# Calculate the number of parameters

Calculate the number of parameters for GDINA estimates. Returned the
total number of parameters, the number of item parameters and the number
parameters of joint attribute distribution.

## Usage

``` r
npar(object, ...)
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
npar(fit)
} # }
```
