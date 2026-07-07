# Score function

Calculate score function for each dichotomous item or each nonzero
category for polytomous items Only applicable to saturated model ofr
joint attribute distribution

## Usage

``` r
score(object, parm = "delta")
```

## Arguments

- object:

  an object of class GDINA

- parm:

  Either `delta` or `prob` indicating score function for delta
  parameters and success probabily parameters

## Value

a list where elements give the score functions for each item or category

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
score(fit)
} # }

```
