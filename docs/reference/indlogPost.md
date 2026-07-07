# Extract log posterior for each individual

Extract individual log posterior.

## Usage

``` r
indlogPost(object, ...)
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
iP <- indlogPost(fit)
iP[1:6,]
} # }
```
