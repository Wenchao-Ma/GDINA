# calculate person (incidental) parameters

Function to calculate various person attribute parameters, including
`"EAP"`, `"MAP"`, and `"MLE"`, for EAP, MAP and MLE estimates of
attribute patterns (see Huebner & Wang, 2011), `"mp"` for marginal
mastery probabilities, and `"HO"` for higher-order ability estimates if
a higher-order model is fitted. See
[`GDINA`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md) for
examples.

## Usage

``` r
personparm(object, what = c("EAP", "MAP", "MLE", "mp", "HO"), digits = 4, ...)
```

## Arguments

- object:

  estimated GDINA object returned from
  [`GDINA`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md)

- what:

  what to extract; It can be `"EAP"`, `"MAP"`, and `"MLE"`, for EAP, MAP
  and MLE estimates of attribute patterns, and `"mp"` for marginal
  mastery probabilities, and `"HO"` for higher-order ability estimates
  if a higher-order model is fitted.

- digits:

  number of decimal places.

- ...:

  additional arguments

## References

Huebner, A., & Wang, C. (2011). A note on comparing examinee
classification methods for cognitive diagnosis models. *Educational and
Psychological Measurement, 71*, 407-419.

## Author

Wenchao Ma, The University of Alabama, <wenchao.ma@ua.edu>  
Jimmy de la Torre, The University of Hong Kong

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
# EAP
head(personparm(fit))
# MAP
head(personparm(fit, what = "MAP"))
} # }
```
