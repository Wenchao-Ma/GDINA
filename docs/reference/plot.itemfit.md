# Item fit plots

Create plots of bivariate heatmap for item fit

## Usage

``` r
# S3 method for class 'itemfit'
plot(x, type = "all", adjusted = TRUE, ...)
```

## Arguments

- x:

  model object of class `itemfit`

- type:

  type of heatmap plot

- adjusted:

  logical; plot adjusted or unadjusted p-values?

- ...:

  additional arguments

## See also

[`GDINA`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md),
[`itemfit`](https://wenchao-ma.github.io/GDINA/reference/itemfit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ

fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
ift <- itemfit(fit)
# plot the adjusted p values for log odds or transformed correlation
plot(ift)
# plot unadjusted p values for log odds
plot(ift,adjusted = FALSE, type = "logOR")
} # }
```
