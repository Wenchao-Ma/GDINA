# Create plots for GDINA estimates

Create various plots for GDINA estimates

## Usage

``` r
# S3 method for class 'GDINA'
plot(
  x,
  what = "IRF",
  item = "all",
  withSE = FALSE,
  SE.type = 2,
  person = 1,
  att.names = NULL,
  ...
)
```

## Arguments

- x:

  model object of class
  [`GDINA`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md)

- what:

  type of plot. Can be `"IRF"` for item/category response function plot,
  or `"mp"` for mastery probabilities for individuals.

- item:

  A scalar or vector specifying the item(s) for IRF plots.

- withSE:

  logical; Add error bar (estimate - SE, estimate + SE) to the IRF
  plots?

- SE.type:

  How is SE estimated. By default, it's based on OPG using incomplete
  information.

- person:

  A scalar or vector specifying the number of individuals for mastery
  plots.

- att.names:

  Optional; a vector for attribute names.

- ...:

  additional arguments

## See also

[`GDINA`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md),
[`autoGDINA`](https://wenchao-ma.github.io/GDINA/reference/autoGDINA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
#plot item response functions for item 10
plot(mod1, item = 10)
plot(mod1, what = "IRF", item = 10,withSE = TRUE)

# plot mastery probabilities for individuals 4 and 10
plot(mod1, what = "mp", person = c(4,10))
plot(mod1, what = "mp", person = c(4,10,15),
att.names = c("addition","subtraction","multiplication"))
} # }
```
