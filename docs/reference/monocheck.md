# This function checks if monotonicity is violated

If mastering an additional attribute lead to a lower probabilities of
success, the monotonicity is violated.

## Usage

``` r
monocheck(object, strict = FALSE)
```

## Arguments

- object:

  object of class
  [`GDINA`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md)

- strict:

  whether a strict monotonicity is checked?

## Value

a logical vector for each item or category indicating whether the
monotonicity is violated (`TRUE`) or not (`FALSE`)

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ


mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
check <- monocheck(mod1)
check
mod2 <- GDINA(dat = dat, Q = Q, model = "GDINA",mono.constraint = check)
check2 <- monocheck(mod2)
check2
} # }
```
