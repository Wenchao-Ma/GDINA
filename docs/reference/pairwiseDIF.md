# Pairwise post hoc DIF analysis

Conduct pairwise Wald follow-up tests after omnibus DIF analysis with
three or more groups. The input can be a `dif` object returned by
[`dif()`](https://wenchao-ma.github.io/GDINA/reference/dif.md), in which
case DIF items are selected using adjusted omnibus p-values, or raw data
plus user-specified DIF items and optional anchor items.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) function can be
used to draw grouped bar chart for DIF items.

## Usage

``` r
pairwiseDIF(
  object = NULL,
  dat = NULL,
  Q = NULL,
  group = NULL,
  model = "GDINA",
  sequential = FALSE,
  dif.items = NULL,
  anchor.items = NULL,
  alpha.level = 0.05,
  p.adjust.methods = "holm",
  SE.type = NULL,
  ...
)
```

## Arguments

- object:

  an object returned by
  [`dif()`](https://wenchao-ma.github.io/GDINA/reference/dif.md).

- dat:

  item responses from two or more groups; missing data need to be coded
  as `NA`.

- Q:

  Q-matrix specifying the association between items and attributes.

- group:

  a factor or a vector indicating the group each individual belongs to.

- model:

  model for each item.

- sequential:

  Logical; whether a sequential model is fit to the data. Default is
  `FALSE`.

- dif.items:

  which items are subject to pairwise DIF follow-up.

- anchor.items:

  optional anchor items. If omitted, all non-DIF items are treated as
  shared across groups.

- alpha.level:

  adjusted omnibus p-value cutoff used to flag DIF items when `object`
  is supplied.

- p.adjust.methods:

  adjusted p-values for pairwise Wald tests within each item.

- SE.type:

  Type of standard error estimation methods for the Wald test.

- ...:

  arguments passed to
  [`GDINA()`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md) for
  the final post hoc refit.

## Value

A `pairwiseDIF` object containing the pairwise Wald test table and final
post hoc fit.

## See also

[`dif`](https://wenchao-ma.github.io/GDINA/reference/dif.md)

## Author

Wenchao Ma, The University of Minnesota, <wma@umn.edu>

## Examples

``` r
if (FALSE) { # \dontrun{

set.seed(123456)
N <- 500
Q <- sim30GDINA$simQ
gs <- matrix(.2,ncol = 2, nrow = nrow(Q))
# By default, individuals are simulated from uniform distribution
# and deltas are simulated randomly
sim1 <- simGDINA(N,Q,gs.parm = gs,model="DINA")
sim2 <- simGDINA(N,Q,gs.parm = gs,model=c(rep("DINA",nrow(Q)-1),"DINO"))
sim3 <- simGDINA(N,Q,gs.parm = gs,model="DINA")
dat <- rbind(extract(sim1,"dat"),extract(sim2,"dat"),extract(sim3,"dat"))
gr <- rep(c("G1","G2","G3"),each=N)

# DIF using Wald test - omnibus test
dif.wald <- dif(dat, Q, group=gr, method = "Wald")
dif.wald
# pairwise comparison between all pair of groups
dif.pair = pairwiseDIF(dif.wald)
dif.pair
# draw plot to show the DIF
plot(dif.pair, withSE = TRUE)
} # }
```
