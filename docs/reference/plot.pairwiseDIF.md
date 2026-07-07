# Grouped bar plots for pairwise DIF posthoc analysis

Create grouped bar charts of group-specific item/category response
probabilities for DIF items identified by
[`pairwiseDIF()`](https://wenchao-ma.github.io/GDINA/reference/pairwiseDIF.md).

## Usage

``` r
# S3 method for class 'pairwiseDIF'
plot(x, item = "all", withSE = FALSE, SE.type = 2, ...)
```

## Arguments

- x:

  model object of class
  [`pairwiseDIF`](https://wenchao-ma.github.io/GDINA/reference/pairwiseDIF.md)

- item:

  A scalar or vector specifying the DIF item(s) to plot.

- withSE:

  logical; Add error bar (estimate - SE, estimate + SE) to the plots?

- SE.type:

  How is SE estimated. By default, it's based on OPG using incomplete
  information.

- ...:

  additional arguments

## See also

[`pairwiseDIF`](https://wenchao-ma.github.io/GDINA/reference/pairwiseDIF.md),
[`dif`](https://wenchao-ma.github.io/GDINA/reference/dif.md)
