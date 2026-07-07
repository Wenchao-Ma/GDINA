# Classification Rate Evaluation

This function evaluates the classification rates for two sets of
attribute profiles

## Usage

``` r
ClassRate(att1, att2)
```

## Arguments

- att1:

  a matrix or data frame of attribute profiles

- att2:

  a matrix or data frame of attribute profiles

## Value

a list with the following components:

- PCA:

  the proportion of correctly classified attributes (i.e., attribute
  level classification rate)

- PCV:

  a vector giving the proportions of correctly classified attribute
  vectors (i.e., vector level classification rate). The fist element is
  the proportion of at least one attribute in the vector are correctly
  identified; the second element is the proportion of at least two
  attributes in the vector are correctly identified; and so forth. The
  last element is the proportion of all elements in the vector are
  correctly identified.

## References

Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive
Diagnosis Modeling. *Journal of Statistical Software, 93(14)*, 1-26.

## Author

Wenchao Ma, The University of Minnesota, <wma@umn.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
N <- 2000
# model does not matter if item parameter is probability of success
Q <- sim30GDINA$simQ
J <- nrow(Q)
gs <- matrix(0.1,J,2)

set.seed(12345)
sim <- simGDINA(N,Q,gs.parm = gs)
GDINA.est <- GDINA(sim$dat,Q)

CR <- ClassRate(sim$attribute,personparm(GDINA.est))
CR
} # }
```
