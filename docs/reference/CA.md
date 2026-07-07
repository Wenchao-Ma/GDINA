# Calculate classification accuracy

This function calculates test-, pattern- and attribute-level
classification accuracy indices based on GDINA estimates from the
`GDINA` function using approaches in Iaconangelo (2017) and Wang, Song,
Chen, Meng, and Ding (2015). It is only applicable for dichotomous
attributes.

## Usage

``` r
CA(GDINA.obj, what = "MAP")
```

## Arguments

- GDINA.obj:

  estimated GDINA object returned from
  [`GDINA`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md)

- what:

  what attribute estimates are used? Default is `"MAP"`.

## Value

a list with elements

- tau:

  estimated overall pattern-level classification accuracy, see
  Iaconangelo (2017)

- tau_l:

  estimated pattern-level classification accuracy, see Iaconangelo
  (2017)

- tau_k:

  estimated attribute-level classification accuracy, see Wang, et al
  (2015)

- CCM:

  Conditional classification matrix, see Iaconangelo (2017)

- gamma_k:

  Wang et al.'s (2015) attribute-level classification consistency
  estimator

- gamma:

  estimated overall pattern-level classification consistency

## References

Iaconangelo, C.(2017). *Uses of Classification Error Probabilities in
the Three-Step Approach to Estimating Cognitive Diagnosis Models.*
(Unpublished doctoral dissertation). New Brunswick, NJ: Rutgers
University.

Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive
Diagnosis Modeling. *Journal of Statistical Software, 93(14)*, 1-26.

Wang, W., Song, L., Chen, P., Meng, Y., & Ding, S. (2015).
Attribute-Level and Pattern-Level Classification Consistency and
Accuracy Indices for Cognitive Diagnostic Assessment. *Journal of
Educational Measurement, 52* , 457-476.

## Author

Wenchao Ma, The University of Minnesota, <wma@umn.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- realdata_ECPE$dat
Q <- realdata_ECPE$Q
fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
fit
CA(fit)



} # }
```
