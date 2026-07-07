# Iterative latent-class analysis

This function implements an iterative latent class analysis (ILCA;
Jiang, 2019) approach to estimating attributes for cognitive diagnosis.

## Usage

``` r
ILCA(dat, Q, seed.num = 5)
```

## Arguments

- dat:

  A required binary item response matrix.

- Q:

  A required binary item and attribute association matrix.

- seed.num:

  seed number; Default = 5.

## Value

Estimated attribute profiles.

## References

Jiang, Z. (2019). Using the iterative latent-class analysis approach to
improve attribute accuracy in diagnostic classification models.
*Behavior research methods*, 1-10.

## Author

Zhehan Jiang, Peking University

## Examples

``` r
if (FALSE) { # \dontrun{
ILCA(sim10GDINA$simdat, sim10GDINA$simQ)
} # }
```
