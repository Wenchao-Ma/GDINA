# Transformation between latent classes and latent groups

This function gives the equivalent latent classes which have the same
category success probabilities for each category or item.

## Usage

``` r
LC2LG(Q, sequential = FALSE, att.str = NULL)
```

## Arguments

- Q:

  A required \\J \times K\\ binary Q-matrix. J represents test length
  and K represents the number of attributes of this test. Entry 1 at row
  j and column k represents the \\k^{th}\\ attribute is measured by item
  \\j\\, and 0 means item \\j\\ does not measure attribute \\k\\.

- sequential:

  logical; whether the Q-matrix is a Qc-matrix for sequential models?

- att.str:

  attribute structure. See `GDINA` for details.

## Value

An item or category by latent class matrix. In the G-DINA model, if item
j measures \\Kj\\ attributes, \\2^K\\ latent classes can be combined
into \\2^{Kj}\\ latent groups. This matrix gives which latent group each
of \\2^K\\ latent classes belongs to for each item.

## Author

Wenchao Ma, The University of Alabama, <wenchao.ma@ua.edu>  
Jimmy de la Torre, The University of Hong Kong

## Examples

``` r
attributepattern(3)
#>      A1 A2 A3
#> [1,]  0  0  0
#> [2,]  1  0  0
#> [3,]  0  1  0
#> [4,]  0  0  1
#> [5,]  1  1  0
#> [6,]  1  0  1
#> [7,]  0  1  1
#> [8,]  1  1  1

q <- matrix(scan(text = "0 1 0 1 0 1 1 1 0"),ncol = 3)
q
#>      [,1] [,2] [,3]
#> [1,]    0    1    1
#> [2,]    1    0    1
#> [3,]    0    1    0
LC2LG(Q = q)
#>      000 100 010 001 110 101 011 111
#> [1,]   1   1   2   3   2   3   4   4
#> [2,]   1   2   1   3   2   4   3   4
#> [3,]   1   1   2   1   2   1   2   2
```
