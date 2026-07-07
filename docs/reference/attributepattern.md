# Generate all possible attribute patterns

This function generates all possible attribute patterns. The Q-matrix
needs to be specified for polytomous attributes.

## Usage

``` r
attributepattern(K, Q)
```

## Arguments

- K:

  number of attributes

- Q:

  Q-matrix; required when Q-matrix is polytomous

## Value

A \\2^K\times K\\ matrix consisting of attribute profiles for \\2^K\\
latent classes

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

q <- matrix(scan(text = "0 1 2 1 0 1 1 2 0"),ncol = 3)
q
#>      [,1] [,2] [,3]
#> [1,]    0    1    1
#> [2,]    1    0    2
#> [3,]    2    1    0
attributepattern(Q=q)
#>       A1 A2 A3
#>  [1,]  0  0  0
#>  [2,]  0  0  1
#>  [3,]  0  0  2
#>  [4,]  0  1  0
#>  [5,]  0  1  1
#>  [6,]  0  1  2
#>  [7,]  1  0  0
#>  [8,]  1  0  1
#>  [9,]  1  0  2
#> [10,]  1  1  0
#> [11,]  1  1  1
#> [12,]  1  1  2
#> [13,]  2  0  0
#> [14,]  2  0  1
#> [15,]  2  0  2
#> [16,]  2  1  0
#> [17,]  2  1  1
#> [18,]  2  1  2

q <- matrix(scan(text = "0 1 1 1 0 1 1 1 0"),ncol = 3)
q
#>      [,1] [,2] [,3]
#> [1,]    0    1    1
#> [2,]    1    0    1
#> [3,]    1    1    0
attributepattern(K=ncol(q),Q=q)
#>      A1 A2 A3
#> [1,]  0  0  0
#> [2,]  1  0  0
#> [3,]  0  1  0
#> [4,]  0  0  1
#> [5,]  1  1  0
#> [6,]  1  0  1
#> [7,]  0  1  1
#> [8,]  1  1  1
```
