# Create a block diagonal matrix

Create a block diagonal matrix

## Usage

``` r
bdiagMatrix(mlist, fill = 0)
```

## Arguments

- mlist:

  a list of matrices

- fill:

  value to fill the non-diagnoal elements

## Value

a block diagonal matrix

## See also

`bdiag` in Matrix

## Examples

``` r
m1 <- bdiagMatrix(list(matrix(1:4, 2), diag(3)))
m2 <- bdiagMatrix(list(matrix(1:4, 2), diag(3)),fill = NA)
```
