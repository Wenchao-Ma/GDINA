# Combine R Objects by Columns

Combine a sequence of vector, matrix or data-frame arguments by columns.
Vector is treated as a column matrix.

## Usage

``` r
cjoint(..., fill = NA)
```

## Arguments

- ...:

  vectors or matrices

- fill:

  a scalar used when these objects have different number of rows.

## Value

a data frame

## See also

[`cbind`](https://rdrr.io/r/base/cbind.html)

## Examples

``` r
cjoint(2,c(1,2,3,4),matrix(1:6,2,3))
#>                
#> 1  2 1  1  3  5
#> 2 NA 2  2  4  6
#> 3 NA 3 NA NA NA
#> 4 NA 4 NA NA NA
cjoint(v1 = 2, v2 = c(3,2), v3 = matrix(1:6,3,2),
       v4 = data.frame(c(3,4,5,6,7),rep("x",5)),fill = 99)
#>                  
#> 1  2  3  1  4 3 x
#> 2 99  2  2  5 4 x
#> 3 99 99  3  6 5 x
#> 4 99 99 99 99 6 x
#> 5 99 99 99 99 7 x
```
