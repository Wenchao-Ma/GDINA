# Unique values in a vector

Unique values in a vector

## Usage

``` r
unique_only(vec)
```

## Arguments

- vec:

  a vector

## Value

sorted unique values

## See also

[unique](https://rdrr.io/r/base/unique.html)

## Examples

``` r
vec <- c(4,2,3,5,4,4,4)
unique_only(vec)
#> [1] 2 3 5
# see the difference from unique
unique(vec)
#> [1] 4 2 3 5

vec <- letters[1:5]
unique_only(vec)
#> [1] "a" "b" "c" "d" "e"


```
