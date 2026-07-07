# Count the frequency of a row vector in a data frame

Count the frequency of a row vector in a data frame

## Usage

``` r
rowMatch(df, vec = NULL)
```

## Arguments

- df:

  a data frame or matrix

- vec:

  the vector for matching

## Value

count the number of vector vec in the data frame

row.no row numbers of the vector vec in the data frame

## Examples

``` r
df <- data.frame(V1=c(1L,2L),V2=LETTERS[1:3],V3=rep(1,12))
rowMatch(df,c(2,"B",1))
#> $count
#> [1] 2
#> 
#> $row.no
#> [1] 2 8
#> 


```
