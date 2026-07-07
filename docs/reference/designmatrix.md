# Generate design matrix

This function generates the design matrix for an item

## Usage

``` r
designmatrix(Kj = NULL, model = "GDINA", Qj = NULL, no.bugs = 0)
```

## Arguments

- Kj:

  Required except for the MS-DINA model; The number of attributes
  required for item j

- model:

  the model associated with the design matrix; It can be
  "GDINA","DINA","DINO", "ACDM","LLM", "RRUM", "MSDINA", "BUGDINO", and
  "SISM". The default is "GDINA". Note that models "LLM" and "RRUM" have
  the same design matrix as the "ACDM".

- Qj:

  the Q-matrix for item j; This is required for "MSDINA", and "SISM"
  models; The number of rows is equal to the number of strategies for
  "MSDINA", and the number of columns is equal to the number of
  attributes.

- no.bugs:

  the number of bugs (or misconceptions). Note that bugs must be given
  in the last no.bugs columns.

## Value

a design matrix (Mj). See de la Torre (2011) for details.

## References

de la Torre, J. (2011). The generalized DINA model framework.
*Psychometrika, 76*, 179-199.

## Examples

``` r
if (FALSE) { # \dontrun{
designmatrix(Kj = 2, model = "GDINA")
designmatrix(Kj = 3, model = "DINA")
msQj <- matrix(c(1,0,0,1,
                 1,1,0,0),nrow=2,byrow=TRUE)
designmatrix(model = "MSDINA",Qj = msQj)
} # }
```
