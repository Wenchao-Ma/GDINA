# Mesa plot for Q-matrix validation

The mesa plot was first proposed by de la Torre and Ma (2016) for
graphically illustrating the best q-vector(s) for each item. The
q-vector on the edge of the mesa is likely to be the best q-vector.

## Usage

``` r
# S3 method for class 'Qval'
plot(
  x,
  item,
  type = "best",
  no.qvector = 10,
  data.label = TRUE,
  eps = "auto",
  original.q.label = FALSE,
  auto.ylim = TRUE,
  ...
)
```

## Arguments

- x:

  model object of class `Qvalidation`

- item:

  a vector specifying which item(s) the plots are drawn for

- type:

  types of the plot. It can be `"best"` or `"all"`. If `"best"`, for all
  q-vectors requiring the same number of attributes, only the one with
  the largest PVAF is plotted, which means \\K_j\\ q-vectors are
  plotted; If `"all"`, all q-vectors will be plotted.

- no.qvector:

  the number of q vectors that need to be plotted when `type="all"`. The
  default is 10, which means the 10 q vectors with the largest PVAFs are
  plotted.

- data.label:

  logical; To show data label or not?

- eps:

  the cutoff for PVAF. If not `NULL`, it must be a value between 0
  and 1. A horizontal line will be drawn accordingly.

- original.q.label:

  logical; print the label showing the original q-vector or not?

- auto.ylim:

  logical; create y range automatically or not?

- ...:

  additional arguments passed to `plot` function

## References

de la Torre, J., & Ma, W. (2016, August). Cognitive diagnosis modeling:
A general framework approach and its implementation in R. A Short Course
at the Fourth Conference on Statistical Methods in Psychometrics,
Columbia University, New York.

## See also

[`Qval`](https://wenchao-ma.github.io/GDINA/reference/Qval.md),
[`autoGDINA`](https://wenchao-ma.github.io/GDINA/reference/autoGDINA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
Q[1,] <- c(0,1,0)
mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
out <- Qval(mod1,eps = 0.9)
item <- c(1,2,10)
plot(out,item=item,data.label=FALSE,type="all")
plot(out,item=10,type="best",eps=0.95)
plot(out,item=10,type="all",no.qvector=6)
} # }
```
