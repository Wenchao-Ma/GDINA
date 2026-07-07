# Q-matrix validation

Q-matrix validation for the (sequential) G-DINA model based on PVAF (de
la Torre & Chiu, 2016; Najera, Sorrel, & Abad, 2019; Najera et al.,
2020), stepwise Wald test (Ma & de la Torre, 2020) or mesa plot (de la
Torre & Ma, 2016). All these methods are suitable for dichotomous and
ordinal response data. If too many modifications are suggested based on
the default PVAF method, you are suggested to try the stepwise Wald test
method, iterative procedures or predicted cutoffs. You should always
check the mesa plots for further examination.

## Usage

``` r
Qval(
  GDINA.obj,
  method = "PVAF",
  iter = "none",
  eps = 0.95,
  digits = 4,
  wald.args = list(),
  iter.args = list(empty.att = FALSE, max.iter = 150, verbose = FALSE)
)

# S3 method for class 'Qval'
extract(object, what = c("sug.Q", "varsigma", "PVAF", "eps", "Q"), ...)

# S3 method for class 'Qval'
summary(object, ...)
```

## Arguments

- GDINA.obj:

  an estimated model object of class `GDINA`

- method:

  which Q-matrix validation method is used? Can be either `"PVAF"` or
  `"wald"`.

- iter:

  implement the method iteratively? Can be `"none"` for non-iterative
  validation (by default), `"test"`, `"test.att"`, or `"item"` (Najera
  et al., 2020).

- eps:

  cutoff value for PVAF from 0 to 1. Default = 0.95. Note that it can
  also be -1, indicating the predicted cutoff based on Najera, Sorrel,
  and Abad (2019).

- digits:

  how many decimal places in each number? The default is 4.

- wald.args:

  a list of arguments for the stepwise Wald test method.

  SE.type

  :   type of covariance matrix for the Wald test

  alpha.level

  :   alpha level for the wald test

  GDI

  :   it can be 0, 1 or 2; 0 means GDI is not used to choose the
      attribute - when more than one attributes are significant, the one
      with the largest p-value will be selected; GDI=1 means the
      attribute with the largest GDI will be selected; GDI=2 means the
      q-vector with the largest GDI will be selected.

  verbose

  :   print detailed information or not?

  stepwise

  :   `TRUE` for stepwise approach and `FALSE` for forward approach

- iter.args:

  a list of arguments for the iterative implementation.

  empty.att

  :   can a Q-matrix with an empty attribute (i.e., measured by no
      items) be provided? Default is FALSE

  max.iter

  :   maximum number of iterations. Default is 150

  verbose

  :   print information after each iteration? Default is FALSE

- object:

  `Qval` objects for S3 methods

- what:

  argument for S3 method `extract` indicating what to extract; It can be
  `"sug.Q"` for suggested Q-matrix, `"Q"` for original Q-matrix,
  `"varsigma"` for varsigma index, and `"PVAF"` for PVAF.

- ...:

  additional arguments

## Value

An object of class `Qval`. Elements that can be extracted using
`extract` method include:

- sug.Q:

  suggested Q-matrix

- Q:

  original Q-matrix

- varsigma:

  varsigma index

- PVAF:

  PVAF

## Methods (by generic)

- `extract(Qval)`: extract various elements from `Qval` objects

- `summary(Qval)`: print summary information

## References

de la Torre, J. & Chiu, C-Y. (2016). A General Method of Empirical
Q-matrix Validation. *Psychometrika, 81*, 253-273.

de la Torre, J., & Ma, W. (2016, August). Cognitive diagnosis modeling:
A general framework approach and its implementation in R. A Short Course
at the Fourth Conference on Statistical Methods in Psychometrics,
Columbia University, New York.

Ma, W., & de la Torre, J. (2020). An empirical Q-matrix validation
method for the sequential G-DINA model. *British Journal of Mathematical
and Statistical Psychology, 73*, 142-163.

Najera, P., Sorrel, M. A., & Abad, F.J. (2019). Reconsidering cutoff
points in the general method of empirical Q-matrix validation.
*Educational and Psychological Measurement, 79*, 727-753.

Najera, P., Sorrel, M. A., de la Torre, J., & Abad, F. J. (2020).
Improving robustness in Q-matrix validation using an iterative and
dynamic procedure. *Applied Psychological Measurement*.

## See also

[`GDINA`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md)

## Author

Wenchao Ma, The University of Minnesota, <wma@umn.edu> Miguel A. Sorrel,
Universidad Autónoma de Madrid Jimmy de la Torre, The University of Hong
Kong

## Examples

``` r
if (FALSE) { # \dontrun{
################################
#
# Binary response
#
################################
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
Q[10,] <- c(0,1,0)

# Fit the G-DINA model
mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")

# Q-validation using de la Torre and Chiu's method
pvaf <- Qval(mod1,method = "PVAF",eps = 0.95)
pvaf
extract(pvaf,what = "PVAF")
#See also:
extract(pvaf,what = "varsigma")
extract(pvaf,what = "sug.Q")

# Draw mesa plots using the function plot

plot(pvaf,item=10)

#The stepwise Wald test
stepwise <- Qval(mod1,method = "wald")
stepwise
extract(stepwise,what = "PVAF")
#See also:
extract(stepwise,what = "varsigma")
extract(stepwise,what = "sug.Q")

#Set eps = -1 to determine the cutoff empirically
pvaf2 <- Qval(mod1,method = "PVAF",eps = -1)
pvaf2

#Iterative procedure (test-attribute level)
pvaf3 <- Qval(mod1, method = "PVAF", eps = -1,
              iter = "test.att", iter.args = list(verbose = 1))
pvaf3

################################
#
# Ordinal response
#
################################
seq.est <- GDINA(sim20seqGDINA$simdat,sim20seqGDINA$simQ, sequential = TRUE)
stepwise <- Qval(seq.est, method = "wald")
} # }
```
