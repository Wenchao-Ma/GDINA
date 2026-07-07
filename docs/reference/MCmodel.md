# Multiple-choice models

This function estimates the multiple-choice DINA model (de la Torre,
2009).

## Usage

``` r
MCmodel(
  dat,
  Qc,
  model = "MCDINA",
  key = NULL,
  group = NULL,
  conv.crit = 0.001,
  maxitr = 2000,
  conv.type = "pr",
  SE = FALSE
)
```

## Arguments

- dat:

  A required \\N \times J\\ data matrix of N examinees to J items.
  Values must be 1, 2,... representing nominal categories. Missing
  values are currently not allowed.

- Qc:

  A required category and attribute association matrix. The first column
  gives the item number, which must be numeric (i.e., 1,2,...) and match
  the number of column in the data. The second column indicates the
  coded category of each item. The number in the second column must
  match with the number in the data, but if a category is not coded, it
  should not be included in the Q-matrix. Entry 1 indicates that the
  attribute is measured by the category, and 0 otherwise. Note that the
  MC-DINA model assumes that the category with the largest number of 1s
  is the key and that the coded distractors should allow to assign
  examinees uniquely.

- model:

  `MCDINA` only currently. Other MC models may be incorporated.

- key:

  a numeric vector giving the key of each item. See `Examples`. `NULL`
  by default indicating the coded category requiring the largest number
  of 1s is the key.

- group:

  Group membership vector if a multiple group model is considered. It
  must only contain whole numbers and start from 1. The length must be
  equal to the number of rows of the data.

- conv.crit:

  The convergence criterion for max absolute change in `conv.type` for
  two consecutive iterations.

- maxitr:

  The maximum iterations allowed.

- conv.type:

  convergence criteria; Can be `pr` or `LL`, indicating category
  response function, or -2 times log-likelihood,respectively.

- SE:

  logical; estimating standard error of item parameters? Default is
  `FALSE`.

## Value

an object of class `MCmodel` with the following components:

- prob.parm:

  A list of success probabilities for each reduced latent class on each
  item (IRF)

- prob.se:

  A list of standard errors of item parameters

- attribute:

  A list of estimated attribute profiles including EAP, MLE and MAP
  estimates.

- testfit:

  A list of test fit statistics including deviance, number of
  parameters, AIC and BIC

- R:

  expected \# of individuals in each latent group choosing each option

- lik:

  posterior probability

- itr:

  Total \# of iterations

## References

De La Torre, J. (2009). A cognitive diagnosis model for cognitively
based multiple-choice options. *Applied Psychological Measurement, 33*,
163–183.

Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive
Diagnosis Modeling. *Journal of Statistical Software, 93(14)*, 1-26.

## See also

[`GDINA`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md) for
G-DINA model

## Author

Wenchao Ma, The University of Minnesota, <wma@umn.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
 # check the format of the data
 # Entry 0 is not allowed
 head(sim10MCDINA$simdat)

 #---------------------------------
 # check the format of the Q-matrix
 #---------------------------------
 # Take item 1 as an example:
 # category 2 has a q-vector (1,0,0)
 # category 1 has a q-vector (0,1,0)
 # category 4 has a q-vector (1,1,0)
 # category 3 is not included in the Q-matrix because it is not coded
 # the order of the coded categories in the Q-matrix doesn't matter

 sim10MCDINA$simQ
  #     Item coded cat A1 A2 A3
  #        1         2  1  0  0
  #        1         1  0  1  0
  #        1         4  1  1  0
  #...
 est <- MCmodel(sim10MCDINA$simdat,sim10MCDINA$simQ)
 est
 est$testfit

 #--------------------------------------
 # Distractors involving more attributes
 #--------------------------------------
 # some distractors may involve attributes that are not invovled by the key option
 # this is not allowed by the "original" MC-DINA (de la Torre, 2009) but is allowed
 # in the current implementation

 # Users need to specify the key for each item to appropriate handle such an issue
 # Note item 1 below: category 1 is the key (as indicated in the key argument below)
 # The distractor (category 4) involves an attribute not included by the key option

 Qc <- matrix(c(1,  1,  0,  1,  1,
                1,  2,  0,  1,  0,
                1,  3,  1,  0,  0,
                1,  4,  1,  0,  1,
                2,  1,  1,  0,  0,
                2,  3,  1,  1,  0,
                2,  2,  1,  1,  1,
                3,  4,  1,  1,  1,
                3,  2,  1,  1,  0,
                3,  3,  0,  1,  1,
                4,  1,  0,  1,  1,
                4,  2,  0,  0,  1,
                5,  1,  1,  0,  0,
                6,  3,  0,  1,  0,
                7,  2,  0,  0,  1,
                8,  4,  1,  0,  0,
                9,  1,  0,  1,  0,
                10, 4,  0,  0,  1),ncol = 5,byrow = TRUE)

 est2 <- MCmodel(sim10MCDINA2$simdat,Qc, key = c(1,2,4,1,1,3,2,4,1,4))
 est2
 est2$prob.parm
 est2$testfit
 est2$attribute


 ###############################
 # a multiple group model
 ###############################
est <- MCmodel(sim10MCDINA$simdat, sim10MCDINA$simQ, group = rep(1:2, 1500))
# log posterior of different groups
est$lik$logprior

} # }
```
