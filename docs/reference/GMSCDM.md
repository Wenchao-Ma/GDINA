# Estimating multiple-strategy cognitive diagnosis models

An (experimental) function for calibrating the multiple-strategy CDMs
for dichotomous response data (Ma & Guo, 2019)

## Usage

``` r
GMSCDM(
  dat,
  msQ,
  model = "ACDM",
  s = 1,
  att.prior = NULL,
  delta = NULL,
  control = list()
)
```

## Arguments

- dat:

  A required binary item response matrix

- msQ:

  A multiple-strategy Q-matrix; the first column gives item numbers and
  the second column gives the strategy number. See examples.

- model:

  CDM used; can be `"DINA"`,`"DINO"`,`"ACDM"`,`"LLM"`, and `"RRUM"`,
  representing the GMS-DINA, GMS-DINO, GMS-ACDM, GMS-LLM and GMS-RRUM in
  Ma & Guo (2019), respectively. It can also be `"rDINA"` and `"rDINO"`,
  representing restricted GMS-DINA and GMS-DINO models where
  \\\delta\_{jm1}\\ are equal for all strategies. Note that only a
  single model can be used for the whole test.

- s:

  strategy selection parameter. It is equal to 1 by default.

- att.prior:

  mixing proportion parameters.

- delta:

  delta parameters in list format.

- control:

  a list of control arguments

## Value

an object of class `GMSCDM` with the following components:

- IRF:

  A matrix of success probabilities for each latent class on each item
  (IRF)

- delta:

  A list of delta parameters

- attribute:

  A list of estimated attribute profiles including EAP, MLE and MAP
  estimates.

- testfit:

  A list of test fit statistics including deviance, number of
  parameters, AIC and BIC

- sIRF:

  strategy-specific item response function

- pjmc:

  Probability of adopting each strategy on each item for each latent
  class

- sprv:

  Strategy pravelence

## References

Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive
Diagnosis Modeling. *Journal of Statistical Software, 93(14)*, 1-26.

Ma, W., & Guo, W. (2019). Cognitive Diagnosis Models for Multiple
Strategies. *British Journal of Mathematical and Statistical
Psychology.*

## See also

[`GDINA`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md) for
MS-DINA model and single strategy CDMs, and
[`DTM`](https://wenchao-ma.github.io/GDINA/reference/DTM.md) for
diagnostic tree model for multiple strategies in polytomous response
data

## Author

Wenchao Ma, The University of Minnesota, <wma@umn.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
##################
#
# data simulation
#
##################
set.seed(123)
msQ <- matrix(
c(1,1,0,1,
1,2,1,0,
2,1,1,0,
3,1,0,1,
4,1,1,1,
5,1,1,1),6,4,byrow = T)
# J x L - 00,10,01,11
LC.prob <- matrix(c(
0.2,0.7727,0.5889,0.8125,
0.1,0.9,0.1,0.9,
0.1,0.1,0.8,0.8,
0.2,0.5,0.4,0.7,
0.2,0.4,0.7,0.9),5,4,byrow=TRUE)
N <- 10000
att <- sample(1:4,N,replace=TRUE)
dat <- 1*(t(LC.prob[,att])>matrix(runif(N*5),N,5))


est <- GMSCDM(dat,msQ)
# item response function
est$IRF
# strategy specific IRF
est$sIRF


################################
#
# Example 14 from GDINA function
#
################################
Q <- matrix(c(1,1,1,1,0,
1,2,0,1,1,
2,1,1,0,0,
3,1,0,1,0,
4,1,0,0,1,
5,1,1,0,0,
5,2,0,0,1),ncol = 5,byrow = TRUE)
d <- list(
  item1=c(0.2,0.7),
  item2=c(0.1,0.6),
  item3=c(0.2,0.6),
  item4=c(0.2,0.7),
  item5=c(0.1,0.8))

  set.seed(123)
sim <- simGDINA(N=1000,Q = Q, delta.parm = d,
               model = c("MSDINA","MSDINA","DINA",
                         "DINA","DINA","MSDINA","MSDINA"))

# simulated data
dat <- extract(sim,what = "dat")
# estimation
# MSDINA need to be specified for each strategy
est <- GDINA(dat,Q,model = c("MSDINA","MSDINA","DINA",
                             "DINA","DINA","MSDINA","MSDINA"),
             control = list(conv.type = "neg2LL",conv.crit = .01))

# Approximate the MS-DINA model using GMS DINA model
est2 <- GMSCDM(dat, Q, model = "rDINA", s = 10,
               control = list(conv.type = "neg2LL",conv.crit = .01))
} # }

```
