# Item-level model comparison using Wald, LR or LM tests

This function evaluates whether the saturated G-DINA model can be
replaced by reduced CDMs without significant loss in model data fit for
each item using the Wald test, likelihood ratio (LR) test or Lagrange
multiplier (LM) test. For Wald test, see de la Torre (2011), de la Torre
and Lee (2013), Ma, Iaconangelo and de la Torre (2016) and Ma & de la
Torre (2018) for details. For LR test and a two-step LR approximation
procedure, see Sorrel, de la Torre, Abad, and Olea (2017), Ma (2017) and
Ma & de la Torre (2019). For LM test, which is only applicable for DINA,
DINO and ACDM, see Sorrel, Abad, Olea, de la Torre, and Barrada (2017).
This function also calculates the dissimilarity between the reduced
models and the G-DINA model, which can be viewed as a measure of effect
size (Ma, Iaconangelo & de la Torre, 2016).

## Usage

``` r
modelcomp(
  GDINA.obj = NULL,
  method = "Wald",
  items = "all",
  p.adjust.methods = "holm",
  models = c("DINA", "DINO", "ACDM", "LLM", "RRUM"),
  decision.args = list(rule = "largestp", alpha.level = 0.05, adjusted = FALSE),
  DS = FALSE,
  Wald.args = list(SE.type = 2, varcov = NULL),
  LR.args = list(LR.approx = FALSE),
  LM.args = list(reducedMDINA = NULL, reducedMDINO = NULL, reducedMACDM = NULL, SE.type =
    2)
)

# S3 method for class 'modelcomp'
extract(
  object,
  what = c("stats", "pvalues", "adj.pvalues", "df", "DS", "selected.model"),
  digits = 4,
  ...
)

# S3 method for class 'modelcomp'
summary(object, ...)
```

## Arguments

- GDINA.obj:

  An estimated model object of class `GDINA`

- method:

  method for item level model comparison; can be `wald`, `LR` or `LM`.

- items:

  a vector of items to specify the items for model comparsion

- p.adjust.methods:

  adjusted p-values for multiple hypothesis tests. This is conducted
  using `p.adjust` function in stats, and therefore all adjustment
  methods supported by `p.adjust` can be used, including `"holm"`,
  `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BH"` and `"BY"`. See
  `p.adjust` for more details. `"holm"` is the default, indicating the
  Holm method.

- models:

  a vector specifying which reduced CDMs are possible reduced CDMs for
  each item. The default is "DINA","DINO","ACDM","LLM",and "RRUM".

- decision.args:

  a list of options for determining the most appropriate models
  including (1) `rule` can be either `"simpler"` or `"largestp"`. See
  details; (2) `alpha.level` for the nominal level of decision; and (3)
  `adjusted` can be either `TRUE` or `FALSE` indicating whether the
  decision is based on p value (`adjusted = FALSE`) or adjusted p
  values.

- DS:

  whether dissimilarity index should be calculated? `FALSE` is the
  default.

- Wald.args:

  a list of options for Wald test including (1) `SE.type` giving the
  type of covariance matrix for the Wald test; by default, it uses outer
  product of gradient based on incomplete information matrix; (2)
  `varcov` for user specified variance-covariance matrix. If supplied,
  it must be a list, giving the variance covariance matrix of success
  probability for each item or category. The default is `NULL`, in which
  case, the estimated variance-covariance matrix from the GDINA function
  is used.

- LR.args:

  a list of options for LR test including for now only `LR.approx`,
  which is either `TRUE` or `FALSE`, indicating whether a two-step LR
  approximation is implemented or not.

- LM.args:

  a list of options for LM test including `reducedMDINA`,
  `reducedMDINO`, and `reducedMACDM` for DINA, DINO and ACDM estimates
  from the `GDINA` function; `SE.type` specifies the type of covariance
  matrix.

- object:

  object of class `modelcomp` for various S3 methods

- what:

  argument for S3 method `extract` indicating what to extract; It can be
  `"wald"` for wald statistics, `"wald.p"` for associated p-values,
  `"df"` for degrees of freedom, and `"DS"` for dissimilarity between
  G-DINA and other CDMs.

- digits:

  How many decimal places in each number? The default is 4.

- ...:

  additional arguments

## Value

an object of class `modelcomp`. Elements that can be extracted using
`extract` method include

- stats:

  Wald or LR statistics

- pvalues:

  p-values associated with the test statistics

- adj.pvalues:

  adjusted p-values

- df:

  degrees of freedom

- DS:

  dissimilarity between G-DINA and other CDMs

## Details

After the test statistics for each reduced CDM were calculated for each
item, the reduced models with p values less than the pre-specified alpha
level were rejected. If all reduced models were rejected for an item,
the G-DINA model was used as the best model; if at least one reduced
model was retained, two diferent rules can be implemented for selecting
the best model specified in argument `decision.args`:

\(1\) when `rule="simpler"`,

If (a) the DINA or DINO model was one of the retained models, then the
DINA or DINO model with the larger p value was selected as the best
model; but if (b) both DINA and DINO were rejected, the reduced model
with the largest p value was selected as the best model for this item.
Note that when the p-values of several reduced models were greater than
0.05, the DINA and DINO models were preferred over the A-CDM, LLM, and
R-RUM because of their simplicity.

\(2\) When `rule="largestp"` (default),

The reduced model with the largest p-values is selected as the most
appropriate model.

## Methods (by generic)

- `extract(modelcomp)`: extract various elements from `modelcomp`
  objects

- `summary(modelcomp)`: print summary information

## References

de la Torre, J., & Lee, Y. S. (2013). Evaluating the wald test for
item-level comparison of saturated and reduced models in cognitive
diagnosis. *Journal of Educational Measurement, 50*, 355-373.

Ma, W., Iaconangelo, C., & de la Torre, J. (2016). Model similarity,
model selection and attribute classification. *Applied Psychological
Measurement, 40*, 200-217.

Ma, W. (2017). *A Sequential Cognitive Diagnosis Model for Graded
Response: Model Development, Q-Matrix Validation,and Model Comparison.
Unpublished doctoral dissertation.* New Brunswick, NJ: Rutgers
University.

Ma, W., & de la Torre, J. (2019). Category-Level Model Selection for the
Sequential G-DINA Model. *Journal of Educational and Behavioral
Statistics*. 44, 61-82.

Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive
Diagnosis Modeling. *Journal of Statistical Software, 93(14)*, 1-26.

Sorrel, M. A., Abad, F. J., Olea, J., de la Torre, J., & Barrada, J. R.
(2017). Inferential Item-Fit Evaluation in Cognitive Diagnosis Modeling.
*Applied Psychological Measurement, 41,* 614-631.

Sorrel, M. A., de la Torre, J., Abad, F. J., & Olea, J. (2017). Two-Step
Likelihood Ratio Test for Item-Level Model Comparison in Cognitive
Diagnosis Models. *Methodology, 13*, 39-47.

## See also

[`GDINA`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md),
[`autoGDINA`](https://wenchao-ma.github.io/GDINA/reference/autoGDINA.md)

## Author

Wenchao Ma, The University of Minnesota, <wma@umn.edu> Miguel A. Sorrel,
Universidad Autonoma de Madrid Jimmy de la Torre, The University of Hong
Kong

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
# --- GDINA model ---#
fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
fit

###################
#
# Wald test
#
###################

w <- modelcomp(fit)
w
# wald statistics
extract(w,"stats")
#p values
extract(w,"pvalues")
# selected models
extract(w,"selected.model")
##########################
#
# LR and Two-step LR test
#
##########################

lr <- modelcomp(fit,method = "LR")
lr
TwostepLR <- modelcomp(fit,items =c(6:10),method = "LR",LR.args = list(LR.approx = TRUE))
TwostepLR

##########################
#
# LM test
#
##########################

dina <- GDINA(dat = dat, Q = Q, model = "DINA")
dino <- GDINA(dat = dat, Q = Q, model = "DINO")
acdm <- GDINA(dat = dat, Q = Q, model = "ACDM")
lm <- modelcomp(method = "LM",LM.args=list(reducedMDINA = dina,
reducedMDINO = dino, reducedMACDM  = acdm))
lm


} # }
```
