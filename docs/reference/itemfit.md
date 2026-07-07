# Item fit statistics

Calculate item fit statistics (Chen, de la Torre, & Zhang, 2013) and
draw heatmap plot for item pairs. For polytomous response data and
sequential models, the log odds is calculated by dichotomizing the data
by a cutoff (default:\>0), transformed correlation is based on the raw
data.

## Usage

``` r
itemfit(
  GDINA.obj,
  person.sim = "post",
  p.adjust.methods = "holm",
  cor.use = "pairwise.complete.obs",
  seq.cut = 0,
  digits = 4,
  N.resampling = NULL,
  randomseed = 123456
)

# S3 method for class 'itemfit'
extract(object, what, ...)

# S3 method for class 'itemfit'
summary(object, ...)
```

## Arguments

- GDINA.obj:

  An estimated model object of class `GDINA`

- person.sim:

  Simulate expected responses from the posterior or based on EAP, MAP
  and MLE estimates.

- p.adjust.methods:

  p-values for the proportion correct (mean response), transformed
  correlation, and log-odds ratio can be adjusted for multiple
  comparisons at test and item level. This is conducted using `p.adjust`
  function in stats, and therefore all adjustment methods supported by
  `p.adjust` can be used, including `"holm"`, `"hochberg"`, `"hommel"`,
  `"bonferroni"`, `"BH"` and `"BY"`. See `p.adjust` for more details.
  `"holm"` is the default.

- cor.use:

  how to deal with missing values when calculating correlations? This
  argument will be passed to `use` when calling
  [`stats::cor`](https://rdrr.io/r/stats/cor.html).

- seq.cut:

  the cutoff for dichotomizing the data for sequential models
  (default=0). A response that is great than `seq.cut` is converted to
  1; otherwise 0.

- digits:

  How many decimal places in each number? The default is 4.

- N.resampling:

  the sample size of resampling. By default, it is the maximum of 1e+5
  and ten times of current sample size.

- randomseed:

  random seed; This is used to make sure the results are replicable. The
  default random seed is 123456.

- object:

  objects of class `itemfit` for various S3 methods

- what:

  argument for S3 method `extract` indicating what to extract; It can be
  `"p"` for proportion correct statistics, `"r"` for transformed
  correlations, `logOR` for log odds ratios and `"maxitemfit"` for
  maximum statistics for each item.

- ...:

  additional arguments

## Value

an object of class `itemfit` consisting of several elements that can be
extracted using method `extract`. Components that can be extracted
include:

- p:

  the proportion correct (mean) statistics, adjusted and unadjusted p
  values for each item

- r:

  the transformed correlations, adjusted and unadjusted p values for
  each item pair

- logOR:

  the log odds ratios, adjusted and unadjusted p values for each item
  pair

- maxitemfit:

  the maximum proportion correct, transformed correlation, and log-odds
  ratio for each item with associated item-level adjusted p-values

## Methods (by generic)

- `extract(itemfit)`: extract various elements from `itemfit` objects

- `summary(itemfit)`: print summary information

## References

Chen, J., de la Torre, J., & Zhang, Z. (2013). Relative and Absolute Fit
Evaluation in Cognitive Diagnosis Modeling. *Journal of Educational
Measurement, 50*, 123-140.

Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive
Diagnosis Modeling. *Journal of Statistical Software, 93(14)*, 1-26.

## Author

Wenchao Ma, The University of Minnesota, <wma@umn.edu> Jimmy de la
Torre, The University of Hong Kong

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ

mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
mod1
itmfit <- itemfit(mod1)

# Print "test-level" item fit statistics
# p-values are adjusted for multiple comparisons
# for proportion correct, there are J comparisons
# for log odds ratio and transformed correlation,
# there are J*(J-1)/2 comparisons

itmfit

# The following gives maximum item fit statistics for
# each item with item level p-value adjustment
# For each item, there are J-1 comparisons for each of
# log odds ratio and transformed correlation
summary(itmfit)

# use extract to extract various components
extract(itmfit,"r")

mod2 <- GDINA(dat,Q,model="DINA")
itmfit2 <- itemfit(mod2)
#misfit heatmap
plot(itmfit2)
itmfit2
} # }
```
