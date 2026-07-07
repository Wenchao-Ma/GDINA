# Calculate Misclassification Matrices

This function computes profile-level and attribute-level
misclassification matrices from a fitted
[`GDINA`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md) model.
This function is only applicable to models with binary attributes.

## Usage

``` r
CM(object, classification = "MAP", matrixtype = "profile")
```

## Arguments

- object:

  An estimated GDINA object returned from
  [`GDINA`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md).

- classification:

  A character string specifying the classification rule. Supported
  values are `"MAP"`, `"MLE"`, and `"EAP"`. Alternatively, a matrix of
  user-supplied classifications can be provided, with one row per
  respondent and one column per attribute.

- matrixtype:

  A character string specifying which matrix to return. Supported values
  are `"profile"`, `"attribute"`, and `"both"`.

## Value

Depending on `matrixtype`, the function returns one of the following:

- `"profile"`:

  A class-by-class matrix of estimated \\P(\hat{\alpha} = c' \mid \alpha
  = c)\\.

- `"attribute"`:

  A list of \\2 \times 2\\ matrices, one for each attribute, with
  estimated \\P(\hat{\alpha}\_k = a' \mid \alpha_k = a)\\.

- `"both"`:

  A list with elements `profile_classification` and
  `att_classification`.

For both pattern-level and attribute-level matrices, rows correspond to
true classes and columns correspond to classified classes.

## Details

The profile-level matrix is computed from posterior latent class
probabilities. The attribute-level matrices are computed from marginal
attribute mastery probabilities and the chosen classifications.

## Author

Wenchao Ma, The University of Minnesota, <wma@umn.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- realdata_ECPE$dat
Q <- realdata_ECPE$Q
fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
CM(fit)
CM(fit, matrixtype = "both")
} # }
```
