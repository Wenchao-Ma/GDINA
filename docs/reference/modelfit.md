# Model fit statistics

Calculate various absolute model-data fit statistics

## Usage

``` r
modelfit(GDINA.obj, CI = 0.9, ItemOnly = FALSE)
```

## Arguments

- GDINA.obj:

  An estimated model object of class `GDINA`

- CI:

  numeric value from 0 to 1 indicating the range of the confidence
  interval for RMSEA. Default returns the 90% interval.

- ItemOnly:

  should joint attribute distribution parameters be considered? Default
  = FALSE. See Ma (2019).

## Details

Various model-data fit statistics including M2 statistic for G-DINA
model with dichotmous responses (Liu, Tian, & Xin, 2016; Hansen, Cai,
Monroe, & Li, 2016) and for sequential G-DINA model with graded
responses (Ma, 2020). It also calculates SRMSR and RMSEA2.

## References

Hansen, M., Cai, L., Monroe, S., & Li, Z. (2016). Limited-information
goodness-of-fit testing of diagnostic classification item response
models. *British Journal of Mathematical and Statistical Psychology.
69,* 225–252.

Liu, Y., Tian, W., & Xin, T. (2016). An Application of M2 Statistic to
Evaluate the Fit of Cognitive Diagnostic Models. *Journal of Educational
and Behavioral Statistics, 41*, 3-26.

Ma, W. (2020). Evaluating the fit of sequential G-DINA model using
limited-information measures. *Applied Psychological Measurement, 44*,
167-181.

Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive
Diagnosis Modeling. *Journal of Statistical Software, 93(14)*, 1-26.

Maydeu-Olivares, A. (2013). Goodness-of-Fit Assessment of Item Response
Theory Models. *Measurement, 11*, 71-101.

## Author

Wenchao Ma, The University of Minnesota, <wma@umn.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
mod1 <- GDINA(dat = dat, Q = Q, model = "DINA")
modelfit(mod1)
} # }
```
