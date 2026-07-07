# Tatsuoka's fraction subtraction data

Fraction Subtraction data (Tatsuoka, 1990, 2002) consists of responses
of 536 examinees to 20 items measuring 8 attributes.

## Usage

``` r
data(realdata_Tatsuoka1990)
```

## Format

A list of responses and Q-matrix with components:

- `dat`:

  responses of 536 examinees to 20 items

- `Q`:

  The \\20 \times 8\\ Q-matrix

## References

Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive
Diagnosis Modeling. *Journal of Statistical Software, 93(14)*, 1-26.
Tatsuoka, K. K. (1990). Toward an integration of item-response theory
and cognitive error diagnosis. In N. Frederiksen, R. Glaser, A. Lesgold,
& M. Shafto (Eds.), Diagnostic monitoring of skill and knowledge
acquisition (pp. 453-488). Hillsdale, NJ: Erlbaum. Tatsuoka, C. (2002).
Data analytic methods for latent partially ordered classification
models. *Journal of the Royal Statistical Society, Series C, Applied
Statistics, 51*, 337-350.

## Author

Wenchao Ma, The University of Minnesota, <wma@umn.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
mod1 <- GDINA(realdata_Tatsuoka1990$dat,realdata_Tatsuoka1990$Q,model="DINA")
mod1
summary(mod1)
# Higher order model
mod2 <- GDINA(realdata_Tatsuoka1990$dat,
realdata_Tatsuoka1990$Q,model="DINA",
att.dist="higher.order")
mod2
anova(mod1,mod2)
} # }
```
