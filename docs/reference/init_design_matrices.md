# Initialize design matrices

Creates design matrices for all items based on model type and Q-matrix.

## Usage

``` r
init_design_matrices(
  ncat,
  model,
  rule,
  Kj,
  reduced.LG,
  Q,
  originalQ,
  no.bugs = NULL
)
```

## Arguments

- ncat:

  Number of categories

- model:

  Vector of model codes

- rule:

  Vector of condensation rules

- Kj:

  Vector of attribute counts per item

- reduced.LG:

  List of reduced latent group matrices

- Q:

  Q-matrix (attributes only, no item/category columns)

- originalQ:

  Original Q-matrix with item/category columns for MS-DINA

- no.bugs:

  Number of bugs for BUGDINO model

## Value

List of design matrices
