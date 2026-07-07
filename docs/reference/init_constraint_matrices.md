# Initialize constraint matrices

Creates constraint matrices for monotonicity constraints across items.

## Usage

``` r
init_constraint_matrices(
  ncat,
  mono.constraint,
  Kj,
  reduced.LG,
  ConstrPairs = NULL
)
```

## Arguments

- ncat:

  Number of categories

- mono.constraint:

  Vector of logicals indicating which items have monotonicity
  constraints

- Kj:

  Vector of attribute counts per item

- reduced.LG:

  List of reduced latent group matrices

- ConstrPairs:

  Optional pre-specified constraint pairs

## Value

List with ConstrType vector and ConstrMatrix list
