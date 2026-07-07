# Shared utility functions for estimation

This file contains helper functions shared between
SingleGroup_Estimation.R and MultipleGroup_Estimation.R to reduce code
duplication. Initialize control parameters for EM estimation

## Usage

``` r
init_control(control, ncat, model, is_mg = FALSE)
```

## Arguments

- control:

  User-provided control list

- ncat:

  Number of categories (nonzero categories in sequential models)

- model:

  Numeric model vector

- is_mg:

  Logical: is this multiple-group estimation?

## Value

Modified control list with defaults applied

## Details

Sets up default control values and validates user-provided overrides.
Handles both single-group and multiple-group estimation contexts.
