# Check EM convergence

Calculates maximum change across convergence criteria and checks whether
convergence has been achieved.

## Usage

``` r
check_em_convergence(dif_parm, conv_type, conv_crit, neg2LL_current)
```

## Arguments

- dif_parm:

  List with convergence differences: ip, prior, neg2LL, delt

- conv_type:

  Character vector of convergence types to check

- conv_crit:

  Numeric convergence criterion

- neg2LL_current:

  Current deviance value

## Value

List with maxchg and converged flag
