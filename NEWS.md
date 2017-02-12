# GDINA 1.2.1
* Fixed     - bugs in model estimation with user specified structures
* Fixed     - bugs in `summary.GDINA` function for multiple group estimation
* Changed   - include the pseudo q-vector 0 for the mesa plot
* Changed   - prior distribution is not re-calculated after likelihood calculation in `GDINA` function
* Changed   - output for `att.prior` in `extract` function
* Changed   - print number of group for `GDINA` function
* Changed   - documents in `simGDINA` and `GDINA`
* Changed   - examples in `GDINA`
* Added     - extract `ngroup` using `extract.GDINA` function

# GDINA 1.2.0
* Fixed   - infinite value issue during M-step optimization
* Fixed   - `itemfit` function for missing data
* Fixed   - adding monotonic constraints for the sequential models
* Fixed   - the maximum likelihood estimation of person attribute through `personparm`
* Changed - package imports, and suggests
* Changed - slsqp as the first optimizer used for monotonic G-DINA 
* Changed - Non zero prior probabilities

# GDINA 1.0.0

* Initial release



