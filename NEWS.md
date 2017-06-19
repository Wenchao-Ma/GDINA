# GDINA 1.4.4
* Changed   - rename function `rowCount` as `FreqTable`
* Changed   - `maxitr` in the `GDINA` function can be a vector, indicating different no. of iterations for different items
* Fixed     - when some observations are removed automatically, the standard errors may not be correctly calculated (thanks to Kazuhiro Yamaguchi)
* Fixed     - Format of SEs for DINO model was incorrect when printing category proabilities
* Fixed     - higher-order ability estimate was fixed; which was caused since v 1.4.1

# GDINA 1.4.2
* Fixed     - a bug in `modelcomp` in the version 1.4.1
* Changed   - adjusted p values are provided for DIF detection

# GDINA 1.4.1
* Changed   - update GUI interface `startGDINA` to report absolute fit statistics
* Changed   - estimation methods for higher-order models were adjusted
* Changed   - priors can be imposed to structural paprameters in `higher-order` models
* Changed   - missing values were identified when calculating likelihood
* Changed   - `GDINA` function arguments were modified
* Changed   - outputs for `npar`, `AIC`,`BIC`,`logLik` and `deviance` functions are modified
* Changed   - `anova` function is changed and multiple models can be comparied
* Changed   - `print` function for the GDINA function is changed
* Changed   - `summary` function is changed
* Changed   - objects returned from `extract` function are not rounded 
* Fixed     - bug when extracting covariance matrix and discrimination index (thanks to Kevin Carl Santos)
* Fixed     - bug in `dif` function when DIF items are specified
* Fixed     - bug for calculating standard errors for sequential models
* Fixed     - bug - recover the random seeds after running `GDINA` and `itemfit`
* Added     - standard errors using complete information matrix
* Added     - output for `initial.catprob` in `extract` function
* Added     - function `LC2LG` to find equivalent latent groups
* Added     - function `score` to find score functions
* Added     - Standard error estimation based on complete information is available 

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



