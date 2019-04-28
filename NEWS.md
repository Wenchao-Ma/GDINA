# GDINA 2.4.3
* Added   - `ICLA()` function for attribute profile estimation
* Changed - data are grouped before single group GDINA analysis 

# GDINA 2.4.1
* Fixed   - `print.GDINA`()` prints valid number of individuals by default now
* Fixed   - covariance matrix cannot be calculated when some individuals are deleted
* Fixed   - `Qval` not work when estimated number of individuals in some latent classes are 0
* Fixed   - `extract` only gives valid data

# GDINA 2.4.0
* Added   - Predicted cutoff based on the data for Q-matrix validation  
* Changed - add .5 for elements 0 when calculating log odds statistic for item fit
* Fixed   - model fit cannot be printed in GUI when M2 cannot be calculated

# GDINA 2.3.2
* Added   - Attribute Hierarchical structure for A-CDM, LLM and RRUM
* Added   - Estimating generalized multiple strategy CDMs using `GMSCDM`
* Changed - `att.str` argument of the `GDINA` function has been updated
* Changed - several examples in `GDINA` function
* Fixed   - The number of parameters is not correct when attribute strcture exists

# GDINA 2.2.0
* Changed - set anchor attributes for multiple group higher-order CDMs
* Changed - change arguments for `dif` function
* Changed - print estimated mixing proportions in R GUI
* Changed - GUI has been largely improved
* Changed - update `modelcomp` function to provide selected models directly
* Fixed   - calculation of the number of parameters when parameters of some items are fixed
* Fixed   - mesa plot with best q-vectors when some q-vectors have the same PVAF
* Fixed   - bug when calculating standard error using complete information based on multiple group models

# GDINA 2.1.15
* Changed - arguments for `plot` function
* Added   - additional example for `ecpe` data
* Changed - graphical interface

# GDINA 2.1.7
* Added   - experimental functions for simulating and estimating diagnostic tree model (to be optimized)
* Fixed   - bugs in `dif` function and `startGDINA` function
* Changed - `logLik` function to be consistent with default S3 methods
* Added   - Q-matrix validation using stepwise Wald test

# GDINA 2.0.7
* bug fixed

# GDINA 2.0
* This is a major update including a large number of new features. The `GDINA` function has been largely rewritten for both flexibility and speed. 
  Users are now allowed to fit `MS-DINA` model, `Bugs` models, and define models by providing design matrix and link functions. For joint attribute
  distribution, in addition to saturated and higher-order models, users are now allowed to fit independent model and loglinear model. Code for 
  joint attribute distribution modelling has been restructured as well. Other major updates include model fit evaluation using M2 statistics and other
  limited information measures, item-level model selection using likelihood ratio test and score test, classification accuracy evaluation indices, 
  bootstrap standard error estimation, etc.
* Note that due to the major updates, some results from this version may be slightly different from those using the previous versions. 

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



