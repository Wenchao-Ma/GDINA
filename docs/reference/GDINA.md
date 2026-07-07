# CDM calibration under the G-DINA model framework

`GDINA` calibrates the generalized deterministic inputs, noisy and gate
(G-DINA; de la Torre, 2011) model for dichotomous responses, and its
extension, the sequential G-DINA model (Ma, & de la Torre, 2016a; Ma,
2017) for ordinal and nominal responses. By setting appropriate
constraints, the deterministic inputs, noisy and gate (DINA; de la
Torre, 2009; Junker & Sijtsma, 2001) model, the deterministic inputs,
noisy or gate (DINO; Templin & Henson, 2006) model, the reduced
reparametrized unified model (R-RUM; Hartz, 2002), the additive CDM
(A-CDM; de la Torre, 2011), the linear logistic model (LLM; Maris,
1999), and the multiple-strategy DINA model (MSDINA; de la Torre &
Douglas, 2008; Huo & de la Torre, 2014) can also be calibrated. Note
that the LLM is equivalent to the C-RUM (Hartz, 2002), a special case of
the GDM (von Davier, 2008), and that the R-RUM is also known as a
special case of the generalized NIDA model (de la Torre, 2011).

In addition, users are allowed to specify design matrix and link
function for each item, and distinct models may be used in a single test
for different items. The attributes can be either dichotomous or
polytomous (Chen & de la Torre, 2013). Joint attribute distribution may
be modelled using independent or saturated model, structured model,
higher-order model (de la Torre & Douglas, 2004), or loglinear model (Xu
& von Davier, 2008). Marginal maximum likelihood method with
Expectation-Maximization (MMLE/EM) alogrithm is used for item parameter
estimation.

To compare two or more `GDINA` objects, use method
[`anova`](https://rdrr.io/r/stats/anova.html).

To calculate structural parameters for item and joint attribute
distributions, use method [`coef`](https://rdrr.io/r/stats/coef.html).

To calculate lower-order incidental (person) parameters use method
[`personparm`](https://wenchao-ma.github.io/GDINA/reference/personparm.md).
To extract other components returned, use
[`extract`](https://wenchao-ma.github.io/GDINA/reference/extract.md). To
plot item/category response function, use
[`plot`](https://rdrr.io/r/graphics/plot.default.html). To check whether
monotonicity is violated, use
[`monocheck`](https://wenchao-ma.github.io/GDINA/reference/monocheck.md).
To conduct anaysis in graphical user interface, use
[`startGDINA`](https://wenchao-ma.github.io/GDINA/reference/startGDINA.md).

## Usage

``` r
GDINA(
  dat,
  Q,
  model = "GDINA",
  sequential = FALSE,
  att.dist = "saturated",
  mono.constraint = FALSE,
  group = NULL,
  linkfunc = NULL,
  design.matrix = NULL,
  no.bugs = 0,
  att.prior = NULL,
  att.str = NULL,
  verbose = 1,
  higher.order = list(),
  loglinear = 2,
  catprob.parm = NULL,
  control = list(),
  item.names = NULL,
  solver = NULL,
  nloptr.args = list(),
  auglag.args = list(),
  solnp.args = list(),
  ...
)

# S3 method for class 'GDINA'
anova(object, ...)

# S3 method for class 'GDINA'
coef(
  object,
  what = c("catprob", "delta", "gs", "itemprob", "LCprob", "rrum", "lambda"),
  withSE = FALSE,
  SE.type = 2,
  digits = 4,
  ...
)

# S3 method for class 'GDINA'
extract(object, what, SE.type = 2, ...)

# S3 method for class 'GDINA'
personparm(object, what = c("EAP", "MAP", "MLE", "mp", "HO"), digits = 4, ...)

# S3 method for class 'GDINA'
logLik(object, ...)

# S3 method for class 'GDINA'
deviance(object, ...)

# S3 method for class 'GDINA'
nobs(object, ...)

# S3 method for class 'GDINA'
vcov(object, ...)

# S3 method for class 'GDINA'
npar(object, ...)

# S3 method for class 'GDINA'
indlogLik(object, ...)

# S3 method for class 'GDINA'
indlogPost(object, ...)

# S3 method for class 'GDINA'
summary(object, ...)
```

## Arguments

- dat:

  A required \\N \times J\\ `matrix` or `data.frame` consisting of the
  responses of \\N\\ individuals to \\J\\ items. Missing values need to
  be coded as `NA`.

- Q:

  A required matrix; The number of rows occupied by a single-strategy
  dichotomous item is 1, by a polytomous item is the number of nonzero
  categories, and by a mutiple-strategy dichotomous item is the number
  of strategies. The number of column is equal to the number of
  attributes if all items are single-strategy dichotomous items, but the
  number of attributes + 2 if any items are polytomous or have multiple
  strategies. For a polytomous item, the first column represents the
  item number and the second column indicates the nonzero category
  number. For a multiple-strategy dichotomous item, the first column
  represents the item number and the second column indicates the
  strategy number. For binary attributes, 1 denotes the attributes are
  measured by the items and 0 means the attributes are not measured. For
  polytomous attributes, non-zero elements indicate which level of
  attributes are needed (see Chen, & de la Torre, 2013). See `Examples`.

- model:

  A vector for each item or nonzero category, or a scalar which will be
  used for all items or nonzero categories to specify the CDMs fitted.
  The possible options include
  `"GDINA"`,`"DINA"`,`"DINO"`,`"ACDM"`,`"LLM"`, `"RRUM"`, `"MSDINA"`,
  `"BUGDINO"`, `"SISM"`, and `"UDF"`. Note that model can also be
  `"logitGDINA"` and `"logGDINA"`, indicating the saturated G-DINA model
  in logit and log link functions. They are equivalent to the identity
  link saturated G-DINA model. The logit G-DINA model is identical to
  the log-linear CDM. When `"UDF"`, indicating user defined function, is
  specified for any item, arguments `design.matrix` and `linkfunc` need
  to be specified.

- sequential:

  logical; `TRUE` if the sequential model is fitted for polytomous
  responses.

- att.dist:

  How is the joint attribute distribution estimated? It can be (1)
  `saturated`, which is the default, indicating that the proportion
  parameter for each permissible latent class is estimated
  separately; (2) `higher.order`, indicating that a higher-order joint
  attribute distribution is assumed (higher-order model can be specified
  in `higher.order` argument); (3) `fixed`, indicating that the weights
  specified in `att.prior` argument are fixed in the estimation process.
  If `att.prior` is not specified, a uniform joint attribute
  distribution is employed initially; (4) `independent`, indicating that
  all attributes are assumed to be independent; and (5) `loglinear`,
  indicating a loglinear model is employed. If different groups have
  different joint attribute distributions, specify `att.dist` as a
  character vector with the same number of elements as the number of
  groups. However, if a higher-order model is used for any group, it
  must be used for all groups.

- mono.constraint:

  logical; `TRUE` indicates that \\P(\boldsymbol{\alpha}\_1)
  \<=P(\boldsymbol{\alpha}\_2)\\ if for all \\k\\, \\\alpha\_{1k} \<
  \alpha\_{2k}\\. Can be a vector for each item or nonzero category or a
  scalar which will be used for all items to specify whether
  monotonicity constraint should be added.

- group:

  a factor or a vector indicating the group each individual belongs to.
  Its length must be equal to the number of individuals.

- linkfunc:

  a vector of link functions for each item/category; It can be
  `"identity"`,`"log"` or `"logit"`. Only applicable when, for some
  items, `model="UDF"`.

- design.matrix:

  a list of design matrices; Its length must be equal to the number of
  items (or nonzero categories for sequential models). If CDM for item j
  is specified as "UDF" in argument `model`, the corresponding design
  matrix must be provided; otherwise, the design matrix can be `NULL`,
  which will be generated automatically.

- no.bugs:

  A numeric scalar (whole numbers only) indicating the number of bugs or
  misconceptions in the Q-matrix. The bugs must be included in the last
  `no.bugs` columns. It can be used along with the `BUGDINO` and `SISM`
  models (see, Kuo, Chen, Yang & Mok, 2016; Kuo, Chen, & de la Torre,
  2018). This argument will be ignored if the model is not specified in
  `model` argument. Note that the `BUGDINO` and `SISM` models are
  reparametrized - see `Details` below. By default, `no.bugs`=0,
  implying that there is no bugs/misconceptions.

- att.prior:

  A vector of length \\2^K\\ for single group model, or a matrix of
  dimension \\2^K\times\\ no. of groups to specify attribute prior
  distribution for \\2^K\\ latent classes for all groups under a
  multiple group model. Only applicable for dichotomous attributes. The
  sum of all elements does not have to be equal to 1; however, it will
  be normalized so that the sum is equal to 1 before calibration. The
  label for each latent class can be obtained by calling
  `attributepattern(K)`. See `examples` for more info.

- att.str:

  Specify attribute structures. `NULL`, by default, means there is no
  structure. Attribute structure needs be specified as a list - which
  will be internally handled by `att.structure` function. See examples.
  It can also be a matrix giving all permissible attribute profiles.

- verbose:

  How to print calibration information after each EM iteration? Can be
  0, 1 or 2, indicating to print no information, information for current
  iteration, or information for all iterations.

- higher.order:

  A list specifying the higher-order joint attribute distribution with
  the following components:

  - `model` - a number indicating the model for higher-order joint
    attribute distribution. Can be 1, 2 or 3, representing the intercept
    only approach, common slope approach and varied slope approach (see
    `Details`).

  - `nquad` - a scalar specifying the number of integral nodes. Default
    = 25.

  - `SlopeRange` - a vector of length two specifying the range of slope
    parameters. Default = \[0.1, 5\].

  - `InterceptRange` - a vector of length two specifying the range of
    intercept parameters. Default = \[-4, 4\].

  - `SlopePrior` - a vector of length two specifying the mean and
    variance of log(slope) parameters, which are assumed normally
    distributed. Default: mean = 0 and sd = 0.25.

  - `InterceptPrior` - a vector of length two specifying the mean and
    variance of intercept parameters, which are assumed normally
    distributed. Default: mean = 0 and sd = 1.

  - `Prior` - logical; indicating whether prior distributions should be
    imposed to slope and intercept parameters. Default is `FALSE`.

- loglinear:

  the order of loglinear smooth for attribute space. It can be either 1
  or 2 indicating the loglinear model with main effect only and with
  main effect and first-order interaction; It can also be a matrix,
  representing the design matrix for the loglinear model.

- catprob.parm:

  A list of initial success probability parameters for each nonzero
  category.

- control:

  A list of control parameters with elements:

  - `maxitr` A vector for each item or nonzero category, or a scalar
    which will be used for all items or nonzero categories to specify
    the maximum number of EM cycles allowed. Default = 2000.

  - `conv.crit` The convergence criterion. Default = 0.0001.

  - `conv.type` How is the convergence criterion evaluated? A vector
    with possible elements: `"ip"`, indicating the maximum absolute
    change in item success probabilities, `"mp"`, representing the
    maximum absolute change in mixing proportion parameters, `"delta"`,
    indicating the maximum absolute change in delta parameters, `neg2LL`
    indicating the absolute change in negative two times loglikeihood,
    or `neg2LL` indicating the relative absolute change in negative two
    times loglikeihood (i.e., the absolute change divided by -2LL of the
    previous iteration). Multiple criteria can be specified. If so, all
    criteria need to be met. Default = c("ip", "mp").

  - `nstarts` how many sets of starting values? Default = 3.

  - `lower.p` A vector for each item or nonzero category, or a scalar
    which will be used for all items or nonzero categories to specify
    the lower bound for success probabilities. Default = .0001.

  - `upper.p` A vector for each item or nonzero category, or a scalar
    which will be used for all items or nonzero categories to specify
    the upper bound for success probabilities. Default = .9999.

  - `lower.prior` The lower bound for mixing proportion parameters
    (latent class sizes). Default = .Machine\$double.eps.

  - `randomseed` Random seed for generating initial item parameters.
    Default = 123456.

  - `smallNcorrection` A numeric vector with two elements specifying the
    corrections applied when the expected number of individuals in some
    latent groups are too small. If the expected no. of examinees is
    less than the second element, the first element and two times the
    first element will be added to the numerator and denominator of the
    closed-form solution of probabilities of success. Only applicable
    for the G-DINA, DINA and DINO model estimation without monotonic
    constraints.

  - `MstepMessage` Integer; Larger number prints more information from
    Mstep optimizer. Default = 1.

- item.names:

  A vector giving the item names. By default, items are named as "Item
  1", "Item 2", etc.

- solver:

  A string indicating which solver should be used in M-step. By default,
  the solver is automatically chosen according to the models specified.
  Possible options include
  [slsqp](https://astamm.github.io/nloptr/reference/slsqp.html),
  [nloptr](https://astamm.github.io/nloptr/reference/nloptr.html),
  [solnp](https://rdrr.io/pkg/Rsolnp/man/solnp.html) and
  [auglag](https://rdrr.io/pkg/alabama/man/auglag.html).

- nloptr.args:

  a list of control parameters to be passed to `opts` argument of
  [nloptr](https://astamm.github.io/nloptr/reference/nloptr.html)
  function.

- auglag.args:

  a list of control parameters to be passed to the alabama::auglag()
  function. It can contain two elements: `control.outer` and
  `control.optim`. See
  [auglag](https://rdrr.io/pkg/alabama/man/auglag.html).

- solnp.args:

  a list of control parameters to be passed to `control` argument of
  [solnp](https://rdrr.io/pkg/Rsolnp/man/solnp.html) function.

- ...:

  additional arguments

- object:

  GDINA object for various S3 methods

- what:

  argument for various S3 methods; For calculating structural parameters
  using [`coef`](https://rdrr.io/r/stats/coef.html), `what` can be

  - `itemprob` - item success probabilities of each reduced attribute
    pattern.

  - `catprob` - category success probabilities of each reduced attribute
    pattern; the same as `itemprob` for dichtomous response data.

  - `LCprob` - item success probabilities of each attribute pattern.

  - `gs` - guessing and slip parameters of each item/category.

  - `delta` - delta parameters of each item/category, see G-DINA formula
    in details.

  - `rrum` - RRUM parameters when items are estimated using `RRUM`.

  - `lambda` - structural parameters for joint attribute distribution.

  For calculating incidental parameters using
  [`personparm`](https://wenchao-ma.github.io/GDINA/reference/personparm.md),
  `what` can be

  - `EAP` - EAP estimates of attribute pattern.

  - `MAP` - MAP estimates of attribute pattern.

  - `MLE` - MLE estimates of attribute pattern.

  - `mp` - marginal mastery probabilities.

  - `HO` - EAP estimates of higher-order ability if a higher-order is
    fitted.

- withSE:

  argument for method [`coef`](https://rdrr.io/r/stats/coef.html);
  estimate standard errors or not?

- SE.type:

  type of standard errors. For now, SEs are calculated based on
  outper-product of gradient. It can be `1` based on item-wise
  information, `2` based on incomplete information and `3` based on
  complete information.

- digits:

  How many decimal places in each number? The default is 4.

## Value

`GDINA` returns an object of class `GDINA`. Methods for `GDINA` objects
include
[`extract`](https://wenchao-ma.github.io/GDINA/reference/extract.md) for
extracting various components,
[`coef`](https://rdrr.io/r/stats/coef.html) for extracting structural
parameters,
[`personparm`](https://wenchao-ma.github.io/GDINA/reference/personparm.md)
for calculating incidental (person) parameters, `summary` for summary
information. `AIC`, `BIC`,`logLik`, `deviance` and `npar` can also be
used to calculate AIC, BIC, observed log-likelihood, deviance and number
of parameters.

## Methods (by generic)

- `anova(GDINA)`: Model comparison using likelihood ratio test

- `coef(GDINA)`: extract structural parameter estimates

- `extract(GDINA)`: extract various elements of GDINA estimates

- `personparm(GDINA)`: calculate person attribute patterns and
  higher-order ability

- `logLik(GDINA)`: calculate log-likelihood

- `deviance(GDINA)`: calculate deviance

- `nobs(GDINA)`: calculate number of observations

- `vcov(GDINA)`: calculate covariance-matrix for delta parameters

- `npar(GDINA)`: calculate the number of parameters

- `indlogLik(GDINA)`: extract log-likelihood for each individual

- `indlogPost(GDINA)`: extract log posterior for each individual

- `summary(GDINA)`: print summary information

## Note

anova function does NOT check whether models compared are nested or not.

## The G-DINA model

The generalized DINA model (G-DINA; de la Torre, 2011) is an extension
of the DINA model. Unlike the DINA model, which collaspes all latent
classes into two latent groups for each item, if item \\j\\ requires
\\K_j^\*\\ attributes, the G-DINA model collapses \\2^K\\ latent classes
into \\2^{K_j^\*}\\ latent groups with unique success probabilities on
item \\j\\, where \\K_j^\*=\sum\_{k=1}^{K}q\_{jk}\\.

Let \\\boldsymbol{\alpha}\_{lj}^\*\\ be the reduced attribute pattern
consisting of the columns of the attributes required by item \\j\\,
where \\l=1,\ldots,2^{K_j^\*}\\. For example, if only the first and the
last attributes are required,
\\\boldsymbol{\alpha}\_{lj}^\*=(\alpha\_{l1},\alpha\_{lK})\\. For
notational convenience, the first \\K_j^\*\\ attributes can be assumed
to be the required attributes for item \\j\\ as in de la Torre (2011).
The probability of success \\P(X\_{j}=1\|\boldsymbol{\alpha}\_{lj}^\*)\\
is denoted by \\P(\boldsymbol{\alpha}\_{lj}^\*)\\. To model this
probability of success, different link functions as in the generalized
linear models are used in the G-DINA model. The item response function
of the G-DINA model using the identity link can be written as \$\$
f\[P(\boldsymbol{\alpha}\_{lj}^\*)\]=\delta\_{j0}+\sum\_{k=1}^{K_j^\*}\delta\_{jk}\alpha\_{lk}+
\sum\_{k'=k+1}^{K_j^\*}\sum\_{k=1}^{K_j^\*-1}\delta\_{jkk'}\alpha\_{lk}\alpha\_{lk'}+\cdots+
\delta\_{j12{\cdots}K_j^\*}\prod\_{k=1}^{K_j^\*}\alpha\_{lk}, \$\$ or in
matrix form, \$\$
f\[\boldsymbol{P}\_j\]=\boldsymbol{M}\_j\boldsymbol{\delta}\_j, \$\$
where \\\delta\_{j0}\\ is the intercept for item \\j\\, \\\delta\_{jk}\\
is the main effect due to \\\alpha\_{lk}\\, \\\delta\_{jkk'}\\ is the
interaction effect due to \\\alpha\_{lk}\\ and \\\alpha\_{lk'}\\,
\\\delta\_{j12{\ldots}K_j^\*}\\ is the interaction effect due to
\\\alpha\_{l1}, \cdots,\alpha\_{lK_j^\*}\\. The log and logit links can
also be employed.

## Other CDMs as special cases

Several widely used CDMs can be obtained by setting appropriate
constraints to the G-DINA model. This section introduces the
parameterization of different CDMs within the G-DINA model framework
very breifly. Readers interested in this please refer to de la
Torre(2011) for details.

- `DINA model`:

  In DINA model, each item has two item parameters - guessing (\\g\\)
  and slip (\\s\\). In traditional parameterization of the DINA model, a
  latent variable \\\eta\\ for person \\i\\ and item \\j\\ is defined as
  \$\$\eta\_{ij}=\prod\_{k=1}^K\alpha\_{ik}^{q\_{jk}}\$\$ Briefly
  speaking, if individual \\i\\ master all attributes required by item
  \\j\\, \\\eta\_{ij}=1\\; otherwise, \\\eta\_{ij}=0\\. Item response
  function of the DINA model can be written by
  \$\$P(X\_{ij}=1\|\eta\_{ij})=(1-s_j)^{\eta\_{ij}}g_j^{1-\eta\_{ij}}\$\$
  To obtain the DINA model from the G-DINA model, all terms in identity
  link G-DINA model except \\\delta_0\\ and
  \\\delta\_{12{\ldots}K_j^\*}\\ need to be fixed to zero, that is,
  \$\$P(\boldsymbol{\alpha}\_{lj}^\*)=
  \delta\_{j0}+\delta\_{j12{\cdots}K_j^\*}\prod\_{k=1}^{K_j^\*}\alpha\_{lk}\$\$
  In this parameterization, \\\delta\_{j0}=g_j\\ and
  \\\delta\_{j0}+\delta\_{j12{\cdots}K_j^\*}=1-s_j\\.

- `DINO model`:

  The DINO model can be given by \$\$P(\boldsymbol{\alpha}\_{lj}^{\*})=
  \delta\_{j0}+\delta\_{j1}I(\boldsymbol{\alpha}\_{lj}^\*\neq
  \boldsymbol{0})\$\$

  where \\I(\cdot)\\ is an indicator variable. The DINO model is also a
  constrained identity link G-DINA model. As shown by de la Torre
  (2011), the appropriate constraint is
  \$\$\delta\_{jk}=-\delta\_{jk^{'}k^{''}}=\cdots=(-1)^{K_j^\*+1}\delta\_{j12{\cdots}K_j^\*},\$\$
  for \\k=1,\cdots,K_j^\*, k^{'}=1,\cdots,K_j^\*-1\\, and
  \\k^{''}\>k^{'},\cdots,K_j^\*\\.

- `Additive models with different link functions`:

  The A-CDM, LLM and R-RUM can be obtained by setting all interactions
  to be zero in identity, logit and log link G-DINA model, respectively.
  Specifically, the A-CDM can be formulated as
  \$\$P(\boldsymbol{\alpha}\_{lj}^\*)=
  \delta\_{j0}+\sum\_{k=1}^{K_j^\*}\delta\_{jk}\alpha\_{lk}\$\$ The item
  response function for LLM can be given by \$\$
  logit\[P(\boldsymbol{\alpha}\_{lj}^\*)\]=\delta\_{j0}+\sum\_{k=1}^{K_j^\*}\delta\_{jk}\alpha\_{lk},\$\$
  and lastly, the RRUM, can be written as
  \$\$log\[P(\boldsymbol{\alpha}\_{lj}^\*)\]=\delta\_{j0}+\sum\_{k=1}^{K_j^\*}\delta\_{jk}\alpha\_{lk}.\$\$
  It should be noted that the LLM is equivalent to the compensatory RUM,
  which is subsumed by the GDM, and that the RRUM is a special case of
  the generalized noisy inputs, deterministic “And" gate model (G-NIDA).

- `Simultaneously identifying skills and misconceptions (SISM)`:

  The SISM can be can be reformulated as
  \$\$P(\boldsymbol{\alpha}\_{lj}^\*)=
  \delta\_{j0}+\delta\_{j1}I\[{mastering all skills}\]+
  \delta\_{j2}I\[{having no bugs}\]+\delta\_{j12}I\[{mastering all
  skills and having no bugs}\].\$\$ As a result,the success probability
  of students who have mastered all the measured skills and possess none
  of the measured misconceptions (\\h_j\\ in Equation 4 of Kuo, et
  al, 2018) is \\\delta\_{j0}+\delta\_{j1}+\delta\_{j2}+\delta\_{j12}\\,
  the success probability of students who have mastered all the measured
  skills but possess some of the measured misconceptions (\\\omega_j\\)
  is \\\delta\_{j0}+\delta\_{j1}\\, the success probability of students
  who have not mastered all the measured skills and possess none of the
  measured misconceptions (\\g_j\\) is \\\delta\_{j0}+\delta\_{j2}\\ and
  success probability of students who have not mastered all the measured
  skills and possess at least one of the measured
  misconceptions(\\\epsilon_j\\) is \\\delta\_{j0}\\.

  By specifying `no.bugs` being equal to the number of attributes, the
  Bug-DINO is obtained, as in
  \$\$P(\boldsymbol{\alpha}\_{lj}^\*)=\delta\_{j0}+\delta\_{j1}I\[{having
  no bugs}\].\$\$

## Joint Attribute Distribution

The joint attribute distribution can be modeled using various methods.
This section mainly focuses on the so-called higher-order approach,
which was originally proposed by de la Torre and Douglas (2004) for the
DINA model. It has been extended in this package for all condensation
rules. Particularly, three approaches are available for the higher-order
attribute structure: intercept only approach, common slope approach and
varied slope approach. For the intercept only approach, the probability
of mastering attribute \\k\\ for individual \\i\\ is defined as
\$\$P(\alpha_k=1\|\theta_i,\lambda\_{0k})=\frac{exp(\theta_i+\lambda\_{0k})}{1+exp(\theta_i+\lambda\_{0k})}\$\$
For the common slope approach, the probability of mastering attribute
\\k\\ for individual \\i\\ is defined as
\$\$P(\alpha_k=1\|\theta_i,\lambda\_{0k},\lambda\_{1})=\frac{exp(\lambda\_{1}\theta_i+\lambda\_{0k})}{1+exp(\lambda\_{1}\theta_i+\lambda\_{0k})}\$\$
For the varied slope approach, the probability of mastering attribute
\\k\\ for individual \\i\\ is defined as
\$\$P(\alpha_k=1\|\theta_i,\lambda\_{0k},\lambda\_{1k})=\frac{exp(\lambda\_{1k}\theta_i+\lambda\_{0k})}{1+exp(\lambda\_{1k}\theta_i+\lambda\_{0k})}\$\$
where \\\theta_i\\ is the ability of examinee \\i\\. \\\lambda\_{0k}\\
and \\\lambda\_{1k}\\ are the intercept and slope parameters for
attribute \\k\\, respectively. The probability of joint attributes can
be written as
\$\$P(\boldsymbol{\alpha}\|\theta_i,\boldsymbol{\lambda})=\prod_k
P(\alpha_k\|\theta_i,\boldsymbol{\lambda})\$\$.

## Model Estimation

The MMLE/EM algorithm is implemented in this package. For G-DINA, DINA
and DINO models, closed-form solutions exist. See de la Torre (2009) and
de la Torre (2011) for details. For ACDM, LLM and RRUM, closed-form
solutions do not exist, and therefore some general optimization
techniques are adopted in M-step (Ma, Iaconangelo & de la Torre, 2016).
The selection of optimization techniques mainly depends on whether some
specific constraints need to be added.

The sequential G-DINA model is a special case of the diagnostic tree
model (DTM; Ma, 2019) and estimated using the mapping matrix accordingly
(See Tutz, 1997; Ma, 2019).

## The Number of Parameters

For dichotomous response models: Assume a test measures \\K\\ attributes
and item \\j\\ requires \\K_j^\*\\ attributes: The DINA and DINO model
has 2 item parameters for each item; if item \\j\\ is ACDM, LLM or RRUM,
it has \\K_j^\*+1\\ item parameters; if it is G-DINA model, it has
\\2^{K_j^\*}\\ item parameters. Apart from item parameters, the
parameters involved in the estimation of joint attribute distribution
need to be estimated as well. When using the saturated attribute
structure, there are \\2^K-1\\ parameters for joint attribute
distribution estimation; when using a higher-order attribute structure,
there are \\K\\, \\K+1\\, and \\2\times K\\ parameters for the intercept
only approach, common slope approach and varied slope approach,
respectively. For polytomous response data using the sequential G-DINA
model, the number of item parameters are counted at category level.

## References

Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation
of item parameters: Application of an EM algorithm. *Psychometrika, 46*,
443-459.

Bock, R. D., & Lieberman, M. (1970). Fitting a response model forn
dichotomously scored items. *Psychometrika, 35*, 179-197.

Bor-Chen Kuo, Chun-Hua Chen, Chih-Wei Yang, & Magdalena Mo Ching Mok.
(2016). Cognitive diagnostic models for tests with multiple-choice and
constructed-response items. *Educational Psychology, 36*, 1115-1133.

Carlin, B. P., & Louis, T. A. (2000). Bayes and empirical bayes methods
for data analysis. New York, NY: Chapman & Hall

de la Torre, J., & Douglas, J. A. (2008). Model evaluation and multiple
strategies in cognitive diagnosis: An analysis of fraction subtraction
data. *Psychometrika, 73*, 595-624.

de la Torre, J. (2009). DINA Model and Parameter Estimation: A Didactic.
*Journal of Educational and Behavioral Statistics, 34*, 115-130.

de la Torre, J. (2011). The generalized DINA model framework.
*Psychometrika, 76*, 179-199.

de la Torre, J., & Douglas, J. A. (2004). Higher-order latent trait
models for cognitive diagnosis. *Psychometrika, 69*, 333-353.

de la Torre, J., & Lee, Y. S. (2013). Evaluating the wald test for
item-level comparison of saturated and reduced models in cognitive
diagnosis. *Journal of Educational Measurement, 50*, 355-373.

Haertel, E. H. (1989). Using restricted latent class models to map the
skill structure of achievement items. *Journal of Educational
Measurement, 26*, 301-321.

Hartz, S. M. (2002). A bayesian framework for the unified model for
assessing cognitive abilities: Blending theory with practicality
(Unpublished doctoral dissertation). University of Illinois at
Urbana-Champaign.

Huo, Y., & de la Torre, J. (2014). Estimating a Cognitive Diagnostic
Model for Multiple Strategies via the EM Algorithm. *Applied
Psychological Measurement, 38*, 464-485.

Junker, B. W., & Sijtsma, K. (2001). Cognitive assessment models with
few assumptions, and connections with nonparametric item response
theory. *Applied Psychological Measurement, 25*, 258-272.

Kuo, B.-C., Chen.-H., & de la Torre,J. (2018). A cognitive diagnosis
model for identifying coexisting skills and misconceptions.*Applied
Psychological Measuremet*, 179–191.

Ma, W., & de la Torre, J. (2016). A sequential cognitive diagnosis model
for polytomous responses. *British Journal of Mathematical and
Statistical Psychology. 69,* 253-275.

Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive
Diagnosis Modeling. *Journal of Statistical Software, 93(14)*, 1-26.

Ma, W. (2019). A diagnostic tree model for polytomous responses with
multiple strategies. *British Journal of Mathematical and Statistical
Psychology, 72*, 61-82.

Ma, W., Iaconangelo, C., & de la Torre, J. (2016). Model similarity,
model selection and attribute classification. *Applied Psychological
Measurement, 40*, 200-217.

Ma, W. (2017). *A Sequential Cognitive Diagnosis Model for Graded
Response: Model Development, Q-Matrix Validation,and Model Comparison.
Unpublished doctoral dissertation.* New Brunswick, NJ: Rutgers
University.

Maris, E. (1999). Estimating multiple classification latent class
models. *Psychometrika, 64*, 187-212.

Tatsuoka, K. K. (1983). Rule space: An approach for dealing with
misconceptions based on item response theory. *Journal of Educational
Measurement, 20*, 345-354.

Templin, J. L., & Henson, R. A. (2006). Measurement of psychological
disorders using cognitive diagnosis models. *Psychological Methods, 11*,
287-305.

Tutz, G. (1997). Sequential models for ordered responses. In W.J. van
der Linden & R. K. Hambleton (Eds.), Handbook of modern item response
theory p. 139-152). New York, NY: Springer.

Xu, X., & von Davier, M. (2008). Fitting the structured general
diagnostic model to NAEP data. ETS research report, RR-08-27.

## See also

See
[`autoGDINA`](https://wenchao-ma.github.io/GDINA/reference/autoGDINA.md)
for Q-matrix validation, item-level model comparison and model
calibration in one run; See
[`modelfit`](https://wenchao-ma.github.io/GDINA/reference/modelfit.md)
and [`itemfit`](https://wenchao-ma.github.io/GDINA/reference/itemfit.md)
for model and item fit analysis,
[`Qval`](https://wenchao-ma.github.io/GDINA/reference/Qval.md) for
Q-matrix validation,
[`modelcomp`](https://wenchao-ma.github.io/GDINA/reference/modelcomp.md)
for item level model comparison and
[`simGDINA`](https://wenchao-ma.github.io/GDINA/reference/simGDINA.md)
for data simulation.
[`GMSCDM`](https://wenchao-ma.github.io/GDINA/reference/GMSCDM.md) for a
series of multiple strategy CDMs for dichotomous data, and
[`DTM`](https://wenchao-ma.github.io/GDINA/reference/DTM.md) for
diagnostic tree model for multiple strategies in polytomous response
data Also see `gdina` in CDM package for the G-DINA model estimation.

## Author

Wenchao Ma, The University of Minnesota, <wma@umn.edu> Jimmy de la
Torre, The University of Hong Kong

## Examples

``` r
if (FALSE) { # \dontrun{
####################################
#        Example 1.                #
#     GDINA, DINA, DINO            #
#    ACDM, LLM and RRUM            #
# estimation and comparison        #
#                                  #
####################################

dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ

#--------GDINA model --------#

mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
mod1
# summary information
summary(mod1)

AIC(mod1) #AIC
BIC(mod1) #BIC
logLik(mod1) #log-likelihood value
deviance(mod1) # deviance: -2 log-likelihood
npar(mod1) # number of parameters


head(indlogLik(mod1)) # individual log-likelihood
head(indlogPost(mod1)) # individual log-posterior

# structural parameters
# see ?coef
coef(mod1) # item probabilities of success for each latent group
coef(mod1, withSE = TRUE) # item probabilities of success & standard errors
coef(mod1, what = "delta") # delta parameters
coef(mod1, what = "delta",withSE=TRUE) # delta parameters
coef(mod1, what = "gs") # guessing and slip parameters
coef(mod1, what = "gs",withSE = TRUE) # guessing and slip parameters & standard errors

# person parameters
# see ?personparm
personparm(mod1) # EAP estimates of attribute profiles
personparm(mod1, what = "MAP") # MAP estimates of attribute profiles
personparm(mod1, what = "MLE") # MLE estimates of attribute profiles

#plot item response functions for item 10
plot(mod1,item = 10)
plot(mod1,item = 10,withSE = TRUE) # with error bars
#plot mastery probability for individuals 1, 20 and 50
plot(mod1,what = "mp", person =c(1,20,50))

# Use extract function to extract more components
# See ?extract

# ------- DINA model --------#
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
mod2 <- GDINA(dat = dat, Q = Q, model = "DINA")
mod2
coef(mod2, what = "gs") # guess and slip parameters
coef(mod2, what = "gs",withSE = TRUE) # guess and slip parameters and standard errors

# Model comparison at the test level via likelihood ratio test
anova(mod1,mod2)

# -------- DINO model -------#
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
mod3 <- GDINA(dat = dat, Q = Q, model = "DINO")
#slip and guessing
coef(mod3, what = "gs") # guess and slip parameters
coef(mod3, what = "gs",withSE = TRUE) # guess and slip parameters + standard errors

# Model comparison at test level via likelihood ratio test
anova(mod1,mod2,mod3)

# --------- ACDM model -------#
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
mod4 <- GDINA(dat = dat, Q = Q, model = "ACDM")
mod4
# --------- LLM model -------#
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
mod4b <- GDINA(dat = dat, Q = Q, model = "LLM")
mod4b
# --------- RRUM model -------#
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
mod4c <- GDINA(dat = dat, Q = Q, model = "RRUM")
mod4c

# --- Different CDMs for different items --- #

dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
models <- c(rep("GDINA",3),"LLM","DINA","DINO","ACDM","RRUM","LLM","RRUM")
mod5 <- GDINA(dat = dat, Q = Q, model = models)
anova(mod1,mod2,mod3,mod4,mod4b,mod4c,mod5)


####################################
#        Example 2.                #
#        Model estimations         #
# With monotonocity constraints    #
####################################

dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
# for item 10 only
mod11 <- GDINA(dat = dat, Q = Q, model = "GDINA",mono.constraint = c(rep(FALSE,9),TRUE))
mod11
mod11a <- GDINA(dat = dat, Q = Q, model = "DINA",mono.constraint = TRUE)
mod11a
mod11b <- GDINA(dat = dat, Q = Q, model = "ACDM",mono.constraint = TRUE)
mod11b
mod11c <- GDINA(dat = dat, Q = Q, model = "LLM",mono.constraint = TRUE)
mod11c
mod11d <- GDINA(dat = dat, Q = Q, model = "RRUM",mono.constraint = TRUE)
mod11d
coef(mod11d,"delta")
coef(mod11d,"rrum")

####################################
#           Example 3a.            #
#        Model estimations         #
# With Higher-order att structure  #
####################################

dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
# --- Higher order G-DINA model ---#
mod12 <- GDINA(dat = dat, Q = Q, model = "DINA",
               att.dist="higher.order",higher.order=list(nquad=31,model = "2PL"))
personparm(mod12,"HO") # higher-order ability
# structural parameters
# first column is slope and the second column is intercept
coef(mod12,"lambda")
# --- Higher order DINA model ---#
mod22 <- GDINA(dat = dat, Q = Q, model = "DINA", att.dist="higher.order",
               higher.order=list(model = "2PL",Prior=TRUE))

####################################
#           Example 3b.            #
#        Model estimations         #
#   With log-linear att structure  #
####################################

# --- DINA model with loglinear smoothed attribute space ---#
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
mod23 <- GDINA(dat = dat, Q = Q, model = "DINA",att.dist="loglinear",loglinear=1)
coef(mod23,"lambda") # intercept and three main effects

####################################
#           Example 3c.            #
#        Model estimations         #
#  With independent att structure  #
####################################

# --- GDINA model with independent attribute space ---#
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
mod33 <- GDINA(dat = dat, Q = Q, att.dist="independent")
coef(mod33,"lambda") # mastery probability for each attribute

####################################
#          Example 4.              #
#        Model estimations         #
#    With fixed att structure      #
####################################

# --- User-specified attribute priors ----#
# prior distribution is fixed during calibration
# Assume each of 000,100,010 and 001 has probability of 0.1
# and each of 110, 101,011 and 111 has probability of 0.15
# Note that the sum is equal to 1
#
prior <- c(0.1,0.1,0.1,0.1,0.15,0.15,0.15,0.15)
# fit GDINA model  with fixed prior dist.
dat <- sim10GDINA$simdat
Q <- sim10GDINA$simQ
modp1 <- GDINA(dat = dat, Q = Q, att.prior = prior, att.dist = "fixed")
extract(modp1, what = "att.prior")

####################################
#          Example 5a.             #
#             G-DINA               #
# with hierarchical att structure  #
####################################

# --- User-specified attribute structure ----#
Q <- sim30GDINA$simQ
K <- ncol(Q)
# divergent structure A1->A2->A3;A1->A4->A5
diverg <- list(c(1,2),
               c(2,3),
               c(1,4),
               c(4,5))
struc <- att.structure(diverg,K)
set.seed(123)
# data simulation
N <- 1000
true.lc <- sample(c(1:2^K),N,replace=TRUE,prob=struc$att.prob)
table(true.lc) #check the sample
true.att <- attributepattern(K)[true.lc,]
 gs <- matrix(rep(0.1,2*nrow(Q)),ncol=2)
 # data simulation
 simD <- simGDINA(N,Q,gs.parm = gs, model = "GDINA",attribute = true.att)
 dat <- extract(simD,"dat")

modp1 <- GDINA(dat = dat, Q = Q, att.str = diverg, att.dist = "saturated")
modp1
coef(modp1,"lambda")

####################################
#          Example 5b.             #
#   Reduced model (e.g.,ACDM)      #
# with hierarchical att structure  #
####################################

# --- User-specified attribute structure ----#
Q <- sim30GDINA$simQ
K <- ncol(Q)
# linear structure A1->A2->A3->A4->A5
linear <- list(c(1,2),
               c(2,3),
               c(3,4),
               c(4,5))
struc <- att.structure(linear,K)
set.seed(123)
# data simulation
N <- 1000
true.lc <- sample(c(1:2^K),N,replace=TRUE,prob=struc$att.prob)
table(true.lc) #check the sample
true.att <- attributepattern(K)[true.lc,]
 gs <- matrix(rep(0.1,2*nrow(Q)),ncol=2)
 # data simulation
 simD <- simGDINA(N,Q,gs.parm = gs, model = "ACDM",attribute = true.att)
 dat <- extract(simD,"dat")

modp1 <- GDINA(dat = dat, Q = Q, model = "ACDM",
               att.str = linear, att.dist = "saturated")
coef(modp1)
coef(modp1,"lambda")

####################################
#           Example 6.             #
# Specify initial values for item  #
# parameters                       #
####################################

 # check initials to see the format for initial item parameters
 initials <- sim10GDINA$simItempar

 dat <- sim10GDINA$simdat
 Q <- sim10GDINA$simQ
 mod.initial <- GDINA(dat,Q,catprob.parm = initials)

 # compare initial item parameters
 Map(rbind, initials,extract(mod.initial,"initial.catprob"))

####################################
#           Example 7a.            #
# Fix item and structure parameters#
# Estimate person attribute profile#
####################################

 # check initials to see the format for initial item parameters
 initials <- sim10GDINA$simItempar
 prior <- c(0.1,0.1,0.1,0.1,0.15,0.15,0.15,0.15)
 dat <- sim10GDINA$simdat
 Q <- sim10GDINA$simQ
 mod.ini <- GDINA(dat,Q,catprob.parm = initials,att.prior = prior,
                  att.dist = "fixed",control=list(maxitr = 0))
 personparm(mod.ini)
 # compare item parameters
 Map(rbind, initials,coef(mod.ini))

####################################
#           Example 7b.            #
#   Fix parameters for some items  #
# Estimate person attribute profile#
####################################

 # check initials to see the format for initial item parameters
 initials <- sim10GDINA$simItempar
 prior <- c(0.1,0.1,0.1,0.1,0.15,0.15,0.15,0.15)
 dat <- sim10GDINA$simdat
 Q <- sim10GDINA$simQ
 # fix parameters of the first 5 items; do not fix mixing proportion parameters
 mod.ini <- GDINA(dat,Q,catprob.parm = initials,
                  att.dist = "saturated",control=list(maxitr = c(rep(0,5),rep(2000,5))))
 personparm(mod.ini)
 # compare item parameters
 Map(rbind, initials,coef(mod.ini))

####################################
#           Example 8.             #
#        polytomous attribute      #
#          model estimation        #
#    see Chen, de la Torre 2013    #
####################################


# --- polytomous attribute G-DINA model --- #
dat <- sim30pGDINA$simdat
Q <- sim30pGDINA$simQ
#polytomous G-DINA model
pout <- GDINA(dat,Q)

# ----- polymous DINA model --------#
pout2 <- GDINA(dat,Q,model="DINA")
anova(pout,pout2)

####################################
#           Example 9.             #
#        Sequential G-DINA model   #
#    see Ma, & de la Torre 2016    #
####################################

# --- polytomous attribute G-DINA model --- #
dat <- sim20seqGDINA$simdat
Q <- sim20seqGDINA$simQ
Q
#    Item Cat A1 A2 A3 A4 A5
#       1   1  1  0  0  0  0
#       1   2  0  1  0  0  0
#       2   1  0  0  1  0  0
#       2   2  0  0  0  1  0
#       3   1  0  0  0  0  1
#       3   2  1  0  0  0  0
#       4   1  0  0  0  0  1
#       ...

#sequential G-DINA model
sGDINA <- GDINA(dat,Q,sequential = TRUE)
sDINA <- GDINA(dat,Q,sequential = TRUE,model = "DINA")
anova(sGDINA,sDINA)
coef(sDINA) # processing function
coef(sDINA,"itemprob") # success probabilities for each item
coef(sDINA,"LCprob") # success probabilities for each category for all latent classes

####################################
#           Example 10a.           #
#    Multiple-Group G-DINA model   #
####################################

Q <- sim10GDINA$simQ
K <- ncol(Q)
# parameter simulation
# Group 1 - female
N1 <- 3000
gs1 <- matrix(rep(0.1,2*nrow(Q)),ncol=2)
# Group 2 - male
N2 <- 3000
gs2 <- matrix(rep(0.2,2*nrow(Q)),ncol=2)

# data simulation for each group
sim1 <- simGDINA(N1,Q,gs.parm = gs1,model = "DINA",att.dist = "higher.order",
                 higher.order.parm = list(theta = rnorm(N1),
                 lambda = data.frame(a=rep(1.5,K),b=seq(-1,1,length.out=K))))
sim2 <- simGDINA(N2,Q,gs.parm = gs2,model = "DINO",att.dist = "higher.order",
                 higher.order.parm = list(theta = rnorm(N2),
                 lambda = data.frame(a=rep(1,K),b=seq(-2,2,length.out=K))))

# combine data - all items have the same item parameters
dat <- rbind(extract(sim1,"dat"),extract(sim2,"dat"))
gr <- rep(c(1,2),c(3000,3000))
# Fit G-DINA model
mg.est <- GDINA(dat = dat,Q = Q,group = gr)
summary(mg.est)
extract(mg.est,"posterior.prob")
coef(mg.est,"lambda")


####################################
#           Example 10b.           #
#    Multiple-Group G-DINA model   #
####################################

Q <- sim30GDINA$simQ
K <- ncol(Q)
# parameter simulation
N1 <- 3000
gs1 <- matrix(rep(0.1,2*nrow(Q)),ncol=2)
N2 <- 3000
gs2 <- matrix(rep(0.2,2*nrow(Q)),ncol=2)

# data simulation for each group
# two groups have different theta distributions
sim1 <- simGDINA(N1,Q,gs.parm = gs1,model = "DINA",att.dist = "higher.order",
                 higher.order.parm = list(theta = rnorm(N1),
                 lambda = data.frame(a=rep(1,K),b=seq(-2,2,length.out=K))))
sim2 <- simGDINA(N2,Q,gs.parm = gs2,model = "DINO",att.dist = "higher.order",
                 higher.order.parm = list(theta = rnorm(N2,1,1),
                 lambda = data.frame(a=rep(1,K),b=seq(-2,2,length.out=K))))

# combine data - different groups have distinct item parameters
# see ?bdiagMatrix
dat <- bdiagMatrix(list(extract(sim1,"dat"),extract(sim2,"dat")),fill=NA)
Q <- rbind(Q,Q)
gr <- rep(c(1,2),c(3000,3000))
mg.est <- GDINA(dat = dat,Q = Q,group = gr)
# Fit G-DINA model
mg.est <- GDINA(dat = dat,Q = Q,group = gr,att.dist="higher.order",
higher.order=list(model = "Rasch"))
summary(mg.est)
coef(mg.est,"lambda")
personparm(mg.est)
personparm(mg.est,"HO")
extract(mg.est,"posterior.prob")


####################################
#           Example 11.            #
#           Bug DINO model         #
####################################

set.seed(123)
Q <- sim10GDINA$simQ # 1 represents misconceptions/bugs
N <- 1000
J <- nrow(Q)

gs <- data.frame(guess=rep(0.1,J),slip=rep(0.1,J))

sim <- simGDINA(N,Q,gs.parm = gs,model = "BUGDINO")
dat <- extract(sim,"dat")
est <- GDINA(dat=dat,Q=Q,model = "BUGDINO")
coef(est)

####################################
#           Example 12.            #
#           SISM model             #
####################################

# The Q-matrix used in Kuo, et al (2018)
# The first four columns are for Attributes 1-4
# The last three columns are for Bugs 1-3
Q <- matrix(c(1,0,0,0,0,0,0,
0,1,0,0,0,0,0,
0,0,1,0,0,0,0,
0,0,0,1,0,0,0,
0,0,0,0,1,0,0,
0,0,0,0,0,1,0,
0,0,0,0,0,0,1,
1,0,0,0,1,0,0,
0,1,0,0,1,0,0,
0,0,1,0,0,0,1,
0,0,0,1,0,1,0,
1,1,0,0,1,0,0,
1,0,1,0,0,0,1,
1,0,0,1,0,0,1,
0,1,1,0,0,0,1,
0,1,0,1,0,1,1,
0,0,1,1,0,1,1,
1,0,1,0,1,1,0,
1,1,0,1,1,1,0,
0,1,1,1,1,1,0),ncol = 7,byrow = TRUE)

J <- nrow(Q)
N <- 1000
gs <- data.frame(guess=rep(0.1,J),slip=rep(0.1,J))

sim <- simGDINA(N,Q,gs.parm = gs,model = "SISM",no.bugs=3)
dat <- extract(sim,"dat")
est <- GDINA(dat=dat,Q=Q,model="SISM",no.bugs=3)
coef(est,"delta")

####################################
#           Example 13a.           #
#     user specified design matrix #
#        LCDM (logit G-DINA)       #
####################################

dat <- sim30GDINA$simdat
Q <- sim30GDINA$simQ

# LCDM
lcdm <- GDINA(dat = dat, Q = Q, model = "logitGDINA", control=list(conv.type="neg2LL"))

#Another way is to find design matrix for each item first => must be a list
D <- lapply(rowSums(Q),designmatrix,model="GDINA")
# for comparison, use change in -2LL as convergence criterion
# LCDM
lcdm2 <- GDINA(dat = dat, Q = Q, model = "UDF", design.matrix = D,
linkfunc = "logit", control=list(conv.type="neg2LL"),solver="slsqp")

# identity link GDINA
iGDINA <- GDINA(dat = dat, Q = Q, model = "GDINA",
control=list(conv.type="neg2LL"),solver="slsqp")

# compare all three models => identical
anova(lcdm,lcdm2,iGDINA)

####################################
#           Example 13b.           #
#     user specified design matrix #
#            RRUM                  #
####################################

dat <- sim30GDINA$simdat
Q <- sim30GDINA$simQ

# specify design matrix for each item => must be a list
# D can be defined by the user
D <- lapply(rowSums(Q),designmatrix,model="ACDM")
# for comparison, use change in -2LL as convergence criterion
# RRUM
logACDM <- GDINA(dat = dat, Q = Q, model = "UDF", design.matrix = D,
linkfunc = "log", control=list(conv.type="neg2LL"),solver="slsqp")

# identity link GDINA
RRUM <- GDINA(dat = dat, Q = Q, model = "RRUM",
              control=list(conv.type="neg2LL"),solver="slsqp")

# compare two models => identical
anova(logACDM,RRUM)

####################################
#           Example 14.            #
#     Multiple-strategy DINA model #
####################################

Q <- matrix(c(1,1,1,1,0,
1,2,0,1,1,
2,1,1,0,0,
3,1,0,1,0,
4,1,0,0,1,
5,1,1,0,0,
5,2,0,0,1),ncol = 5,byrow = TRUE)
d <- list(
  item1=c(0.2,0.7),
  item2=c(0.1,0.6),
  item3=c(0.2,0.6),
  item4=c(0.2,0.7),
  item5=c(0.1,0.8))

  set.seed(12345)
sim <- simGDINA(N=1000,Q = Q, delta.parm = d,
               model = c("MSDINA","MSDINA","DINA",
                         "DINA","DINA","MSDINA","MSDINA"))

# simulated data
dat <- extract(sim,what = "dat")
# estimation
# MSDINA need to be specified for each strategy
est <- GDINA(dat,Q,model = c("MSDINA","MSDINA","DINA",
                             "DINA","DINA","MSDINA","MSDINA"))
coef(est,"delta")
} # }
```
