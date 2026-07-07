# The Generalized DINA Model Framework

For conducting CDM analysis within the G-DINA model framework

## Details

This package (Ma & de la Torre, 2020a) provides a framework for a series
of cognitively diagnostic analyses for dichotomous and polytomous
responses.

Various cognitive diagnosis models (CDMs) can be calibrated using the
[`GDINA`](https://wenchao-ma.github.io/GDINA/reference/GDINA.md)
function, including the G-DINA model (de la Torre, 2011), the
deterministic inputs, noisy and gate (DINA; de la Torre, 2009; Junker &
Sijtsma, 2001) model, the deterministic inputs, noisy or gate (DINO;
Templin & Henson, 2006) model, the reduced reparametrized unified model
(R-RUM; Hartz, 2002), the additive CDM (A-CDM; de la Torre, 2011), and
the linear logistic model (LLM; Maris, 1999), the multiple-strategy DINA
model (de la Torre, & Douglas, 2008) and models defined by users under
the G-DINA framework using different link functions and design matrices
(de la Torre, 2011). Note that the LLM is also called compensatory RUM
and the RRUM is equivalent to the generalized NIDA model.

For ordinal and nominal responses, the sequential G-DINA model (Ma, & de
la Torre, 2016) can be fitted and most of the aforementioned CDMs can be
used as the processing functions (Ma, & de la Torre, 2016) at the
category level. Different CDMs can be assigned to different items within
a single assessment. Item parameters are estimated using the MMLE/EM
algorithm. Details about the estimation algorithm can be found in Ma and
de la Torre (2020). The joint attribute distribution can be modeled
using an independent model, a higher-order IRT model (de la Torre, &
Douglas, 2004), a loglinear model (Xu & von Davier, 2008), a saturated
model or a hierarchical structures (e.g., linear, divergent).
Monotonicity constraints for item/category success probabilities can
also be specified.

In addition, to handle multiple strategies, generalized
multiple-strategy CDMs for dichotomous response (Ma & Guo, 2019) can be
fitted using
[`GMSCDM`](https://wenchao-ma.github.io/GDINA/reference/GMSCDM.md)
function and diagnostic tree model (Ma, 2019) can also be estimated
using [`DTM`](https://wenchao-ma.github.io/GDINA/reference/DTM.md)
function for polytomous responses. Note that these functions are
experimental, and are expected to be further extended in the future.
Other diagnostic approaches include the multiple-choice model (de la
Torre, 2009) and an iterative latent class analysis (ILCA; Jiang, 2019).

Various Q-matrix validation methods (de la Torre, & Chiu, 2016; de la
Torre & Ma, 2016; Ma & de la Torre, 2020b; Najera, Sorrel, & Abad, 2019;
see [`Qval`](https://wenchao-ma.github.io/GDINA/reference/Qval.md)),
model-data fit statistics (Chen, de la Torre, & Zhang, 2013; Hansen,
Cai, Monroe, & Li, 2016; Liu, Tian, & Xin, 2016; Ma, 2020; see
[`modelfit`](https://wenchao-ma.github.io/GDINA/reference/modelfit.md)
and
[`itemfit`](https://wenchao-ma.github.io/GDINA/reference/itemfit.md)),
model comparison at test and item level (de la Torre, 2011; de la Torre,
& Lee, 2013; Ma, Iaconangelo, & de la Torre, 2016; Ma & de la Torre,
2019; Sorrel, Abad, Olea, de la Torre, & Barrada, 2017; Sorrel, de la
Torre, Abad, & Olea, 2017; see
[`modelcomp`](https://wenchao-ma.github.io/GDINA/reference/modelcomp.md)),
and differential item functioning (Hou, de la Torre, & Nandakumar, 2014;
Ma, Terzi, Lee, & de la Torre, 2017; see
[`dif`](https://wenchao-ma.github.io/GDINA/reference/dif.md)) can also
be conducted.

To use the graphical user interface, check
[`startGDINA`](https://wenchao-ma.github.io/GDINA/reference/startGDINA.md).

## See also

CDM for estimating G-DINA model and a set of other CDMs; ACTCD and NPCD
for nonparametric CDMs; dina for DINA model in Bayesian framework

## Author

Wenchao Ma, The University of Minnesota, <wma@umn.edu> Jimmy de la
Torre, The University of Hong Kong
