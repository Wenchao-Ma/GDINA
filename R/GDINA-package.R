#' The Generalized DINA Model Framework
#'
#' For conducting CDM analysis within the G-DINA model framework
#'
#' This package (Ma & de la Torre, 2020a) provides a framework for a series of cognitively diagnostic analyses
#' for dichotomous and polytomous responses.
#'
#' Various cognitive
#' diagnosis models (CDMs) can be calibrated using the \code{\link{GDINA}}
#' function, including the G-DINA model (de la Torre, 2011), the deterministic inputs,
#' noisy and gate (DINA; de la Torre, 2009; Junker & Sijtsma, 2001) model,
#' the deterministic inputs, noisy or gate (DINO; Templin & Henson, 2006)
#' model, the reduced reparametrized unified model (R-RUM; Hartz, 2002),
#' the additive CDM (A-CDM; de la Torre, 2011), and the linear logistic
#' model (LLM; Maris, 1999), the multiple-strategy DINA model (de la Torre, & Douglas, 2008) and models defined
#' by users under the G-DINA framework using different link functions and design
#' matrices (de la Torre, 2011). Note that the LLM is also called
#' compensatory RUM and the RRUM is equivalent to the generalized NIDA model.
#'
#' For ordinal and nominal responses,
#' the sequential G-DINA model (Ma, & de la Torre, 2016) can be fitted and most of the
#' aforementioned CDMs can be used as the processing functions (Ma, & de la Torre, 2016) at the category level.
#' Different CDMs can be assigned to different items within a single assessment.
#' Item parameters are estimated using the MMLE/EM algorithm. Details about the estimation algorithm
#' can be found in Ma and de la Torre (2020).
#' The joint attribute distribution can be modeled using an independent model,
#' a higher-order IRT model (de la Torre, & Douglas, 2004), a loglinear model (Xu & von Davier, 2008),
#' a saturated model or a hierarchical structures (e.g., linear, divergent). Monotonicity constraints for item/category success
#' probabilities can also be specified.
#'
#' In addition, to handle multiple strategies, generalized multiple-strategy CDMs for dichotomous response (Ma & Guo, 2019) can be fitted using \code{\link{GMSCDM}} function and
#' diagnostic tree model (Ma, 2019) can also be estimated using \code{\link{DTM}} function for polytomous responses. Note that these functions are experimental, and are expected to be further extended
#' in the future. Other diagnostic approaches include the multiple-choice model (de la Torre, 2009) and an iterative latent class analysis (ILCA; Jiang, 2019).
#'
#' Various Q-matrix validation methods (de la Torre, & Chiu, 2016; de la Torre & Ma, 2016; Ma & de la Torre, 2020b; Najera, Sorrel, & Abad, 2019; see \code{\link{Qval}}),
#' model-data fit statistics (Chen, de la Torre, & Zhang, 2013; Hansen, Cai, Monroe, & Li, 2016; Liu, Tian, & Xin, 2016; Ma, 2020; see \code{\link{modelfit}} and \code{\link{itemfit}}),
#' model comparison at test and item level (de la Torre, 2011; de la Torre, & Lee, 2013;
#' Ma, Iaconangelo, & de la Torre, 2016; Ma & de la Torre, 2019; Sorrel, Abad, Olea, de la Torre, & Barrada, 2017; Sorrel, de la Torre, Abad, & Olea, 2017; see \code{\link{modelcomp}}),
#' and differential item functioning (Hou, de la Torre, & Nandakumar, 2014; Ma, Terzi, Lee, & de la Torre, 2017;
#' see \code{\link{dif}}) can also be conducted.
#'
#' To use the graphical user interface, check \code{\link{startGDINA}}.
#'
#' @name GDINA-package
#' @aliases GDINA-package
#'
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @seealso \pkg{CDM} for estimating G-DINA model and a set of other CDMs;
#' \pkg{ACTCD} and \pkg{NPCD}
#' for nonparametric CDMs; \pkg{dina} for DINA model in Bayesian framework
"_PACKAGE"
#'
#' @references
#'
#' Chen, J., & de la Torre, J. (2013). A General Cognitive Diagnosis Model for Expert-Defined Polytomous Attributes.
#' \emph{Applied Psychological Measurement, 37}, 419-437.
#'
#' Chen, J., de la Torre, J., & Zhang, Z. (2013). Relative and Absolute Fit Evaluation in Cognitive Diagnosis Modeling.
#' \emph{Journal of Educational Measurement, 50}, 123-140.
#'
#' de la Torre, J. (2009). DINA Model and Parameter Estimation: A Didactic. \emph{Journal of Educational and Behavioral Statistics, 34}, 115-130.
#'
#' de la Torre, J. (2011). The generalized DINA model framework. \emph{Psychometrika, 76}, 179-199.
#'
#' de la Torre, J. & Chiu, C-Y. (2016). A General Method of Empirical Q-matrix Validation. \emph{Psychometrika, 81}, 253-273.
#'
#' de la Torre, J., & Douglas, J. A. (2004). Higher-order latent trait models for cognitive diagnosis.
#' \emph{Psychometrika, 69}, 333-353.
#'
#' de La Torre, J., & Douglas, J. A. (2008). Model evaluation and multiple strategies in cognitive diagnosis: An analysis of fraction subtraction data. \emph{Psychometrika, 73}, 595.
#'
#' de la Torre, J., & Lee, Y. S. (2013). Evaluating the wald test for item-level comparison of
#' saturated and reduced models in cognitive diagnosis. \emph{Journal of Educational Measurement, 50}, 355-373.
#'
#' de la Torre, J., & Ma, W. (2016, August). Cognitive diagnosis modeling: A general framework approach and its implementation in R. A Short Course at the Fourth Conference on Statistical Methods in Psychometrics, Columbia University, New York.
#'
#' Haertel, E. H. (1989). Using restricted latent class models to map the skill structure of achievement items.
#' \emph{Journal of Educational Measurement, 26}, 301-321.
#'
#' Hartz, S. M. (2002). A bayesian framework for the unified model for assessing cognitive abilities:
#' Blending theory with practicality (Unpublished doctoral dissertation). University of Illinois at Urbana-Champaign.
#'
#' Hou, L., de la Torre, J., & Nandakumar, R. (2014). Differential item functioning assessment in cognitive diagnostic modeling: Application of the Wald test to
#' investigate DIF in the DINA model. \emph{Journal of Educational Measurement, 51}, 98-125.
#'
#' Junker, B. W., & Sijtsma, K. (2001). Cognitive assessment models with few assumptions, and connections with nonparametric
#' item response theory. \emph{Applied Psychological Measurement, 25}, 258-272.
#'
#' Ma, W. (2019). A diagnostic tree model for polytomous responses with multiple strategies. \emph{British Journal of Mathematical and Statistical Psychology, 72}, 61-82.
#'
#' Ma, W. (2020). Evaluating the fit of sequential G-DINA model using limited-information measures. \emph{Applied Psychological Measurement, 44}, 167-181.
#'
#' Ma, W., & de la Torre, J. (2016). A sequential cognitive diagnosis model for polytomous responses. \emph{British Journal of Mathematical and Statistical Psychology. 69,} 253-275.
#'
#' Ma, W., & de la Torre, J. (2019). Category-Level Model Selection for the Sequential G-DINA Model. \emph{Journal of Educational and Behavioral Statistics}. 44, 61-82.
#'
#' Ma, W., & de la Torre, J. (2019). Digital Module 05: Diagnostic measurement-The G-DINA framework.\emph{ Educational Measurement: Issues and Practice, 39}, 114-115.
#'
#' Ma, W., & de la Torre, J. (2020a). GDINA: An R Package for Cognitive Diagnosis Modeling. \emph{Journal of Statistical Software, 93(14)}, 1-26.
#'
#' Ma, W., & de la Torre, J. (2020b). An empirical Q-matrix validation method for the sequential G-DINA model. \emph{British Journal of Mathematical and Statistical Psychology, 73}, 142-163.
#'
#' Ma, W., & Guo, W. (2019). Cognitive diagnosis models for multiple strategies. \emph{British Journal of Mathematical and Statistical Psychology, 72}, 370-392.
#'
#' Ma, W., Iaconangelo, C., & de la Torre, J. (2016). Model similarity, model selection and attribute classification.
#' \emph{Applied Psychological Measurement, 40}, 200-217.
#'
#' Ma, W., Terzi, R., Lee, S., & de la Torre, J. (2017, April). Multiple group cognitive diagnosis models and their applications in detecting differential item functioning. Paper presented at the Annual Meeting ofthe American Educational Research Association, San Antonio, TX.
#'
#' Maris, E. (1999). Estimating multiple classification latent class models. \emph{Psychometrika, 64}, 187-212.
#'
#' Najera, P., Sorrel, M., & Abad, P. (2019). Reconsidering Cutoff Points in the General Method of Empirical Q-Matrix Validation. \emph{Educational and Psychological Measurement}.
#'
#' Sorrel, M. A., Abad, F. J., Olea, J., de la Torre, J., & Barrada, J. R. (2017). Inferential Item-Fit Evaluation in Cognitive Diagnosis Modeling. \emph{Applied Psychological Measurement, 41,} 614-631.
#'
#' Sorrel, M. A., de la Torre, J., Abad, F. J., & Olea, J. (2017). Two-Step Likelihood Ratio Test for Item-Level Model Comparison in Cognitive Diagnosis Models. \emph{Methodology, 13}, 39-47.

#' Xu, X., & von Davier, M. (2008). Fitting the structured general diagnostic model to NAEP data. ETS research report, RR-08-27.
#'
#'
#' @keywords package
NULL
