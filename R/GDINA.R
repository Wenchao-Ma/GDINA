#' @title CDM calibration under the G-DINA model framework
#'
#' @description
#' \code{GDINA} calibrates the generalized deterministic inputs, noisy and
#' gate (G-DINA; de la Torre, 2011) model for dichotomous responses, and its extension, the sequential
#' G-DINA model (Ma, & de la Torre, 2016a; Ma, 2017) for ordinal and nominal responses.
#' By setting appropriate constraints, the deterministic inputs,
#' noisy and gate (DINA; de la Torre, 2009; Junker & Sijtsma, 2001) model,
#' the deterministic inputs, noisy or gate (DINO; Templin & Henson, 2006)
#' model, the reduced reparametrized unified model (R-RUM; Hartz, 2002),
#' the additive CDM (A-CDM; de la Torre, 2011), the linear logistic
#' model (LLM; Maris, 1999), and the multiple-strategy DINA model (MSDINA; de la Torre & Douglas, 2008; Huo & de la Torre, 2014)
#' can also be calibrated. Note that the LLM is equivalent to
#' the C-RUM (Hartz, 2002), a special case of the GDM (von Davier, 2008), and that the R-RUM
#' is also known as a special case of the generalized NIDA model (de la Torre, 2011).
#'
#' In addition, users are allowed to specify design matrix and link function for each item, and
#' distinct models may be used in a single test for different items.
#' The attributes can be either dichotomous or polytomous
#' (Chen & de la Torre, 2013). Joint attribute distribution may be modelled using independent or saturated model,
#' structured model, higher-order model (de la Torre & Douglas, 2004), or loglinear model (Xu & von Davier, 2008).
#' Marginal maximum likelihood method with Expectation-Maximization (MMLE/EM) alogrithm is used for item parameter estimation.
#'
#' @section The G-DINA model:
#'
#' The generalized DINA model (G-DINA; de la Torre, 2011) is an extension of the DINA model.
#' Unlike the DINA model, which collaspes all latent classes into two latent groups for
#' each item, if item \eqn{j} requires \eqn{K_j^*}
#' attributes, the G-DINA model collapses \eqn{2^K} latent classes into \eqn{2^{K_j^*}}
#' latent groups with unique success probabilities on item \eqn{j}, where
#' \eqn{K_j^*=\sum_{k=1}^{K}q_{jk}}.
#'
#' Let \eqn{\bm{\alpha}_{lj}^*} be the reduced attribute
#' pattern consisting of the columns of the attributes required by item \eqn{j}, where
#' \eqn{l=1,\ldots,2^{K_j^*}}. For example, if only the first and the last attributes are
#' required, \eqn{\bm{\alpha}_{lj}^*=(\alpha_{l1},\alpha_{lK})}. For notational
#' convenience, the first \eqn{K_j^*} attributes can be assumed to be the required attributes
#' for item \eqn{j} as in de la Torre (2011). The probability of success \eqn{P(X_{j}=1|\bm{\alpha}_{lj}^*)} is denoted
#' by \eqn{P(\bm{\alpha}_{lj}^*)}. To model this probability of success, different link functions
#' as in the generalized linear models are used in the G-DINA model. The item response
#' function of the G-DINA model using the identity link can be written as
#' \deqn{
#' f[P(\bm{\alpha}_{lj}^*)]=\delta_{j0}+\sum_{k=1}^{K_j^*}\delta_{jk}\alpha_{lk}+
#' \sum_{k'=k+1}^{K_j^*}\sum_{k=1}^{K_j^*-1}\delta_{jkk'}\alpha_{lk}\alpha_{lk'}+\cdots+
#' \delta_{j12{\cdots}K_j^*}\prod_{k=1}^{K_j^*}\alpha_{lk},
#' }
#' or in matrix form,
#' \deqn{
#' f[\bm{P}_j]=\bm{M}_j\bm{\delta}_j,
#' }
#' where \eqn{\delta_{j0}} is the intercept for item \eqn{j}, \eqn{\delta_{jk}} is the main effect
#' due to \eqn{\alpha_{lk}}, \eqn{\delta_{jkk'}} is the interaction effect due to
#' \eqn{\alpha_{lk}} and \eqn{\alpha_{lk'}}, \eqn{\delta_{j12{\ldots}K_j^*}} is the interaction
#' effect due to \eqn{\alpha_{l1}, \cdots,\alpha_{lK_j^*}}. The log and logit links can also
#' be employed.
#'
#' @section Other CDMs as special cases:
#'
#' Several widely used CDMs can be obtained by setting appropriate constraints to the G-DINA model.
#' This section introduces the parameterization
#' of different CDMs within the G-DINA model framework very breifly. Readers interested in this please refer to
#' de la Torre(2011) for details.
#'
#' \describe{
#'   \item{\code{DINA model}}{
#'        In DINA model, each item has two item parameters - guessing (\eqn{g}) and slip (\eqn{s}). In traditional
#'        parameterization of the DINA model, a latent variable \eqn{\eta} for person \eqn{i} and
#'        item \eqn{j} is defined as
#'        \deqn{\eta_{ij}=\prod_{k=1}^K\alpha_{ik}^{q_{jk}}}
#'        Briefly speaking, if individual \eqn{i} master all attributes required by item \eqn{j},
#'        \eqn{\eta_{ij}=1}; otherwise, \eqn{\eta_{ij}=0}.
#'        Item response function of the DINA model can be written by
#'        \deqn{P(X_{ij}=1|\eta_{ij})=(1-s_j)^{\eta_{ij}}g_j^{1-\eta_{ij}}}
#'        To obtain the DINA model from the G-DINA model,
#'        all terms in identity link G-DINA model except \eqn{\delta_0} and \eqn{\delta_{12{\ldots}K_j^*}}
#'        need to be fixed to zero, that is,
#'        \deqn{ P(\bm{\alpha}_{lj}^*)=\delta_{j0}+\delta_{j12{\cdots}K_j^*}\prod_{k=1}^{K_j^*}\alpha_{lk}}
#'        In this parameterization, \eqn{\delta_{j0}=g_j} and \eqn{\delta_{j0}+\delta_{j12{\cdots}K_j^*}=1-s_j}.
#'
#'    }
#' \item{\code{DINO model}}{
#'        The DINO model can be given by
#'        \deqn{P(\bm{\alpha}_{lj}^*)=\delta_{j0}+\delta_{j1}I(\bm{\alpha}_{lj}^*\neq \bm{0})}
#'
#'        where \eqn{I(\cdot)} is an indicator variable. The DINO model is also a constrained identity
#'        link G-DINA model. As shown by de la Torre (2011), the appropriate constraint is
#'        \deqn{\delta_{jk}=-\delta_{jk^{'}k^{''}}=\cdots=(-1)^{K_j^*+1}\delta_{j12{\cdots}K_j^*},} for
#'        \eqn{k=1,\cdots,K_j^*, k^{'}=1,\cdots,K_j^*-1$, and $k^{''}>k^{'},\cdots,K_j^*}.
#'    }
#' \item{\code{Additive models with different link functions}}{
#'        The A-CDM, LLM and R-RUM can be obtained by setting all interactions to be zero in
#'        identity, logit and log link G-DINA model, respectively. Specifically, the A-CDM can be formulated as
#'        \deqn{P(\bm{\alpha}_{lj}^*)=\delta_{j0}+\sum_{k=1}^{K_j^*}\delta_{jk}\alpha_{lk}.}
#'        The item response function for
#'        LLM can be given by
#'        \deqn{ logit[P(\bm{\alpha}_{lj}^*)]=\delta_{j0}+\sum_{k=1}^{K_j^*}\delta_{jk}\alpha_{lk},}
#'        and lastly, the RRUM, can be written as
#'        \deqn{log[P(\bm{\alpha}_{lj}^*)]=\delta_{j0}+\sum_{k=1}^{K_j^*}\delta_{jk}\alpha_{lk}.} It should be
#'        noted that the LLM is equivalent to the compensatory RUM, which is subsumed by the GDM, and that
#'        the RRUM is a special case of the generalized noisy inputs, deterministic ``And" gate model (G-NIDA).
#'    }
#'    }
#'
#' @section Joint Attribute Distribution:
#'
#' The joint attribute distribution can be modeled using various methods. This section mainly focuses on the so-called
#' higher-order approach, which was originally proposed by de la Torre
#' and Douglas (2004) for the DINA model. It has been extended in this package for all condensation rules.
#' Particularly, three IRT models are available for the higher-order attribute structure:
#' Rasch model (Rasch), one parameter logistic model (1PL) and two parameter logistic model (2PL).
#' For the Rasch model, the probability of mastering attribute \eqn{k} for individual \eqn{i} is defined as
#' \deqn{P(\alpha_k=1|\theta_i,\lambda_{0k})=\frac{exp(\theta_i+\lambda_{0k})}{1+exp(\theta_i+\lambda_{0k})}}
#' For the 1PL model, the probability of mastering attribute \eqn{k} for individual \eqn{i} is defined as
#' \deqn{P(\alpha_k=1|\theta_i,\lambda_{0k},\lambda_{1})=\frac{exp(\lambda_{1}\theta_i+\lambda_{0k})}{1+exp(\lambda_{1}\theta_i+\lambda_{0k})}}
#' For the 2PL model, the probability of mastering attribute \eqn{k} for individual \eqn{i} is defined as
#' \deqn{P(\alpha_k=1|\theta_i,\lambda_{0k},\lambda_{1k})=\frac{exp(\lambda_{1k}\theta_i+\lambda_{0k})}{1+exp(\lambda_{1k}\theta_i+\lambda_{0k})}}
#' where \eqn{\theta_i} is the ability of examinee \eqn{i}. \eqn{\lambda_{0k}} and \eqn{\lambda_{1k}} are the intercept
#' and slope parameters for attribute \eqn{k}, respectively. In the Rasch model, \eqn{\lambda_{1k}=1 \forall k};
#' whereas in the 1PL model, a common slope parameter \eqn{\lambda_{1}} is estimated.
#' The probability of joint attributes can be written as
#'  \deqn{P(\strong{\alpha}|\theta_i,\strong{\lambda})=\prod_k P(\alpha_k|\theta_i,\strong{\lambda})}.
#'
#'
#'@section Model Estimation:
#'
#' The MMLE/EM algorithm is implemented in this package. For G-DINA, DINA and DINO models, closed-form solutions exist.
#' See de la Torre (2009) and de la Torre (2011) for details.
#' For ACDM, LLM and RRUM, closed-form solutions do not exist, and therefore some general optimization techniques are
#' adopted in M-step (Ma, Iaconangelo & de la Torre, 2016). The selection of optimization techniques mainly depends on whether
#' some specific constraints need to be added.
#'
#' The sequential G-DINA model is a special case of the diagnostic tree model (DTM; Ma, in press) and estimated using
#' the mapping matrix accordingly (See Tutz, 1997; Ma, in press).
#'
#'
#' @section The Number of Parameters:
#'
#' For dichotomous response models:
#' Assume a test measures \eqn{K} attributes and item \eqn{j} requires \eqn{K_j^*} attributes:
#' The DINA and DINO model has 2 item parameters for each item;
#' if item \eqn{j} is ACDM, LLM or RRUM, it has \eqn{K_j^*+1} item parameters; if it is G-DINA model, it has \eqn{2^{K_j^*}} item parameters.
#' Apart from item parameters, the parameters involved in the estimation of joint attribute distribution need to be estimated as well.
#' When using the saturated attribute structure, there are \eqn{2^K-1} parameters for joint attribute distribution estimation; when
#' using a higher-order attribute structure, there are \eqn{K}, \eqn{K+1}, and \eqn{2\times K} parameters for the Rasch model,
#' 1PL model and 2PL model, respectively.
#' For polytomous response data using the sequential G-DINA model, the number of item parameters
#' are counted at category level.
#'
#' @param dat A required \eqn{N \times J} \code{matrix} or \code{data.frame} consisting of the
#'     responses of \eqn{N} individuals to \eqn{J} items. Missing values need to be coded as \code{NA}.
#' @param Q A required matrix; The number of rows occupied by a single-strategy dichotomous item is 1, by a polytomous item is
#' the number of nonzero categories, and by a mutiple-strategy dichotomous item is the number of strategies.
#' The number of column is equal to the number of attributes if all items are single-strategy dichotomous items, but
#' the number of attributes + 2 if any items are polytomous or have multiple strategies.
#' For a polytomous item, the first column represents the item number and the second column indicates the nonzero category number.
#' For a multiple-strategy dichotomous item, the first column represents the item number and the second column indicates the strategy number.
#' For binary attributes, 1 denotes the attributes are measured by the items and 0 means the attributes are not
#'    measured. For polytomous attributes, non-zero elements indicate which level
#'    of attributes are needed (see Chen, & de la Torre, 2013).  See \code{Examples}.
#' @param model A vector for each item or nonzero category, or a scalar which will be used for all
#'    items or nonzero categories to specify the CDMs fitted. The possible options
#'    include \code{"GDINA"},\code{"DINA"},\code{"DINO"},\code{"ACDM"},\code{"LLM"}, \code{"RRUM"}, \code{"MSDINA"} and \code{"UDF"}.
#'    When \code{"UDF"}, indicating user defined function, is specified for any item, arguments \code{design.matrix} and \code{linkfunc} need to be defined.
#' @param sequential logical; \code{TRUE} if the sequential model is fitted for polytomous responses.
#' @param group a numerical vector with integer 1, 2, ..., # of groups indicating the group each individual belongs to. It must start from 1 and its
#'    length must be equal to the number of individuals.
#' @param att.dist How is the joint attribute distribution estimated? It can be (1) \code{saturated}, which is the default, indicating that
#'    the proportion parameter for each permissible latent class is estimated separately; (2) \code{higher.order}, indicating
#'    that a higher-order joint attribute distribution is assumed (higher-order model can be specified in \code{higher.order} argument);
#'    (3) \code{fixed}, indicating that the weights specified in \code{att.prior} argument are fixed in the estimation process.
#'    If \code{att.prior} is not specified, a uniform joint attribute distribution is employed initially;  (4) \code{independent}, indicating
#'    that all attributes are assumed to be independent; and (5) \code{loglinear}, indicating a loglinear model is employed.
#'    If different groups have different joint attribute distributions,
#'    specify \code{att.dist} as a character vector with the same number of elements as the number of groups. However, if a higher-order model is used for any group,
#'    it must be used for all groups.
#' @param mono.constraint logical; \code{TRUE} indicates that \eqn{P(\bm{\alpha}_1) <=P(\bm{\alpha}_2)} if
#'    for all \eqn{k}, \eqn{\alpha_{1k} < \alpha_{2k}}. Can be a vector for each item or nonzero category or a scalar which will be used for all
#'    items to specify whether monotonicity constraint should be added.
#' @param catprob.parm A list of initial success probability parameters for each nonzero category.
#' @param item.names A vector giving the item names. By default, items are named as "Item 1", "Item 2", etc.
#' @param higher.order A list specifying the higher-order joint attribute distribution with the following components:
#'  \itemize{
#'    \item \code{model} - a character indicating the IRT model for higher-order joint attribute distribution. Can be
#'    \code{"2PL"}, \code{"1PL"} or \code{"Rasch"}, representing two parameter logistic IRT model,
#'    one parameter logistic IRT model and Rasch model, respectively. For \code{"1PL"} model, a common slope parameter is
#'    estimated. \code{"Rasch"} is the default model when \code{att.dist = "higher.order"}. Note that slope-intercept form
#'    is used for parameterizing the higher-order IRT model (see \code{Details}).
#'    \item \code{nquad} - a scalar specifying the number of integral nodes. Default = 25.
#'    \item \code{SlopeRange} - a vector of length two specifying the range of slope parameters. Default = [0.1, 5].
#'    \item \code{InterceptRange} - a vector of length two specifying the range of intercept parameters. Default = [-4, 4].
#'    \item \code{SlopePrior} - a vector of length two specifying the mean and variance of log(slope) parameters, which are assumed normally distributed. Default: mean = 0 and sd = 0.25.
#'    \item \code{InterceptPrior} - a vector of length two specifying the mean and variance of intercept parameters, which are assumed normally distributed. Default: mean = 0 and sd = 1.
#'    \item \code{Prior} - logical; indicating whether prior distributions should be imposed to slope and intercept parameters. Default is \code{FALSE}.
#'    }
#' @param verbose How to print calibration information
#'     after each EM iteration? Can be 0, 1 or 2, indicating to print no information,
#'     information for current iteration, or information for all iterations.
#' @param att.prior A vector of length \eqn{2^K} for single group model, or a matrix of dimension \eqn{2^K\times} no. of groups to specify
#'    attribute prior distribution for \eqn{2^K} latent classes for all groups under a multiple group model. Only applicable for dichotomous attributes.
#'    The sum of all elements does not have to be equal to 1; however, it will be normalized so that the sum is equal to 1
#'    before calibration.
#'    The label for each latent class can be obtained by calling \code{attributepattern(K)}. See \code{examples} for more info.
#' @param latent.var A string indicating the nature of the latent variables. It is \code{"att"} (by default) if the latent variables are attributes,
#'    and \code{"bugs"} if the latent variables are misconceptions. When \code{"bugs"} is specified, only the DINA, DINO or G-DINA model can be
#'    specified in \code{model} argument (Kuo, Chen, Yang & Mok, 2016).
#' @param att.str logical; are attributes structured? If yes, \code{att.prior} must be specified where impossible latent classes have prior weights 0.
#'    If attributes are structured, only the DINA, DINO or G-DINA model can be specified in \code{model} argument.
#' @param loglinear the order of loglinear smooth for attribute space. It can be either 1 or 2 indicating the loglinear model with main effect only
#'    and with main effect and first-order interaction.
#' @param control A list of control parameters with elements:
#' \itemize{
#'      \item \code{maxitr} A vector for each item or nonzero category, or a scalar which will be used for all
#'    items or nonzero categories to specify the maximum number of EM cycles allowed. Default = 2000.
#'     \item \code{conv.crit} The convergence criterion for max absolute change in item parameters or deviance. Default = 0.0001.
#'     \item \code{conv.type} How is the convergence criterion evaluated? A vector with possible elements: \code{"ip"}, indicating
#'    the maximum absolute change in item success probabilities, \code{"mp"}, representing
#'    the maximum absolute change in mixing proportion parameters, \code{"delta"}, indicating the maximum absolute change in delta
#'    parameters or \code{neg2LL} indicating the absolute change in negative two times loglikeihood. Multiple criteria can be specified.
#'    If so, all criteria need to be met. Default = c("ip", "mp").
#'     \item \code{nstarts} how many sets of starting values? Default = 3.
#'     \item \code{lower.p} A vector for each item or nonzero category,
#'    or a scalar which will be used for all items or nonzero categories to specify the lower bound for success probabilities.
#'    Default = .0001.
#'     \item \code{upper.p} A vector for each item or nonzero category, or a scalar which will be used for all
#'    items or nonzero categories to specify the upper bound for success probabilities. Default = .9999.
#'     \item \code{lower.prior} The lower bound for mixing proportion parameters (latent class sizes). Default = .Machine$double.eps.
#'     \item \code{randomseed} Random seed for generating initial item parameters. Default = 123456.
#'     \item \code{smallNcorrection} A numeric vector with two elements specifying the corrections applied when the expected number of
#' individuals in some latent groups are too small. If the expected no. of examinees is less than the second element,
#' the first element and two times the first element will be added to the numerator and denominator of the closed-form solution of
#' probabilities of success. Only applicable for the G-DINA, DINA and DINO model estimation without monotonic constraints.
#'     \item \code{MstepMessage} Integer; Larger number prints more information from Mstep optimizer. Default = 1.
#'  }
#' @param linkfunc a vector of link functions for each item/category; It can be \code{"identity"},\code{"log"} or \code{"logit"}. Only applicable
#'    when, for some items, \code{model="UDF"}.
#' @param design.matrix a list of design matrices; Its length must be equal to the number of items (or nonzero categories for sequential models).
#'    If CDM for item j is specified as "UDF" in argument \code{model}, the corresponding design matrix must be provided; otherwise, the design matrix can be \code{NULL},
#'    which will be generated automatically.
#' @param solver A string indicating which solver should be used in M-step. By default, the solver is automatically chosen according to the models specified.
#'    Possible options include \link[nloptr]{slsqp}, \link[nloptr]{nloptr}, \link[Rsolnp]{solnp} and \link[alabama]{auglag}.
#' @param auglag.args a list of control parameters to be passed to the alabama::auglag() function. It can contain two elements:
#'    \code{control.outer} and \code{control.optim}. See \link[alabama]{auglag}.
#' @param nloptr.args a list of control parameters to be passed to \code{opts} argument of \link[nloptr]{nloptr} function.
#' @param solnp.args  a list of control parameters to be passed to \code{control} argument of \link[Rsolnp]{solnp} function.
#' @param ... additional arguments
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @seealso See \code{\link{autoGDINA}} for Q-matrix validation, item-level model comparison and model calibration
#' in one run; See \code{\link{modelfit}} and \code{\link{itemfit}} for model and item fit analysis, \code{\link{Qval}} for Q-matrix validation,
#' \code{\link{modelcomp}} for item level model comparison and \code{\link{simGDINA}} for data simulation.
#' Also see \code{gdina} in \pkg{CDM} package for the G-DINA model estimation.
#'
#' @return \code{GDINA} returns an object of class \code{GDINA}. Methods for \code{GDINA} objects
#'  include \code{\link{extract}} for extracting various components, \code{\link{coef}}
#'  for extracting structural parameters, \code{\link{personparm}}
#'  for calculating incidental (person) parameters, \code{summary} for summary information.
#'  \code{AIC}, \code{BIC},\code{logLik}, \code{deviance} and \code{npar} can also be used to
#'  calculate AIC, BIC, observed log-likelihood, deviance and number of parameters.
#'
#'
#'
#' @references
#' Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of item parameters: Application of an EM algorithm. \emph{Psychometrika, 46}, 443-459.
#'
#' Bock, R. D., & Lieberman, M. (1970). Fitting a response model forn dichotomously scored items. \emph{Psychometrika, 35}, 179-197.
#'
#' Bor-Chen Kuo, Chun-Hua Chen, Chih-Wei Yang, & Magdalena Mo Ching Mok. (2016). Cognitive diagnostic models for tests with multiple-choice and constructed-response items. \emph{Educational Psychology, 36}, 1115-1133.
#'
#' Carlin, B. P., & Louis, T. A. (2000). Bayes and empirical bayes methods for data analysis. New York, NY: Chapman & Hall
#'
#' de la Torre, J., & Douglas, J. A. (2008). Model evaluation and multiple strategies in cognitive diagnosis: An analysis of fraction subtraction data. \emph{Psychometrika, 73}, 595-624.
#'
#' de la Torre, J. (2009). DINA Model and Parameter Estimation: A Didactic. \emph{Journal of Educational and Behavioral Statistics, 34}, 115-130.
#'
#' de la Torre, J. (2011). The generalized DINA model framework. \emph{Psychometrika, 76}, 179-199.
#'
#' de la Torre, J., & Douglas, J. A. (2004). Higher-order latent trait models for cognitive diagnosis. \emph{Psychometrika, 69}, 333-353.
#'
#' de la Torre, J., & Lee, Y. S. (2013). Evaluating the wald test for item-level comparison of
#' saturated and reduced models in cognitive diagnosis. \emph{Journal of Educational Measurement, 50}, 355-373.
#'
#' Haertel, E. H. (1989). Using restricted latent class models to map the skill structure of achievement items.
#' \emph{Journal of Educational Measurement, 26}, 301-321.
#'
#' Hartz, S. M. (2002). A bayesian framework for the unified model for assessing cognitive abilities:
#' Blending theory with practicality (Unpublished doctoral dissertation). University of Illinois at Urbana-Champaign.
#'
#' Huo, Y., & de la Torre, J. (2014). Estimating a Cognitive Diagnostic Model for Multiple Strategies via the EM Algorithm. \emph{Applied Psychological Measurement, 38}, 464-485.
#'
#' Junker, B. W., & Sijtsma, K. (2001). Cognitive assessment models with few assumptions, and connections with nonparametric
#' item response theory. \emph{Applied Psychological Measurement, 25}, 258-272.
#'
#' Ma, W., & de la Torre, J. (2016). A sequential cognitive diagnosis model for polytomous responses. \emph{British Journal of Mathematical and Statistical Psychology. 69,} 253-275.
#'
#' Ma, W. (2018). A Diagnostic Tree Model for Polytomous Responses with Multiple Strategies. \emph{British Journal of Mathematical and Statistical Psychology.}
#'
#' Ma, W., Iaconangelo, C., & de la Torre, J. (2016). Model similarity, model selection and attribute classification. \emph{Applied Psychological Measurement, 40}, 200-217.
#'
#' Ma, W. (2017). \emph{A Sequential Cognitive Diagnosis Model for Graded Response: Model Development, Q-Matrix Validation,and Model Comparison. Unpublished doctoral dissertation.} New Brunswick, NJ: Rutgers University.
#'
#' Maris, E. (1999). Estimating multiple classification latent class models. \emph{Psychometrika, 64}, 187-212.
#'
#' Tatsuoka, K. K. (1983). Rule space: An approach for dealing with misconceptions based on item response theory. \emph{Journal of Educational Measurement, 20}, 345-354.
#'
#' Templin, J. L., & Henson, R. A. (2006). Measurement of psychological disorders using cognitive diagnosis models. \emph{Psychological Methods, 11}, 287-305.
#'
#' Tutz, G. (1997). Sequential models for ordered responses. In W.J. van der Linden & R. K. Hambleton (Eds.), Handbook of modern item response theory p. 139-152). New York, NY: Springer.
#'
#' Xu, X., & von Davier, M. (2008). Fitting the structured general diagnostic model to NAEP data. ETS research report, RR-08-27.
#'
#' @importFrom nloptr nloptr slsqp nl.grad nl.jacobian
#' @importFrom Rcpp sourceCpp
#' @importFrom numDeriv hessian jacobian
#' @importFrom alabama auglag
#' @importFrom MASS ginv
#' @importFrom Rsolnp solnp
#' @importFrom graphics plot axis points text abline
#' @importFrom utils combn modifyList
#' @import stats
#' @importFrom ggplot2 aes ggplot geom_tile scale_fill_gradient theme_bw labs geom_bar geom_errorbar ylim
#' @import shiny
#' @import shinydashboard
#'
#' @useDynLib GDINA
#' @export
#' @examples
#' \dontrun{
#'####################################
#'#        Example 1.                #
#'#     GDINA, DINA, DINO            #
#'#    ACDM, LLM and RRUM            #
#'# estimation and comparison        #
#'#                                  #
#'####################################
#'
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#'
#' #--------GDINA model --------#
#'
#' mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' mod1
#' # summary information
#' summary(mod1)
#'
#' AIC(mod1) #AIC
#' BIC(mod1) #BIC
#' logLik(mod1) #log-likelihood value
#' deviance(mod1) # deviance: -2 log-likelihood
#' npar(mod1) # number of parameters
#'
#'
#' head(indlogLik(mod1)) # individual log-likelihood
#' head(indlogPost(mod1)) # individual log-posterior
#'
#' # structural parameters
#' # see ?coef
#' coef(mod1) # item probabilities of success for each latent group
#' coef(mod1, withSE = TRUE) # item probabilities of success & standard errors
#' coef(mod1, what = "delta") # delta parameters
#' coef(mod1, what = "delta",withSE=TRUE) # delta parameters
#' coef(mod1, what = "gs") # guessing and slip parameters
#' coef(mod1, what = "gs",withSE = TRUE) # guessing and slip parameters & standard errors
#'
#' # person parameters
#' # see ?personparm
#' personparm(mod1) # EAP estimates of attribute profiles
#' personparm(mod1, what = "MAP") # MAP estimates of attribute profiles
#' personparm(mod1, what = "MLE") # MLE estimates of attribute profiles
#'
#' #plot item response functions for item 10
#' plot(mod1,item = 10)
#' plot(mod1,item = 10,withSE = TRUE) # with error bars
#' #plot mastery probability for individuals 1, 20 and 50
#' plot(mod1,what = "mp", person =c(1,20,50))
#'
#' # Use extract function to extract more components
#' # See ?extract
#'
#' # ------- DINA model --------#
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' mod2 <- GDINA(dat = dat, Q = Q, model = "DINA")
#' mod2
#' coef(mod2, what = "gs") # guess and slip parameters
#' coef(mod2, what = "gs",withSE = TRUE) # guess and slip parameters and standard errors
#'
#' # Model comparison at the test level via likelihood ratio test
#' anova(mod1,mod2)
#'
#' # -------- DINO model -------#
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' mod3 <- GDINA(dat = dat, Q = Q, model = "DINO")
#' #slip and guessing
#' coef(mod3, what = "gs") # guess and slip parameters
#' coef(mod3, what = "gs",withSE = TRUE) # guess and slip parameters + standard errors
#'
#' # Model comparison at test level via likelihood ratio test
#' anova(mod1,mod2,mod3)
#'
#' # --------- ACDM model -------#
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' mod4 <- GDINA(dat = dat, Q = Q, model = "ACDM")
#' mod4
#' # --------- LLM model -------#
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' mod4b <- GDINA(dat = dat, Q = Q, model = "LLM")
#' mod4b
#' # --------- RRUM model -------#
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' mod4c <- GDINA(dat = dat, Q = Q, model = "RRUM")
#' mod4c
#'
#' # --- Different CDMs for different items --- #
#'
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' models <- c(rep("GDINA",3),"LLM","DINA","DINO","ACDM","RRUM","LLM","RRUM")
#' mod5 <- GDINA(dat = dat, Q = Q, model = models)
#' anova(mod1,mod2,mod3,mod4,mod4b,mod4c,mod5)
#'
#'
#'####################################
#'#        Example 2.                #
#'#        Model estimations         #
#'# With monotonocity constraints    #
#'####################################
#'
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' # for item 10 only
#' mod11 <- GDINA(dat = dat, Q = Q, model = "GDINA",mono.constraint = c(rep(FALSE,9),TRUE))
#' mod11
#' mod11a <- GDINA(dat = dat, Q = Q, model = "DINA",mono.constraint = TRUE)
#' mod11a
#' mod11b <- GDINA(dat = dat, Q = Q, model = "ACDM",mono.constraint = TRUE)
#' mod11b
#' mod11c <- GDINA(dat = dat, Q = Q, model = "LLM",mono.constraint = TRUE)
#' mod11c
#' mod11d <- GDINA(dat = dat, Q = Q, model = "RRUM",mono.constraint = TRUE)
#' mod11d
#' coef(mod11d,"delta")
#' coef(mod11d,"rrum")
#'
#'####################################
#'#           Example 3a.            #
#'#        Model estimations         #
#'# With Higher-order att structure  #
#'####################################
#'
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' # --- Higher order G-DINA model ---#
#' mod12 <- GDINA(dat = dat, Q = Q, model = "DINA",
#'                att.dist="higher.order",higher.order=list(nquad=31,model = "2PL"))
#' personparm(mod12,"HO") # higher-order ability
#' # structural parameters
#' # first column is slope and the second column is intercept
#' coef(mod12,"lambda")
#' # --- Higher order DINA model ---#
#' mod22 <- GDINA(dat = dat, Q = Q, model = "DINA", att.dist="higher.order",
#'                higher.order=list(model = "2PL",Prior=TRUE))
#'
#'####################################
#'#           Example 3b.            #
#'#        Model estimations         #
#'#   With log-linear att structure  #
#'####################################
#'
#' # --- DINA model with loglinear smoothed attribute space ---#
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' mod23 <- GDINA(dat = dat, Q = Q, model = "DINA",att.dist="loglinear",loglinear=1)
#' coef(mod23,"lambda") # intercept and three main effects
#'
#'####################################
#'#           Example 3c.            #
#'#        Model estimations         #
#'#  With independent att structure  #
#'####################################
#'
#' # --- GDINA model with independent attribute space ---#
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' mod33 <- GDINA(dat = dat, Q = Q, att.dist="independent")
#' coef(mod33,"lambda") # mastery probability for each attribute
#'
#'####################################
#'#          Example 4.              #
#'#        Model estimations         #
#'# With user-specified att structure#
#'####################################
#'
#' # --- User-specified attribute priors ----#
#' # prior distribution is fixed during calibration
#' # Assume each of 000,100,010 and 001 has probability of 0.1
#' # and each of 110, 101,011 and 111 has probability of 0.15
#' # Note that the sum is equal to 1
#' #
#' prior <- c(0.1,0.1,0.1,0.1,0.15,0.15,0.15,0.15)
#' # fit GDINA model  with fixed prior dist.
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' modp1 <- GDINA(dat = dat, Q = Q, att.prior = prior, att.dist = "fixed")
#' # See the posterior weights
#' extract(modp1,what = "posterior.prob")
#' extract(modp1,what = "att.prior")
#' # ----Linear structure of attributes -----#
#' # Assuming A1 -> A2 -> A3
#'   Q <- matrix(c(1,0,0,
#'                 1,0,0,
#'                 1,1,0,
#'                 1,1,0,
#'                 1,1,1,
#'                 1,1,1,
#'                 1,0,0,
#'                 1,0,0,
#'                 1,1,0,
#'                 1,1,0,
#'                 1,1,1,
#'                 1,1,1),ncol=3,byrow=TRUE)
#'  # item parameters for DINA model (guessing and slip)
#'  gs <- matrix(rep(0.1,24),ncol=2)
#'  N <- 5000
#'  # attribute simulation
#'  att <- rbind(matrix(0,nrow=500,ncol=3),
#'               matrix(rep(c(1,0,0),1000),ncol=3,byrow=TRUE),
#'               matrix(rep(c(1,1,0),1000),ncol=3,byrow=TRUE),
#'               matrix(rep(c(1,1,1),2500),ncol=3,byrow=TRUE))
#'  # data simulation
#'  simD <- simGDINA(N,Q,gs.parm = gs, model = "DINA",attribute = att)
#'  dat <- simD$dat
#'  # setting structure: A1 -> A2 -> A3
#'  # note: latent classes with prior 0 are assumed impossible
#'  prior <- c(0.1,0.2,0,0,0.2,0,0,0.5)
#'  out <- GDINA(dat, Q, att.prior = prior,att.str = TRUE, att.dist = "fixed", model = "DINA")
#'  # check posterior dist.
#'  extract(out,what = "posterior.prob")
#'  extract(out,what = "att.prior")
#'
#'  out2 <- GDINA(dat, Q, att.prior = prior,att.str = TRUE, att.dist = "saturated",model = "DINA")
#'  # check posterior dist.
#'  extract(out2,what = "posterior.prob")
#'  extract(out2,what = "att.prior")
#'
#'####################################
#'#          Example 5.              #
#'#        Model estimations         #
#'# With user-specified att structure#
#'####################################
#'
#' # --- User-specified attribute structure ----#
#' Q <- sim30GDINA$simQ
#' K <- ncol(Q)
#' # divergent structure A1->A2->A3;A1->A4->A5
#' diverg <- list(c(1,2),
#'                c(2,3),
#'                c(1,4),
#'                c(4,5))
#'struc <- att.structure(diverg,K)
#'
#' # data simulation
#' N <- 1000
#' true.lc <- sample(c(1:2^K),N,replace=TRUE,prob=struc$att.prob)
#' table(true.lc) #check the sample
#' true.att <- attributepattern(K)[true.lc,]
#'  gs <- matrix(rep(0.1,2*nrow(Q)),ncol=2)
#'  # data simulation
#'  simD <- simGDINA(N,Q,gs.parm = gs, model = "DINA",attribute = true.att)
#'  dat <- extract(simD,"dat")
#'
#' modp1 <- GDINA(dat = dat, Q = Q, att.prior = struc$att.prob,
#'                att.str = TRUE, att.dist = "saturated")
#' modp1
#' # prior dist.
#' extract(modp1,what = "att.prior")
#' # Posterior weights were slightly different
#' extract(modp1,what = "posterior.prob")
#' modp2 <- GDINA(dat = dat, Q = Q, att.prior = struc$att.prob,
#'                att.str = TRUE, att.dist = "fixed")
#' modp2
#' extract(modp2,what = "att.prior")
#' extract(modp2,what = "posterior.prob")
#'
#'
#'####################################
#'#           Example 6.             #
#'#        Model estimations         #
#'# With user-specified initial pars #
#'####################################
#'
#'  # check initials to see the format for initial item parameters
#'  initials <- sim10GDINA$simItempar
#'  dat <- sim10GDINA$simdat
#'  Q <- sim10GDINA$simQ
#'  mod.ini <- GDINA(dat,Q,catprob.parm = initials)
#'  extract(mod.ini,"initial.catprob")
#'
#'####################################
#'#           Example 7.             #
#'#        Model estimation          #
#'#          Without M-step          #
#'####################################
#'
#'  # -----------Fix User-specified item parameters
#'  # Item parameters are not estimated
#'  # Only person attributes are estimated
#'  # attribute prior distribution matters if interested in the marginalized likelihood
#'  dat <- frac20$dat
#'  Q <- frac20$Q
#'  # estimation- only 20 iterations for illustration purposes
#'  mod.initial <- GDINA(dat,Q,control = list(maxitr=20))
#'  par <- coef(mod.initial,digits=8)
#'  weights <- extract(mod.initial,"posterior.prob",digits=8) #posterior weights
#'  # use the weights as the priors
#'  mod.fix <- GDINA(dat,Q,catprob.parm = par,
#'                   att.prior=c(weights),control = list(maxitr = 0)) # re-estimation
#'  anova(mod.initial,mod.fix) # very similar - good approximation most of time
#'  # prior used for the likelihood calculation for the last step
#'  priors <- extract(mod.initial,"att.prior")
#'  # use the priors as the priors
#'  mod.fix2 <- GDINA(dat,Q,catprob.parm = par,
#'                    att.prior=priors, control = list(maxitr=0)) # re-estimation
#'  anova(mod.initial,mod.fix2) # identical results
#'
#'####################################
#'#           Example 8.             #
#'#        polytomous attribute      #
#'#          model estimation        #
#'#    see Chen, de la Torre 2013    #
#'####################################
#'
#'
#' # --- polytomous attribute G-DINA model --- #
#' dat <- sim30pGDINA$simdat
#' Q <- sim30pGDINA$simQ
#' #polytomous G-DINA model
#' pout <- GDINA(dat,Q)
#'
#' # ----- polymous DINA model --------#
#' pout2 <- GDINA(dat,Q,model="DINA")
#' anova(pout,pout2)
#'
#'####################################
#'#           Example 9.             #
#'#        Sequential G-DINA model   #
#'#    see Ma, & de la Torre 2016    #
#'####################################
#'
#' # --- polytomous attribute G-DINA model --- #
#' dat <- sim20seqGDINA$simdat
#' Q <- sim20seqGDINA$simQ
#' Q
#'#    Item Cat A1 A2 A3 A4 A5
#'#       1   1  1  0  0  0  0
#'#       1   2  0  1  0  0  0
#'#       2   1  0  0  1  0  0
#'#       2   2  0  0  0  1  0
#'#       3   1  0  0  0  0  1
#'#       3   2  1  0  0  0  0
#'#       4   1  0  0  0  0  1
#'#       ...
#'
#' #sequential G-DINA model
#' sGDINA <- GDINA(dat,Q,sequential = TRUE)
#' sDINA <- GDINA(dat,Q,sequential = TRUE,model = "DINA")
#' anova(sGDINA,sDINA)
#' coef(sDINA) # processing function
#' coef(sDINA,"itemprob") # success probabilities for each item
#' coef(sDINA,"LCprob") # success probabilities for each category for all latent classes
#'
#'####################################
#'#           Example 10a.           #
#'#    Multiple-Group G-DINA model   #
#'####################################
#'
#' Q <- sim10GDINA$simQ
#' K <- ncol(Q)
#' # parameter simulation
#' # Group 1 - female
#' N1 <- 3000
#' gs1 <- matrix(rep(0.1,2*nrow(Q)),ncol=2)
#' # Group 2 - male
#' N2 <- 3000
#' gs2 <- matrix(rep(0.2,2*nrow(Q)),ncol=2)
#'
#' # data simulation for each group
#' sim1 <- simGDINA(N1,Q,gs.parm = gs1,model = "DINA",att.dist = "higher.order",
#'                  higher.order.parm = list(theta = rnorm(N1),
#'                  lambda = data.frame(a=rep(1.5,K),b=seq(-1,1,length.out=K))))
#' sim2 <- simGDINA(N2,Q,gs.parm = gs2,model = "DINO",att.dist = "higher.order",
#'                  higher.order.parm = list(theta = rnorm(N2),
#'                  lambda = data.frame(a=rep(1,K),b=seq(-2,2,length.out=K))))
#'
#' # combine data
#' # see ?bdiagMatrix
#' dat <- bdiagMatrix(list(extract(sim1,"dat"),extract(sim2,"dat")),fill=NA)
#' Q <- rbind(Q,Q)
#' gr <- rep(c(1,2),c(3000,3000))
#' # Fit G-DINA model
#' mg.est <- GDINA(dat = dat,Q = Q,group = gr,att.dist="higher.order",
#' higher.order=list(model = "1PL",Prior=TRUE))
#' summary(mg.est)
#' extract(mg.est,"posterior.prob")
#'
#'
#'
#'####################################
#'#           Example 10b.           #
#'#    Multiple-Group G-DINA model   #
#'####################################
#'
#' Q <- sim30GDINA$simQ
#' K <- ncol(Q)
#' # parameter simulation
#' N1 <- 3000
#' gs1 <- matrix(rep(0.1,2*nrow(Q)),ncol=2)
#' N2 <- 3000
#' gs2 <- matrix(rep(0.2,2*nrow(Q)),ncol=2)
#'
#' # data simulation for each group
#' # two groups have different theta distributions
#' sim1 <- simGDINA(N1,Q,gs.parm = gs1,model = "DINA",att.dist = "higher.order",
#'                  higher.order.parm = list(theta = rnorm(N1),
#'                  lambda = data.frame(a=rep(1,K),b=seq(-2,2,length.out=K))))
#' sim2 <- simGDINA(N2,Q,gs.parm = gs2,model = "DINO",att.dist = "higher.order",
#'                  higher.order.parm = list(theta = rnorm(N2,1,1),
#'                  lambda = data.frame(a=rep(1,K),b=seq(-2,2,length.out=K))))
#'
#' # combine data
#' # see ?bdiagMatrix
#' dat <- bdiagMatrix(list(extract(sim1,"dat"),extract(sim2,"dat")),fill=NA)
#' Q <- rbind(Q,Q)
#' gr <- rep(c(1,2),c(3000,3000))
#' # Fit G-DINA model
#' mg.est <- GDINA(dat = dat,Q = Q,group = gr,att.dist="higher.order",
#' higher.order=list(model = "Rasch",Prior=FALSE,Type = "SameLambda"))
#' summary(mg.est)
#' coef(mg.est,"lambda")
#'
#'
#'####################################
#'#           Example 11.            #
#'#           Bug DINO model         #
#'####################################
#'
#' set.seed(123)
#' Q <- sim10GDINA$simQ # 1 represents misconceptions/bugs
#' ip <- list(
#' c(0.8,0.2),
#' c(0.7,0.1),
#' c(0.9,0.2),
#' c(0.9,0.1,0.1,0.1),
#' c(0.9,0.1,0.1,0.1),
#' c(0.9,0.1,0.1,0.1),
#' c(0.9,0.1,0.1,0.1),
#' c(0.9,0.1,0.1,0.1),
#' c(0.9,0.1,0.1,0.1),
#' c(0.9,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
#' sim <- simGDINA(N=1000,Q=Q,catprob.parm = ip)
#' dat <- extract(sim,"dat")
#' # use latent.var to specify a bug model
#' est <- GDINA(dat=dat,Q=Q,latent.var="bugs",model="DINO")
#' coef(est)
#'
#'####################################
#'#           Example 12.            #
#'#           Bug DINA model         #
#'####################################
#'
#' set.seed(123)
#' Q <- sim10GDINA$simQ # 1 represents misconceptions/bugs
#' ip <- list(
#' c(0.8,0.2),
#' c(0.7,0.1),
#' c(0.9,0.2),
#' c(0.9,0.9,0.9,0.1),
#' c(0.9,0.9,0.9,0.1),
#' c(0.9,0.9,0.9,0.1),
#' c(0.9,0.9,0.9,0.1),
#' c(0.9,0.9,0.9,0.1),
#' c(0.9,0.9,0.9,0.1),
#' c(0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.1))
#' sim <- simGDINA(N=1000,Q=Q,catprob.parm = ip)
#' dat <- extract(sim,"dat")
#' # use latent.var to specify a bug model
#' est <- GDINA(dat=dat,Q=Q,latent.var="bugs",model="DINA")
#' coef(est)
#'
#'####################################
#'#           Example 13a.           #
#'#     user specified design matrix #
#'#        LCDM (logit G-DINA)       #
#'####################################
#'
#' dat <- sim30GDINA$simdat
#' Q <- sim30GDINA$simQ
#'
#' #find design matrix for each item => must be a list
#' D <- lapply(rowSums(Q),designmatrix,model="GDINA")
#' # for comparison, use change in -2LL as convergence criterion
#' # LCDM
#' lcdm <- GDINA(dat = dat, Q = Q, model = "UDF", design.matrix = D,
#' linkfunc = "logit", control=list(conv.type="neg2LL"),solver="slsqp")
#'
#' # identity link GDINA
#' iGDINA <- GDINA(dat = dat, Q = Q, model = "GDINA",
#' control=list(conv.type="neg2LL"),solver="slsqp")
#'
#' # compare two models => identical
#' anova(lcdm,iGDINA)
#'
#'####################################
#'#           Example 13b.           #
#'#     user specified design matrix #
#'#            RRUM                  #
#'####################################
#'
#' dat <- sim30GDINA$simdat
#' Q <- sim30GDINA$simQ
#'
#' # specify design matrix for each item => must be a list
#' # D can be defined by the user
#' D <- lapply(rowSums(Q),designmatrix,model="ACDM")
#' # for comparison, use change in -2LL as convergence criterion
#' # RRUM
#' logACDM <- GDINA(dat = dat, Q = Q, model = "UDF", design.matrix = D,
#' linkfunc = "log", control=list(conv.type="neg2LL"),solver="slsqp")
#'
#' # identity link GDINA
#' RRUM <- GDINA(dat = dat, Q = Q, model = "RRUM",
#'               control=list(conv.type="neg2LL"),solver="slsqp")
#'
#' # compare two models => identical
#' anova(logACDM,RRUM)
#'
#'####################################
#'#           Example 14.            #
#'#     Multiple-strategy DINA model #
#'####################################
#'
#' Q <- matrix(c(1,1,1,1,0,
#' 1,2,0,1,1,
#' 2,1,1,0,0,
#' 3,1,0,1,0,
#' 4,1,0,0,1,
#' 5,1,1,0,0,
#' 5,2,0,0,1),ncol = 5,byrow = TRUE)
#' d <- list(
#'   item1=c(0.2,0.7),
#'   item2=c(0.1,0.6),
#'   item3=c(0.2,0.6),
#'   item4=c(0.2,0.7),
#'   item5=c(0.1,0.8))
#'
#'   set.seed(12345)
#'sim <- simGDINA(N=1000,Q = Q, delta.parm = d,
#'                model = c("MSDINA","MSDINA","DINA",
#'                          "DINA","DINA","MSDINA","MSDINA"))
#'
#'# simulated data
#'dat <- extract(sim,what = "dat")
#' # estimation
#' # MSDINA need to be specified for each strategy
#' est <- GDINA(dat,Q,model = c("MSDINA","MSDINA","DINA",
#'                              "DINA","DINA","MSDINA","MSDINA"))
#' coef(est,"delta")
#'}
#'
GDINA <-
  function(dat, Q, model = "GDINA", sequential = FALSE, att.dist = "saturated", mono.constraint = FALSE,
           group = NULL,linkfunc = NULL, design.matrix = NULL,
           latent.var = "att", att.prior = NULL, att.str = FALSE, verbose = 1,
           higher.order = list(),loglinear = 2,catprob.parm = NULL,control=list(),
           item.names = NULL, solver = NULL,
           nloptr.args=list(),auglag.args=list(),solnp.args = list(),...)
  {

    s1 <- Sys.time()
    GDINAcall <- match.call()
    if (exists(".Random.seed", .GlobalEnv)){
      oldseed <- .GlobalEnv$.Random.seed
      on.exit(.GlobalEnv$.Random.seed <- oldseed)
    }else{
      on.exit(rm(".Random.seed", envir = .GlobalEnv))
    }
    if(missing(dat)) missingMsg(dat)
    if(missing(Q)) missingMsg(Q)

    saturated <- dots("saturated",list(type=1,prior=0),...)
    item.prior <- dots("item.prior",list(on = FALSE, type="auto",beta=c(1.001,1.001),normal=c(0,5)),...)
    constr.pairs <- dots("constr.pairs",NULL,...)


    ret <- Est(dat = dat, Q = Q, model = model, sequential = sequential,
             att.dist = att.dist, att.prior = att.prior,saturated=saturated,
             att.str = att.str, mono.constraint=mono.constraint,
             group = group, latent.var=latent.var, verbose = verbose,
             catprob.parm = catprob.parm,loglinear = loglinear,
             item.names = item.names, control = control, item.prior = item.prior,
             nloptr_args = nloptr.args,auglag_args=auglag.args,solnp_args = solnp.args,
             linkfunc = linkfunc,higher.order = higher.order, solver = solver,
             DesignMatrices = design.matrix,ConstrPairs = constr.pairs)


    s2 <- Sys.time()

    ret$extra <- list(timeused = s2 - s1, start.time = s1, end.time = s2, call = GDINAcall)
    class(ret) <- "GDINA"
    invisible(ret)
  }
