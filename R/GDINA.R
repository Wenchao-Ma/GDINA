#' Calibrate dichotomous and polytomous responses
#'
#' \code{GDINA} calibrates the generalized deterministic inputs, noisy and
#' gate (G-DINA; de la Torre, 2011) model for dichotomous responses, and its extension, the sequential
#' G-DINA model (Ma, & de la Torre, 2016a) for ordinal and nominal responses.
#' By setting appropriate constraints, the deterministic inputs,
#' noisy and gate (DINA; de la Torre, 2009; Junker & Sijtsma, 2001) model,
#' the deterministic inputs, noisy or gate (DINO; Templin & Henson, 2006)
#' model, the reduced reparametrized unified model (R-RUM; Hartz, 2002),
#' the additive CDM (A-CDM; de la Torre, 2011), and the linear logistic
#' model (LLM; Maris, 1999) can also be calibrated. Note that the LLM is equivalent to
#' the C-RUM (Hartz, 2002), a special case of the GDM (von Davier, 2008), and that the R-RUM
#' is also known as a special case of the generalized NIDA model (de la Torre, 2011).
#' Different models can be fitted to different
#' items in a single test. The attributes can be either dichotomous or polytomous
#' (Chen & de la Torre, 2013). Joint attribute distribution can be saturated, structured or higher-order
#' (de la Torre & Douglas, 2004) when attributes are binary.
#' Marginal maximum likelihood method with Expectation-Maximization (MMLE/EM) alogrithm
#' is used for item parameter estimation.
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
#' \deqn{P(\bm{\alpha}_{lj}^*)=\delta_{j0}+\sum_{k=1}^{K_j^*}\delta_{jk}\alpha_{lk}+
#' \sum_{k'=k+1}^{K_j^*}\sum_{k=1}^{K_j^*-1}\delta_{jkk'}\alpha_{lk}\alpha_{lk'}+\cdots+
#' \delta_{j12{\cdots}K_j^*}\prod_{k=1}^{K_j^*}\alpha_{lk},
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
#'
#'@section Model Estimation:
#'
#' The MMLE/EM algorithm is implemented in this package. For G-DINA, DINA and DINO models, closed-form solutions can be found.
#' Specifically, for the G-DINA model, \deqn{P(\alpha_{lj}^*)=R_{jl}/N_{jl}} where \eqn{R_{jl}} is the expected number of examinees with attribute pattern \eqn{\alpha_{lj}^*}
#' answering item \eqn{j} correctly and \eqn{N_{jl}} is the expected number of examinees with attribute pattern \eqn{\alpha_{lj}^*}.
#' For DINA or DINO model, \eqn{R_{jl}} and \eqn{N_{jl}} are collapsed for latent classes having the same probability of success.
#' See de la Torre (2009) and de la Torre (2011) for details.
#'
#' For ACDM, LLM and RRUM, closed-form solutions do not exist, and therefore some general optimization techniques are
#' adopted in M-step. See Ma, Iaconangelo and de la Torre (2016) for details.
#' The selection of optimization techniques mainly depends on whether
#' some specific constraints need to be added. It should
#' be noted that adding monotone constraints to the G-DINA model may dramatically increase running time especially when the number of required
#' attributes are large.
#'
#' The sequential G-DINA model can be estimated as in Ma & de la Torre (2016a) using optimization techniques. However,
#' Ma & de la Torre (2016b) found that the sequential G-DINA, DINA and DINO models can be estimated using
#' close-form solutions, which can be implemented in a straightforward
#' manner using the observation-coding (Tutz, 1997).
#'
#' For estimating the joint attribute
#' distribution, by default, an empirical Bayes method (\code{saturated}; Carlin & Louis, 2000) is adopted, which is referred to as
#' the saturated attribute structure. Specifically,
#' the prior distribution of joint attributes is uniform at the beginning, and then updated after
#' each EM iteration based on the posterior distribution.
#'
#' The joint attribute distribution can also be modeled using some higher-order IRT models, which is referred to as
#' higher-order attribute structure. The higher-order attribute structure was originally proposed by de la Torre
#' and Douglas (2004) for the DINA model. It has been extended in this package for the G-DINA model, DINA, DINO, A-CDM, LLM and RRUM.
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
#' responses of \eqn{N} individuals to \eqn{J} items. Missing values need to be coded as \code{NA}.
#' @param Q A required \eqn{J \times K} item or category and attribute association matrix, wher \eqn{J} represents the number of
#'    items or nonzero categories and \eqn{K} represents the number of attributes. For binary attributes,
#'    entry 1 indicates that the attribute is measured by the item, and 0 otherwise.
#'    For polytomous attributes, non-zero elements indicate the level
#'    of attributes that are needed for an individual to answer the item correctly (see Chen, & de la Torre, 2013).
#'    Note that for polytomous items, the sequential G-DINA
#'    model is used and either restricted or unrestricted category-level Q-matrix is needed.
#'    In the category-level Q-matrix, the first column gives the item number, which must be numeric and match the number of column in the data.
#'    The second column indicates the category number. See \code{Examples}.
#' @param model A vector for each item or nonzero category, or a scalar which will be used for all
#'    items or nonzero categories to specify the CDMs fitted. The possible options
#'    include \code{"GDINA"},\code{"DINA"},\code{"DINO"},\code{"ACDM"},\code{"LLM"}, and \code{"RRUM"}.
#'    It is also possible to specify CDMs using numbers. Particularly, 0,1,2,3,4 and 5 represents
#'    \code{"GDINA"},\code{"DINA"},\code{"DINO"},\code{"ACDM"},\code{"LLM"}, and \code{"RRUM"}, respectively.
#' @param sequential logical; \code{TRUE} if the sequential model is fitted for polytomous responses.
#' @param group a scalar indicating which column in \code{dat} is group indicator or
#'    a numerical vector indicating the group each individual belongs to. If it is a vector,
#'    its length must be equal to the number of individuals. Only at most two groups can be handled currently.
#' @param att.dist How is the joint attribute distribution estimated? It can be \code{saturated}, indicating that
#'    the proportion parameter for each permissible latent class is estimated separately; \code{higher.order}, indicating
#'    that a higher-order joint attribute distribution is assumed (higher-order model can be specified in \code{higher.order} argument);
#'    or \code{fixed}, indicating that the weights specified in \code{att.prior} argument are fixed in the estimation process.
#'    If \code{att.prior} is not specified, a uniform joint attribute distribution is employed initially.
#'    If different groups have different joint attribute distributions, specify \code{att.dist} as a character vector with the same
#'    number of elements as the number of groups.
#' @param item.names A vector giving the item names. By default, items are named as "Item 1", "Item 2", etc.
#' @param higher.order A list specifying the higher-order joint attribute distribution with the following components:
#'    (1) \code{model} - a character indicating the IRT model for higher-order joint attribute distribution. Can be
#'    \code{"2PL"}, \code{"1PL"} or \code{"Rasch"}, representing two parameter logistic IRT model,
#'    one parameter logistic IRT model and Rasch model,
#'    respectively. For \code{"1PL"} model, a common slope parameter is
#'    estimated (see \code{Details}). \code{"Rasch"} is the default model when \code{att.dist = "higher.order"}.
#'    (2) \code{method} - a character indicating the algorithm for the higher-order structural parameter estimation;
#'    Can be either \code{"BL"}, \code{"MMLE"}, which is the default, or \code{"BMLE"}, which allows parameter priors to be imposed.
#'    (3) \code{nquad} - a scalar specifying the number of integral nodes. (4) \code{type} - a character specifying
#'    whether all higher-order structural parameters are estimated at the same time (i.e., \code{type="testwise"}) or
#'    estimated attribute by attribute (i.e., \code{type="attwise"}, only applicable when \code{method="MMLE"} or \code{method="BMLE"}).
#'    (5) \code{slope.range} - a vector of length two specifying
#'    the range of slope parameters. (6) \code{intercept.range} - a vector of length two specifying
#'    the range of intercept parameters. (7) \code{slope.prior} - a vector of length two specifying
#'    the mean and variance of log(slope) parameters, which are assumed normally distributed. (8) \code{intercept.prior} -
#'    a vector of length two specifying the mean and variance of intercept parameters, which are assumed normally distributed.
#' @param mono.constraint logical; \code{TRUE} indicates that \eqn{P(\bm{\alpha}_1) <=P(\bm{\alpha}_2)} if
#'    for all \eqn{k}, \eqn{\alpha_{1k} < \alpha_{2k}}. Can be a vector for each item or nonzero category or a scalar which will be used for all
#'    items to specify whether monotonicity constraint should be added.
#' @param catprob.parm A list of initial success probability parameters for each nonzero category.
#' @param verbose How to print calibration information
#'     after each EM iteration? Can be 0, 1 or 2, indicating to print no information,
#'     information for current iteration, or information for all iterations.
#' @param att.prior A vector of length \eqn{2^K} for single group estimation, or a matrix of dimension \eqn{2^K\times} no. of groups to specify
#'    attribute prior distribution for \eqn{2^K} latent classes for all groups. Only applicable for dichotomous attributes.
#'    The sum of all elements does not have to be equal to 1; however, it will be transformed so that the sum is equal to 1
#'    before model calibration.
#'    The label for each latent class can be obtained by calling \code{attributepattern(K)}. See \code{examples} for more info.
#' @param att.str logical; are attributes structured?
#' @param nstarts how many sets of starting values? The default is 1.
#' @param conv.crit The convergence criterion for max absolute change in item parameters or deviance.
#' @param conv.type How is the convergence criterion evaluated? Can be \code{"max.p.change"}, indicating
#'    the maximum absolute change in success probabilities, or \code{"dev.change"}, representing
#'    the absolute change in deviance.
#' @param maxitr A vector for each item or nonzero category, or a scalar which will be used for all
#'    items or nonzero categories to specify the maximum number of EM cycles allowed.
#' @param lower.p A vector for each item or nonzero category, or a scalar which will be used for all
#'    items or nonzero categories to specify the lower bound for success probabilities. The default is \code{1e-4} for all items.
#' @param upper.p A vector for each item or nonzero category, or a scalar which will be used for all
#'    items or nonzero categories to specify the upper bound for success probabilities. The default is 0.9999 for all items.
#' @param lower.prior The lower bound for prior weights. Only applicable for nonstructured attributes.
#'    The default value is -1, which means the lower bound is \eqn{1/2^K/100}.
#' @param randomseed Random seed for generating initial item parameters. The default random seed is 123456.
#' @param smallNcorrection A numeric vector with two elements specifying the corrections applied when the expected number of
#' individuals in some latent groups are too small. If the expected no. of examinees is less than the second element,
#' the first element and two times the first element will be added to the numerator and denominator of the closed-form solution of
#' probabilities of success. Only applicable for the G-DINA, DINA and DINO model estimation without monotonic constraints.
#' @param Mstep.warning Logical; Should the warning message in Mstep, if any, be output immediately.
#' @param diagnosis Run in diagnostic mode? If it is 1 or 2, some intermediate results obtained in each iteration can be extracted.
#' @param optimizer A string indicating which optimizer should be used in M-step.
#' @param optim.control Control options for optimizers in the M-step. Only available when \code{optimizer} is one specific optimization
#' method, including \code{BFGS} from \link[stats]{optim}, \link[nloptr]{slsqp}, \link[Rsolnp]{solnp} and \link[alabama]{auglag}.
#' For the \link[alabama]{auglag} method, \code{optim.control} specifies \code{control.outer}.
#'
#' @author {Wenchao Ma, Rutgers University, \email{wenchao.ma@@rutgers.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @seealso See \code{\link{autoGDINA}} for Q-matrix validation, item level model comparison and model calibration
#' in one run; See \code{\link{itemfit}} for item fit analysis, \code{\link{Qval}} for Q-matrix validation,
#' \code{\link{modelcomp}} for item level model comparison and \code{\link{simGDINA}} for data simulation.
#' Also see \code{gdina} in \pkg{CDM} package for the G-DINA model estimation.
#'
#' @return \code{GDINA} returns an object of class \code{GDINA}. Methods for \code{GDINA} objects
#'  include \code{\link{extract}} for extracting various components, \code{\link{itemparm}}
#'  for extracting item parameters, \code{\link{personparm}}
#'  for calculating person parameters, \code{summary} for summary information.
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
#' Carlin, B. P., & Louis, T. A. (2000). Bayes and empirical bayes methods for data analysis. New York, NY: Chapman & Hall
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
#' Junker, B. W., & Sijtsma, K. (2001). Cognitive assessment models with few assumptions, and connections with nonparametric
#' item response theory. \emph{Applied Psychological Measurement, 25}, 258-272.
#'
#' Ma, W., & de la Torre, J. (2016a). A sequential cognitive diagnosis model for polytomous responses. \emph{British Journal of Mathematical and Statistical Psychology. 69,} 253-275.
#'
#' Ma, W., & de la Torre, J. (2016b, July). A Q-matrix validation method for the sequential G-DINA model. Paper presented at the 80th International Meeting of the Psychometric Society, Asheville, NC.
#'
#' Ma, W., Iaconangelo, C., & de la Torre, J. (2016). Model similarity, model selection and attribute classification.
#' \emph{Applied Psychological Measurement, 40}, 200-217.
#'
#' Maris, E. (1999). Estimating multiple classification latent class models. \emph{Psychometrika, 64}, 187-212.
#'
#' Tatsuoka, K. K. (1983). Rule space: An approach for dealing with misconceptions based on
#' item response theory. \emph{Journal of Educational Measurement, 20}, 345-354.
#'
#' Templin, J. L., & Henson, R. A. (2006). Measurement of psychological disorders using cognitive diagnosis models.
#' \emph{Psychological Methods, 11}, 287-305.
#'
#' Tutz, G. (1997). Sequential models for ordered responses. In W.J. van der Linden & R. K. Hambleton (Eds.), Handbook of modern item response theory p. 139-152). New York, NY: Springer.
#'
#'
#' @export
#' @importFrom  nloptr slsqp
#' @useDynLib GDINA
#' @importFrom Rcpp sourceCpp
#' @import numDeriv
#' @importFrom alabama auglag
#' @importFrom graphics plot axis points text abline
#' @import ggplot2
#' @import stats
#' @importFrom utils combn
#'
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
#' # item parameters
#' # see ?itemparm
#' itemparm(mod1) # item probabilities of success for each latent group
#' itemparm(mod1, withSE = TRUE) # item probabilities of success & standard errors
#' itemparm(mod1, what = "delta") # delta parameters
#' itemparm(mod1, what = "delta",withSE=TRUE) # delta parameters
#' itemparm(mod1, what = "gs") # guessing and slip parameters
#' itemparm(mod1, what = "gs",withSE = TRUE) # guessing and slip parameters & standard errors
#'
#' # person parameters
#' # see ?personparm
#' personparm(mod1) # EAP estimates of attribute profiles
#' personparm(mod1, what = "MAP") # MAP estimates of attribute profiles
#' personparm(mod1, what = "MLE") # MLE estimates of attribute profiles
#'
#' #plot item response functions for item 10
#' plotIRF(mod1,item = 10)
#' plotIRF(mod1,item = 10,errorbar = TRUE) # with error bars
#' plotIRF(mod1,item = c(6,10))
#'
#' # Use extract function to extract more components
#' # See ?extract
#'
#' # ------- DINA model --------#
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' mod2 <- GDINA(dat = dat, Q = Q, model = "DINA")
#' mod2
#' itemparm(mod2, what = "gs") # guess and slip parameters
#' itemparm(mod2, what = "gs",withSE = TRUE) # guess and slip parameters and standard errors
#'
#' # Model comparison at the test level via likelihood ratio test
#' anova(mod1,mod2)
#'
#' # -------- DINO model -------#
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' mod3 <- GDINA(dat = dat, Q = Q, model = "DINO")
#' #slip and guessing
#' itemparm(mod3, what = "gs") # guess and slip parameters
#' itemparm(mod3, what = "gs",withSE = TRUE) # guess and slip parameters + standard errors
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
#' itemparm(mod11d,"delta")
#' itemparm(mod11d,"rrum")
#'
#'####################################
#'#           Example 3.             #
#'#        Model estimations         #
#'# With Higher-order att structure  #
#'####################################
#'
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' # --- Higher order G-DINA model ---#
#' mod12 <- GDINA(dat = dat, Q = Q, model = "DINA",
#'                att.dist="higher.order",higher.order=list(method="MMLE",nquad=31))
#' hoest=hoparm(mod12) # extract higher-order parameters
#' hoest$theta # ability
#' hoest$lambda # structural parameters
#' # --- Higher order DINA model ---#
#' mod22 <- GDINA(dat = dat, Q = Q, model = "DINA",
#'                att.dist="higher.order",higher.order=list(method="BMLE"))
#' # --- Higher order DINO model ---#
#' mod23 <- GDINA(dat = dat, Q = Q, model = "DINO",att.dist="higher.order")
#' # --- Higher order ACDM model ---#
#' mod24 <- GDINA(dat = dat, Q = Q, model = "ACDM",
#'                att.dist="higher.order",higher.order=list(model="1PL"))
#' # --- Higher order LLM model ---#
#' mod25 <- GDINA(dat = dat, Q = Q, model = "LLM",att.dist="higher.order")
#' # --- Higher order RRUM model ---#
#' mod26 <- GDINA(dat = dat, Q = Q, model = "RRUM",att.dist="higher.order")
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
#' # divergent structure A1->A2->A3;A1->A4->A5;A1->A4->A6
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
#'  simD <- simGDINA(N,Q,gs.parm = gs,
#'                    model = "DINA",attribute = true.att)
#'  dat <- extract(simD,"dat")
#'
#' modp1 <- GDINA(dat = dat, Q = Q, att.prior = struc$att.prob, att.str = TRUE, att.dist = "saturated")
#' modp1
#' # Note that fixed priors were used for all iterations
#' extract(modp1,what = "att.prior")
#' # Posterior weights were slightly different
#' extract(modp1,what = "posterior.prob")
#' modp2 <- GDINA(dat = dat, Q = Q, att.prior = struc$att.prob, att.str = TRUE, att.dist = "fixed")
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
#'  mod.initial <- GDINA(dat,Q,maxit=20) # estimation- only 10 iterations for illustration purposes
#'  par <- itemparm(mod.initial,digits=8)
#'  weights <- extract(mod.initial,"posterior.prob",digits=8) #posterior weights
#'  # use the weights as the priors
#'  mod.fix <- GDINA(dat,Q,catprob.parm = par,att.prior=c(weights),maxitr=0) # re-estimation
#'  anova(mod.initial,mod.fix) # very similar - good approximation most of time
#'  # prior used for the likelihood calculation for the last step
#'  priors <- extract(mod.initial,"att.prior")
#'  # use the priors as the priors
#'  mod.fix2 <- GDINA(dat,Q,catprob.parm = par,att.prior=priors,maxitr=0) # re-estimation
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
#' itemparm(sDINA) # processing function
#' itemparm(sDINA,"itemprob") # success probabilities for each item
#' itemparm(sDINA,"LCprob") # success probabilities for each category for all latent classes
#'
#'####################################
#'#           Example 10.            #
#'#    Multiple-Group G-DINA model   #
#'####################################
#' Q <- sim10GDINA$simQ
#'
#' # parameter simulation
#' # Group 1 - female
#' N1 <- 2000
#' gs1 <- matrix(rep(0.1,2*nrow(Q)),ncol=2)
#' # Group 2 - male
#' N2 <- 2000
#' gs2 <- matrix(rep(0.2,2*nrow(Q)),ncol=2)
#'
#' # data simulation for each group
#' sim1 <- simGDINA(N1,Q,gs.parm = gs1,model = "DINA")
#' sim2 <- simGDINA(N2,Q,gs.parm = gs2,model = "DINO")
#'
#' # combine data
#' # see ?bdiagMatrix
#' dat <- bdiagMatrix(list(extract(sim1,"dat"),extract(sim2,"dat")),fill=NA)
#' Q <- rbind(Q,Q)
#' gr <- rep(c("female","male"),each=2000)
#' # Fit G-DINA model
#' mg.est <- GDINA(dat = dat,Q = Q,group = gr)
#' summary(mg.est)
#' extract(mg.est,"posterior.prob")
#'
#' # Fit G-DINA model with different joint attribute dist.
#' mg.est2 <- GDINA(dat = dat,Q = Q,group = gr,
#' att.dist = c("saturated","fixed"))
#' summary(mg.est2)
#'
#' }
#'
GDINA <-
  function(dat, Q, model = "GDINA", sequential = FALSE, att.dist = "saturated", att.prior = NULL,
           att.str = FALSE, mono.constraint = FALSE, group = NULL,
           verbose = 1, catprob.parm = NULL,
           lower.p = 0.0001,upper.p = 0.9999, item.names = NULL,
           nstarts = 1, conv.crit = 0.001, lower.prior = -1,
           conv.type = "max.p.change", maxitr = 1000,
           digits = 4,diagnosis = 0, Mstep.warning = FALSE,optimizer = "all",
           randomseed = 123456, smallNcorrection = c(0.0005,0.001),
           higher.order = list(model="Rasch",method="MMLE",nquad=19,type = "testwise",
                               slope.range=c(0.1,5),intercept.range=c(-3,3),
                               slope.prior=c(0,0.25),intercept.prior=c(0,1)),
           optim.control=list())
  {
    s1 <- Sys.time()
    GDINAcall <- match.call()

    if (exists(".Random.seed", .GlobalEnv)) oldseed <- .GlobalEnv$.Random.seed else  oldseed <- NULL

    if (is.null(group)){
      no.mg <- 1
      gr <- rep(1L,nrow(dat))
      gr.label <- "all data"
    }else {
      # multiple groups
      # scalar -> group variable column
      if (length(group) == 1L) {
        if(!is.positiveInteger(group)) stop("group specification is not correct.",call. = FALSE)
        gr <- dat[, group] # group indicator variable - numeric
        dat <- dat[, -group] # responses
      } else{
        if (nrow(dat) != length(group)) stop("The length of group variable must be equal to the number of individuals.", call. = FALSE)
        gr <- group # group indicator variable
      }

      gr.label <- unique(gr) # group labels
      no.mg <- length(gr.label) # the number of groups
      if (!is.numeric(gr) || max(gr) > no.mg) {
        for (g in 1L:no.mg) gr[gr == gr.label[g]] <- g
        gr <- as.numeric(gr) # numeric variable
      }
if(no.mg>1&&length(att.dist)==1) att.dist <- rep(att.dist,no.mg)
    }
    lambda <- NULL
    lam <- list() # lambda for multiple groups
if(any(att.dist=="higher.order")) {

  higher.order$model <- ifelse(is.null(higher.order$model),"Rasch",higher.order$model)
  higher.order$method <- ifelse(is.null(higher.order$method),"MMLE",higher.order$method)
  higher.order$nquad <- ifelse(is.null(higher.order$nquad),19,higher.order$nquad)
  higher.order$type <- ifelse(is.null(higher.order$type),"testwise",higher.order$type)
  if(is.null(higher.order$slope.range)) higher.order$slope.range <- c(0.1,5)
  if(is.null(higher.order$intercept.range)) higher.order$intercept.range <- c(-3,3)
  if(is.null(higher.order$slope.prior)) higher.order$slope.prior <- c(0,0.25)
  if(is.null(higher.order$intercept.prior)) higher.order$intercept.prior <- c(0,1)
  if(higher.order$method=="MMLE") higher.order$slope.prior <- higher.order$intercept.prior <- c(NA,NA)
  if(higher.order$type=="attwise"&&higher.order$model=="1PL") warning("Estimating 1PL higher-order model using type='attwise' manner is not recommended.",call. = FALSE)

}

     inputcheck(dat = dat, Q = Q, model = model, sequential = sequential, att.dist = att.dist,
               verbose = verbose, catprob.parm = catprob.parm, mono.constraint = mono.constraint,
               att.prior = att.prior,
               lower.p = lower.p, upper.p = upper.p, att.str = att.str, nstarts = nstarts,
               conv.crit = conv.crit, maxitr = maxitr,
               digits = digits, diagnosis = diagnosis)
     if(any(is.na(dat))){
       # some missings individuals with one or fewer valid response are
       # removed
       del.ind <- which(rowSums(1L - is.na(dat)) < 2L, arr.ind = TRUE)
       if (length(del.ind) > 0L) {
         warning(length(del.ind), " individuals with one or fewer valid responses are removed.")
         dat <- dat[-del.ind, ]
         gr <- gr[-del.ind]
       }
     }
    # polytomous responses -> dichotomous responses
    # copy Q and data
    Qc <- Q
    dat_d <- dat
    if (sequential) {
      dat <- seq_coding(dat, Q)
      if (is.null(item.names)) item.names <- paste("Item", Qc[, 1], "Cat", Qc[, 2])

      Q <- Q[, -c(1, 2)]
    } else {
      if (max(dat, na.rm = TRUE) > 1) stop("Maximum response is greater than 1 - set sequential = TRUE to fit a sequential model.", call. = FALSE)
      if (is.null(item.names)) {
        if(is.null(colnames(dat))) item.names <- paste("Item", 1:nrow(Q)) else item.names <- colnames(dat)
      }

    }


    N <- nrow(dat)

    J <- ncol(dat)

    K <- ncol(Q)

    Kj <- rowSums(Q>0)  # vector with length of J

    Lj <- 2^Kj

    L <- no_LC(Q)  # The number of latent groups

    M <- c("GDINA", "DINA", "DINO", "ACDM", "LLM", "RRUM")

    model <- model.transform(model, J)
    model.names <- M[model + 1]
    names(model.names) <- item.names
    if (length(lower.p) == 1) {
      lower.p <- rep(lower.p, J)
    } else {
      if (length(lower.p) != J)
        stop("lower.p must have length of 1 or J.", call. = FALSE)
    }
    if (length(upper.p) == 1) {
      upper.p <- rep(upper.p, J)
    } else {
      if (length(upper.p) != J)
        stop("upper.p must have length of 1 or J.", call. = FALSE)
    }

    if(length(maxitr)==1) {
      vmaxitr <- rep(maxitr,J)
    }else if(length(maxitr)!=J){
      warning("Length of maxitr must be equal to 1 or the number of nonzero categories.",call. = FALSE)
    }else{
      vmaxitr <- maxitr
      maxitr <- max(maxitr)
    }
    # location of item parameters
    parloc <- eta.loc(Q)  #J x L

    # Generate the log prior - a L x no.mg matrix
    if (is.null(att.prior)) {
      att.prior <- matrix(rep(1/L, L),ncol = no.mg)
      logprior <- matrix(log(att.prior),nrow = L,ncol = no.mg)
    } else {
      if(is.vector(att.prior)) {
        att.prior <- matrix(att.prior,ncol = no.mg) # vector -> matrix
      }
      if(is.matrix(att.prior)){
        if (nrow(att.prior) != L||ncol(att.prior)!=no.mg)
          stop("Joint attribute distribution priors must be a vector of length 2^K or a matrix with 2^K rows if specified.", call. = FALSE)
      }

      if (any(att.prior < 0)) stop("Joint attribute distribution prior can only contain numerical values between 0 and 1.", call. = FALSE)
        logprior <- log(att.prior/matrix(colSums(att.prior),nrow = L,ncol = no.mg,byrow = TRUE))
    }

    # generating initial item parameters
    if (is.null(catprob.parm)) {
      item.parm <- initials(Q, nstarts, randomseed = randomseed)  #a list with nstarts matrices of size J x 2^Kjmax
      ###### Multiple starting values
      if (nstarts > 1L) {
        neg2LL <- NULL
        for (i in 1L:nstarts) {
          LC.Prob <- uP(as.matrix(parloc), as.matrix(item.parm[[i]]))
          # calculate likelihood and posterior
          neg2LL <- c(neg2LL, -2 * Lik(mP = LC.Prob,
                                       mX = as.matrix(dat),
                                       vlogPrior = as.matrix(logprior),
                                       vgroup = gr)$LL)
        }
        item.parm <- item.parm[[which.min(neg2LL)]]

      }
    } else {
      item.parm <- l2m(catprob.parm)
    }

    initial.parm <- item.parm

    dif <- 10L
    itr <- 0L
    diag.itemprob <- diag.RN <- diag.likepost <- diag.lambda <- diag.opts <- vector("list", maxitr)
    diag.itrmaxchange <- NULL
    deltas <- calc_delta(item.parm, model, Kj)
    # print(deltas)

    # cat('\nIteration maxchange -2LL\n') E-M algorithm
    dif.p <- dif.LL <- LL.1 <- LL.2 <- 0

    ############### START of while loop ##########
    designMlist <- vector("list",J)
    for(j in 1:J) designMlist[[j]] <- designmatrix(Kj[j],model[j])


    while (dif > conv.crit && itr < maxitr)
    {
      # length of J indicating whether a nonzero cat should  be est. (TRUE) or not (FALSE)
      est.bin <- vmaxitr>itr
      item.parm_copy <- item.parm
      # --------E step probability of getting each score for each latent
      # class
      # print(item.parm)
      LC.Prob <- uP(as.matrix(parloc), as.matrix(item.parm))
      LC.Prob[is.nan(LC.Prob)] <- 0
      # calculate likelihood and posterior
      # vlogPrior: col vector for single group; matrix with g columns for g groups;
      # g: group indicator from 1 to g;
      # print(logprior)
      likepost <- Lik(mP = LC.Prob,
                      mX = as.matrix(dat),
                      vlogPrior = as.matrix(logprior),
                      vgroup = gr)

      # expected # of examinees of answering items correctly in each latent
      # class expected # of examinees for each latent class
      RN <- NgRg(mlogPost = likepost$logpost,
                 mX = as.matrix(dat),
                 mloc = as.matrix(parloc))
      # print(RN)
      LL.1 <- -2*likepost$LL
      dif.LL <- abs(LL.1-LL.2)
      if(tolower(conv.type)=="dev.change") {
        dif <- dif.LL
        if (dif.LL<=conv.crit) break
      }
      LL.2 <- LL.1



      if (att.str) smallNcorrection <- c(-1,-1)
      optims <- Mstep(Kj = Kj, RN = RN, model = model, itmpar = item.parm, delta = deltas,
                      constr = mono.constraint, correction = smallNcorrection, lower.p = lower.p,
                      upper.p = upper.p, itr = itr + 1, warning.immediate = Mstep.warning,
                      optimizer = optimizer, designmatrices = designMlist,optim.control = optim.control,est.bin = est.bin)
      item.parm <- optims$item.parm
      deltas <- optims$delta

      dif.p <- max(abs(item.parm - item.parm_copy),na.rm = TRUE)
      if(tolower(conv.type)=="max.p.change") dif <- dif.p

      if (diagnosis>=1) {
        diag.itemprob[[itr+1]] <- item.parm_copy
        diag.itrmaxchange <- rbind(diag.itrmaxchange,c(dif,-2*likepost$LL))
        if (diagnosis>=2) {
          diag.likepost[[itr+1]] <- likepost
          diag.RN[[itr+1]] <- RN
          diag.opts[[itr+1]] <- optims$optims
        }

      }

      itr <- itr + 1
      if(verbose==1) {
        cat('\rIteration =',itr,' Max change =',formatC(dif,digits = 4, format = "f"),
            ' Deviance =',formatC(-2*likepost$LL,digits = 4, format = "f"),'                                                                 ')
      }else if (verbose==2) {
        cat('Iteration =',itr,' Max change =',formatC(dif,digits = 4, format = "f"),
            ' Deviance =',formatC(-2*likepost$LL,digits = 4, format = "f"),"                                                                \n")
      }
      # update prior for each group

      for(g in 1:no.mg){
        if (att.dist[g]=="saturated"){
          if(att.str) {
            logprior[,g] <- likepost$logprior[,g]
          }else{
            lower.prior <- ifelse(lower.prior==-1,1/L/1000,lower.prior)
            priori <- exp(likepost$logprior[,g])
            priori[which(priori<lower.prior)] <- lower.prior
            priori <- priori/sum(priori)
            logprior[,g] <- log(priori)
          }
        }else if (att.dist[g]=="higher.order")
        {

          Rl <- c(colSums(exp(likepost$logpost[which(gr==g),])))
          a <- b <- NULL
          if (!is.null(lambda)) a <- lambda[,1];b <- lambda[,2]
          HO.out <- HO.est(Rl = Rl, K = K, N = N,
                           IRTmodel = higher.order$model,nnodes = higher.order$nquad, type = higher.order$type,
                           a.bound = higher.order$slope.range,b.bound = higher.order$intercept.range,
                           loga.prior = higher.order$slope.prior,b.prior = higher.order$intercept.prior,
                           method = higher.order$method, a = a, b = b)
          logprior[,g] <- HOlogprior <- HO.out$logprior
          lambda <- HO.out$lambda
          colnames(lambda) <- c("slope", "intercept")
          rownames(lambda) <- paste("Attribute", 1:K, sep = " ")
          lam[[g]] <- lambda
          # print(HO.out$lambda)
          if (diagnosis >= 1) {
            if(no.mg==1) diag.lambda[[itr+1]] <- lambda else diag.lambda[[itr+1]] <- lam
          }
        }else if (att.dist[g]=="fixed")
        {
          logprior[,g] <- log(att.prior[,g])
        }
      }

}
    ################ END of while loop ##########
    if (diagnosis >= 1 & itr < maxitr) {
      diag.itemprob[(itr + 1):maxitr] <- diag.RN[(itr + 1):maxitr] <-
        diag.likepost[(itr + 1):maxitr] <- diag.lambda[(itr + 1):maxitr] <- diag.opts[(itr + 1):maxitr] <- NULL
      names(diag.itemprob) <- names(diag.likepost) <- names(diag.lambda) <-
        names(diag.likepost) <- names(diag.opts) <- names(diag.RN) <- c("Initial",  paste("Iter", 1:(itr - 1)))
      rownames(diag.itrmaxchange) <- paste("Iter", 1:itr)
      colnames(diag.itrmaxchange) <- c("Max change", "Deviance")

    }


    # ---------update posterior probability of getting each score for each
    # latent class
    LC.Prob <- uP(as.matrix(parloc), as.matrix(item.parm))
    LC.Prob[is.nan(LC.Prob)] <- 0

    likepost <- Lik(mP = LC.Prob,
                    mX = as.matrix(dat),
                    vlogPrior = as.matrix(logprior),
                    vgroup = gr)
    # print(likepost)
    logprior0 <- logprior # old log priors
    logprior <- likepost$logprior # marginal weights
    RN <- NgRg(mlogPost = likepost$logpost,
               mX = as.matrix(dat),
               mloc = as.matrix(parloc))
    # -----------------Test Fit information----------------#

    neg2LL <- -2 * likepost$LL
# print(neg2LL)
    npar <- 0
    for (j in 1:J) {
      if (model[j] == 0) {
        # GDINA
        npar <- npar + Lj[j]
      } else if (model[j] == 1) {
        # DINA
        npar <- npar + 2
      } else if (model[j] == 2) {
        # DINO
        npar <- npar + 2
      } else {
        # ACDM,LLM or RRUM
        npar <- npar + Kj[j] + 1
      }
    }
    item.npar <- npar  #item parameters
    for(g in 1:no.mg){
      if (att.dist[g]=="saturated") {
        if (!att.str) {
          npar <- npar + L - 1
        } else {
          npar <- npar + sum(is.finite(logprior[,g])) - 1
        }
      } else if (att.dist[g]=="higher.order") {
        if (higher.order$model == "Rasch") {
          npar <- npar + K
        } else if (higher.order$model == "2PL") {
          npar <- npar + 2 * K
        } else if (higher.order$model == "1PL") {
          npar <- npar + 1 + K
        }

      }
    }

    AIC <- 2 * npar + neg2LL

    BIC <- neg2LL + npar * log(N)

# cat("AIC=",AIC)
    item.prob <- vector("list",J)
initial.parm <- m2l(initial.parm)
    for (j in 1:J){
      item.prob[[j]] <- item.parm[j,1:Lj[j]]
      names(initial.parm[[j]]) <- names(item.prob[[j]]) <- paste0("P(",apply(alpha(Kj[j]),1,paste0,collapse = ""),")")
    }
    postP <- exp(t(logprior))
  if(sequential){ # transform LC.prob
    p <- LC.Prob
    for (j in 1:J){
      locj <- which(Qc[,1]==j)
      Qj <- Qc[locj,-c(1:2),drop=FALSE]
      if(nrow(Qj)>1){ # polytomous items
        pj <- sj <- rbind(LC.Prob[locj,],0)
        for (s in 1:(nrow(sj)-1)){
          pj[s,] <- apply(sj[1:s,,drop=FALSE],2,prod)*(1-sj[s+1,])
        }
        p[locj,] <- pj[-nrow(pj),]
      }else{ #dichotomous items
        p[locj,] <- LC.Prob[locj,]
      }

    }
    LC.Prob <- p
  }
    names(item.prob) <- names(initial.parm) <- rownames(LC.Prob) <- names(deltas) <- rownames(item.parm) <- item.names
    colnames(LC.Prob) <- colnames(postP) <- colnames(likepost$loglik) <-
      colnames(likepost$logpost) <- apply(alpha(K,T,Q),1,paste0,collapse = "")

    if(!is.null(group)) rownames(postP) <- paste("Group",gr.label)

    if(no.mg==1) {
      att.prior = c(exp(logprior0))
    }else{
      att.prior = exp(logprior0)
    }
    if(no.mg>1) lambda <- lam

    if (!is.null(oldseed)) .GlobalEnv$.Random.seed <- oldseed  else  rm(".Random.seed", envir = .GlobalEnv)

    s2 <- Sys.time()

    output <-
      list(
        catprob.parm = item.prob, delta.parm = deltas, catprob.matrix = item.parm,
        higher.order.struc.parm = lambda, model = model.names,LC.prob = LC.Prob,
        testfit = list(npar = npar,Deviance = neg2LL,AIC = AIC,BIC = BIC),
        posterior.prob = postP,
        technicals = list(logposterior.i = likepost$logpost, loglikelihood.i = likepost$loglik,
                          expectedCorrect = RN$Rg, expectedTotal = RN$Ng,initial.parm = initial.parm),
        options = list(timeused = s2 - s1,
                       start.time = s1, end.time = s2, dat = dat_d, Q = Qc, model = model,
                       itr = itr, dif.LL = dif.LL,dif.p=dif.p, npar = npar, item.npar= item.npar,
                       att.dist=att.dist, higher.order=higher.order,higher.order.model = higher.order$model,
                       mono.constraint = mono.constraint, item.names = item.names,
                       att.prior = att.prior, att.str=att.str,
                       nstarts = nstarts, conv.crit = conv.crit, maxitr = maxitr,
                       higher.order.method = higher.order$method, seq.dat = dat, no.group = no.mg,
                       group.label = gr.label,
                       verbose = verbose, catprob.parm = catprob.parm,sequential = sequential,
                       higher.order.struc.parm = higher.order$lambda, digits = digits,
                       diagnosis = diagnosis,call=GDINAcall),
        diagnos = list(itemprob.matrix = diag.itemprob, RN = diag.RN, likepost=diag.likepost,
                       changelog = diag.itrmaxchange, HO.parm = diag.lambda)
      )
    class(output) <- "GDINA"
    invisible(output)
  }
