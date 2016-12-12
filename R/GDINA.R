#' Calibrate dichotomous and polytomous response data
#'
#' \code{GDINA} calibrates the generalized deterministic inputs, noisy and
#' gate (G-DINA; de la Torre, 2011) model for dichotomous responses, and its extension, the sequential
#' G-DINA model (Ma, & de la Torre, 2016a), for ordinal and nominal responses.
#' By setting appropriate constraints, the deterministic inputs,
#' noisy and gate (DINA; de la Torre, 2009; Junker & Sijtsma, 2001) model,
#' the deterministic inputs, noisy or gate (DINO; Templin & Henson, 2006)
#' model, the reduced reparametrized unified model (R-RUM; Hartz, 2002),
#' the additive CDM (A-CDM; de la Torre, 2011), and the linear logistic
#' model (LLM; Maris, 1999) can also be calibrated. Different models can be fitted to different
#' items/categories in a single test. The attributes can be dichotomous or polytomous
#' (Chen & de la Torre, 2013). A higher-order structure
#' (de la Torre & Douglas, 2004) can be
#' assumed for all aforementioned models when attributes are binary.
#' Marginal maximum likelihood method with Expectation-Maximization (MMLE/EM) alogrithm
#' is used for item parameter estimation. Examinees' attributes are
#' estimated via the \code{MLE}, \code{EAP} and \code{MAP} methods.
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
#' distribution, by default, an empirical Bayes method (Carlin & Louis, 2000) is adopted, which is referred to as
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
#'  \deqn{P(\strong{\alpha}|\theta_i,\strong{\lambda})=\prod_k P(\alpha_k|\theta_i,\strong{\lambda})}
#' To estimate the parameters for higher order IRT model, either Bock and Aitkin's (1981) MMLE/EM algorithm or
#' Bock and Lieberman's (BL; 1970) method can be used.
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
#' responses of \eqn{N} examinees to \eqn{J} items. Missing values need to be coded as \code{NA}.
#' @param Q A required \eqn{J \times K} item/category and attribute association matrix, wher J represents the number of
#'    items/categories and K represents the number of attributes. For binary attributes,
#'    1 denotes attributes are measured by items and 0 means attributes are not
#'    necessary. For polytomous attributes, non-zero elements indicate which level
#'    of attributes are needed. Note that for polytomous items, the sequential G-DINA
#'    model is used and either restricted or unrestricted category-level Q-matrix is needed.
#'    The first column represents the item number and
#'    the second column indicates the category number. See \code{Examples}.
#' @param model A vector for each item/category or a scalar which will be used for all
#'    items/categories to specify which model is fitted to each item/category. The possible options
#'    include \code{"GDINA"},\code{"DINA"},\code{"DINO"},\code{"ACDM"},\code{"LLM"}, and \code{"RRUM"}.
#'    If \code{model} is a scalar, the specified model is fitted to all items. Different
#'    models can be assigned to different items or categories.
#'    It is also possible to specify models using numbers. Particularly, 0,1,2,3,4 and 5 represents
#'    \code{"GDINA"},\code{"DINA"},\code{"DINO"},\code{"ACDM"},\code{"LLM"}, and \code{"RRUM"}, respectively.
#' @param sequential logical; whether a sequential model is fitted for polytomous responses?
#' @param group a scalar indicating which column in \code{dat} is group indicator or
#'    a numerical vector indicating the group each individual belongs to. If it is a vector,
#'    its length must be equal to the number of individuals. More than two groups cannot be handled.
#' @param item.names A vector giving the name of items. If NULL (default), items are named as "Item 1", "Item 2", etc.
#' @param higher.order logical; \code{TRUE} indicates a higher order structure of attributes
#'    is assumed and higher order parameters will be estimated. Higher order model
#'    needs to be specified in \code{higher.order.model}. The default is \code{FALSE}.
#' @param higher.order.model An IRT model for higher order attribute structure; It can be
#'    \code{"2PL"}, \code{"1PL"} or \code{"Rasch"}, representing two parameter logistic IRT model,
#'    one parameter logistic IRT model and Rasch model,
#'    respectively. Please note that in \code{"1PL"} model, a common slope parameter will be
#'    estimated (see \code{Details}). \code{"Rasch"} is the default.
#' @param higher.order.method The algorithm for estimation the higher order parameters; it can be either
#'    \code{"MMLE"} using marginal maximum likelihood estimation, or \code{"BL"} based on the Bock and
#'    Lieberman's (1970) approach. \code{"BL"} is suitable when the number of attributes is few. It is not
#'    sensitive to sample size but can be very slow if the number of attributes is large.
#'    \code{"MMLE"}, which is the default, is suitable for most conditions but might be slow if
#'    sample size is extremely large.
#' @param mono.constraint logical; \code{TRUE} indicates that \eqn{P(\bm{\alpha}_1) <=P(\bm{\alpha}_2)} if
#'    for all \eqn{k}, \eqn{\alpha_{1k} <= \alpha_{2k}}. It can be a vector for each item or a scalar which will be used for all
#'    items to specify whether monotonicity constraint should be added for each item.
#' @param catprob.parm initial category probability parameters; It must be a list giving probability of success
#'    of each reduced latent classes for each non-zero category of each item.
#' @param verbose Print the max changes in item parameters and log-likelihood
#'     after each EM iteration or not? It can be 0, 1 or 2. If \code{verbose=0}, no information will
#'     be printed; if \code{verbose=1}, only information for current iteration will be shown; if \code{verbose=2},
#'     information for all iterations will be printed.
#' @param empirical Logical; whether empirical bayes is adopted or not? \code{TRUE} is
#'    the default when higher order attribute structure is not assumed. If estimating
#'    higher order structure, it will be \code{FALSE}.
#' @param att.prior attribute prior distribution for \eqn{2^K} latent classes. Only available for dichotomous attributes.
#'    It can be used to specify the hierarchical structure of attributes.
#'    Its length must be equal to \eqn{2^K}. Element 0 specifies which latent class does not exist.
#'    The sum of all elements does not have to be equal to 1; however, it will be standardized so that the sum is equal to 1
#'    before model calibration. When any latent class is given prior 0, standard errors for item parameters are not available.
#'    (1) If \code{empirical=FALSE} and \code{higher.order=FALSE}, the attribute prior distribution is fixed during model
#'    calibration; (2) if \code{empirical=TRUE} and \code{higher.order=FALSE}, the distribution for all latent classes
#'    with non-zero priors is updated using the empirical bayes method;
#'    (3) if \code{empirical=FALSE} and \code{higher.order=TRUE}, the distribution for all latent classes
#'    with non-zero priors is updated using the higher order model.
#'    The label for each latent class can be obtained by calling \code{attributepattern(K)}. See \code{examples} for more info.
#' @param att.str logical; whether attributes have any structure?
#' @param nstarts how many sets of starting values? The default is 1.
#' @param conv.crit The convergence criterion for max absolute change in item parameters.
#' @param conv.type How is the convergence criterion evaluated? It can be \code{"max.p.change"}, indicating
#'    the maximum absolute change in success probabilities is examined. It can also be \code{"dev.change"}, representing
#'    the absolute change in deviance is evaluated.
#' @param maxitr The maximum iterations of EM cycles allowed.
#' @param higher.order.struc.parm A matrix or data frame providing higher order structural parameters. If supplied, it must be of dimension \eqn{K\times 2}.
#'    The first column is the slope parameters and the second column is the intercept.
#' @param diagnosis Run in advanced mode? If it is 1 or 2, some intermediate results obtained in each iteration will be saved in \code{diagnos}.
#' @param Mstep.warning Logical; indicating whether the warning message in Mstep, if any, should be output immediately.
#' @param optimizer Advanced option: a string indicating which optimizer should be used in M-step.
#' @param lower.p fixed lower bound for success probabilities. If supplied as a single value, it will be assigned to all items;
#'   It can also be specified as a numeric vector corresponding to each item. The default is 1e-4 for all items.
#' @param upper.p fixed upper bound for success probabilities. If supplied as a single value, it will be assigned to all items;
#'   It can also be specified as a numeric vector corresponding to each item. The default is 0.9999 for all items.
#' @param smallNcorrection a numeric vector with two elements specifying the corrections applied when the expected number of
#' examinees in some latent groups are too small. Particularly, if the expected no. of examinees is less than the second element,
#' the first element and two times the first element will be added to the numerator and denominator of the closed-form solution of
#' probabilities of success.
#' @param randomseed random seed for generating initial item parameters. The default random seed is 123456.
#'
#' @author {Wenchao Ma, Rutgers University, \email{wenchao.ma@@rutgers.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @seealso See \code{\link{autoGDINA}} for Q-matrix validation, item level model comparison and model calibration
#' in one run; See \code{\link{itemfit}} for item fit analysis, \code{\link{Qval}} for Q-matrix validation,
#' \code{\link{modelcomp}} for item level model comparison and \code{\link{simGDINA}} for data simulation
#'
#' @return \code{GDINA} returns an object of class \code{GDINA}. S3 methods for \code{GDINA} objects
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
#' @importFrom graphics plot axis points text
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
#' itemparm(mod1, withSE = TRUE) # item probabilities of success + standard errors
#' itemparm(mod1, what = "delta") # delta parameters
#' itemparm(mod1, what = "delta",withSE=TRUE) # delta parameters
#' itemparm(mod1, what = "gs") # guess and slip parameters
#' itemparm(mod1, what = "gs",withSE = TRUE) # guess and slip parameters + standard errors
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
#' mod3 <- GDINA(dat = dat, Q = Q, model = "DINO")
#' #slip and guessing
#' itemparm(mod3, what = "gs") # guess and slip parameters
#' itemparm(mod3, what = "gs",withSE = TRUE) # guess and slip parameters + standard errors
#'
#' # Model comparison at test level via likelihood ratio test
#' anova(mod1,mod3)
#'
#' # --------- ACDM model -------#
#' mod4 <- GDINA(dat = dat, Q = Q, model = "ACDM")
#' mod4
#' # --------- LLM model -------#
#' mod4b <- GDINA(dat = dat, Q = Q, model = "LLM")
#' mod4b
#' # --------- RRUM model -------#
#' mod4c <- GDINA(dat = dat, Q = Q, model = "RRUM")
#' mod4c
#'
#' # --- Different CDMs for different items --- #
#'
#' models <- c(rep("GDINA",3),"LLM","DINA","DINO","ACDM","RRUM","LLM","RRUM")
#' mod5 <- GDINA(dat = dat, Q = Q, model = models)
#' anova(mod1,mod5)
#'
#'
#'####################################
#'#        Example 2.                #
#'#        Model estimations         #
#'# With monotonocity constraints    #
#'####################################
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
#'# With Higher order att structure  #
#'####################################
#'
#' # --- Higher order G-DINA model ---#
#' mod12 <- GDINA(dat = dat, Q = Q, model = "GDINA",
#'                higher.order = TRUE,higher.order.method="BL")
#' hoest=hoparm(mod12) # extract higher-order parameters
#' hoest$theta # ability
#' hoest$lambda # structural parameters
#' # --- Higher order DINA model ---#
#' mod22 <- GDINA(dat = dat, Q = Q, model = "DINA",
#'                higher.order = TRUE,higher.order.method="MMLE")
#' # --- Higher order DINO model ---#
#' mod23 <- GDINA(dat = dat, Q = Q, model = "DINO",higher.order = TRUE)
#' # --- Higher order ACDM model ---#
#' mod24 <- GDINA(dat = dat, Q = Q, model = "ACDM",
#'                higher.order = TRUE,higher.order.model="1PL")
#' # --- Higher order LLM model ---#
#' mod25 <- GDINA(dat = dat, Q = Q, model = "LLM",higher.order = TRUE)
#' # --- Higher order RRUM model ---#
#' mod26 <- GDINA(dat = dat, Q = Q, model = "RRUM",higher.order = TRUE)
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
#' # fit GDINA model  - empirical must be FALSE
#' modp1 <- GDINA(dat = dat, Q = Q, att.prior = prior, att.str = TRUE, empirical = FALSE)
#' # See the posterior weights
#' extract(modp1,what = "posterior.prob")
#'
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
#'  simD <- simGDINA(N,Q,gs.parm = gs,
#'                    model = "DINA",attribute = att)
#'  dat <- simD$dat
#'  # setting structure: A1 -> A2 -> A3
#'  # note: groups with prior 0 are assumed impossible;
#'  # and therefore their posteriors are not updated
#'  prior <- c(0.1,0.2,0,0,0.2,0,0,0.5)
#'  out <- GDINA(dat,Q,att.prior=prior,att.str = TRUE, model="DINA")
#'  # check posterior dist.
#'  extract(out,what = "posterior.prob")
#'
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
#' modp1 <- GDINA(dat = dat, Q = Q, att.prior = struc$att.prob, att.str = TRUE, empirical = FALSE)
#' modp1
#' # Note that calculated posterior below is fixed to the prior
#' extract(modp1,what = "posterior.prob")
#'
#' modp2 <- GDINA(dat = dat, Q = Q, att.prior = struc$att.prob, att.str = TRUE, empirical = TRUE)
#' modp2
#' # Note that calculated posterior below is similar but not fixed to the prior
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
#'
#'####################################
#'#           Example 7.             #
#'#        Model estimations         #
#'#          Without M-step          #
#'####################################
#'
#'  # -----------Fix User specified item parameters
#'  # Item parameters are not estimated
#'  # Only person attributes are estimated
#'  initials <- sim10GDINA$simItempar
#'  dat <- sim10GDINA$simdat
#'  Q <- sim10GDINA$simQ
#'  mod.fix <- GDINA(dat,Q,catprob.parm = initials,maxit = 0)
#'  itemparm(mod.fix)
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
#'
#' # Fit G-DINA model
#' mg.est <- GDINA(dat = dat,Q = Q,group = rep(c("female","male"),each=2000))
#' extract(mg.est,"posterior.prob")
#' }
#'
GDINA <-
  function(dat, Q, model = "GDINA", sequential = FALSE, item.names = NULL,
           higher.order = FALSE, higher.order.model =
             "Rasch", higher.order.method = "MMLE",
           verbose = 1, catprob.parm = NULL, higher.order.struc.parm = NULL,
           mono.constraint = FALSE, group = NULL,
           empirical = !higher.order, att.prior = NULL, att.str = FALSE,
           lower.p = 0.0001,upper.p = 0.9999, smallNcorrection = c(0.0005,0.001),
           nstarts = 1, conv.crit = 0.001, conv.type = "max.p.change",maxitr = 1000,
           digits = 4,diagnosis = 0,Mstep.warning = FALSE,optimizer = "all",
           randomseed = 123456)
  {
    s1 <- Sys.time()
    GDINAcall <- match.call()
    if(!is.null(group)){ # multiple groups
      # scalar -> group variable column
      if (length(group)==1){
        gr <- dat[,group] # group indicator variable - numeric
        dat <- dat[,-group] # responses
      }else{
        if (nrow(dat)!=length(group))stop("The length of group variable must be equal to the number of individuals.",call. = FALSE)
        gr <- group # group indicator variable
      }

      gr.label <- unique(gr) # group labels
      no.mg <- length(gr.label) # the number of groups
      if(!is.numeric(gr)||max(gr)>no.mg){
        for(g in 1:no.mg) gr[gr==gr.label[g]] <- g
        gr <- as.numeric(gr) # numeric variable
      }

    }else{
      no.mg <- 1L
      gr <- rep(1L,nrow(dat))
      gr.label <- "all data"
    }


    inputcheck(dat = dat, Q = Q, model = model, sequential = sequential, higher.order = higher.order,
               higher.order.model = higher.order.model, higher.order.method = higher.order.method,
               verbose = verbose, catprob.parm = catprob.parm, mono.constraint = mono.constraint,
               empirical = empirical, att.prior = att.prior,
               lower.p = lower.p, upper.p = upper.p, att.str = att.str, nstarts = nstarts,
               conv.crit = conv.crit, maxitr = maxitr, higher.order.parm = higher.order.struc.parm,
               digits = digits, diagnosis = diagnosis)

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
      if (is.null(item.names)) item.names <- paste("Item", 1:nrow(Q))

    }

    # missing value indicator 1 - nonmissing; 0 - missing
    missInd <- 1 - is.na(dat)
    # print(dim(missInd))
    if (any(missInd == 0)) {
      # some missings individuals with one or fewer valid response are
      # removed
      del.ind <- which(rowSums(missInd) < 2, arr.ind = TRUE)
      if (length(del.ind) > 0) {
        warning(length(del.ind), " individuals with one or fewer valid responses are removed.")
        dat <- dat[-del.ind, ]
        missInd <- missInd[-del.ind, ]
        gr <- gr[-del.ind]
      }
    }


    if (att.str && higher.order) stop("Higher order models cannot be estimated with structured attributes.", call. = FALSE)
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

    # location of item parameters
    parloc <- eta.loc(Q)  #J x L

    # Generate the log prior - a L x no.mg matrix
    if (is.null(att.prior)) {
      att.prior <- rep(1/L, L)
      logprior <- matrix(log(att.prior),nrow = L,ncol = no.mg)
    } else {
      att.prior <- matrix(att.prior,ncol = no.mg)
      if (nrow(att.prior) != L)
        stop("Length of joint attribute distribution prior must be equal to 2^K if it is specified.",
             call. = FALSE)
      if (any(att.prior < 0))
        stop("Joint attribute distribution prior can only contain 0 or positive numbers.", call. = FALSE)
        logprior <- log(att.prior/matrix(colSums(att.prior),nrow = L,ncol = no.mg,byrow = TRUE))
    }

    # generating initial item parameters
    if (is.null(catprob.parm)) {
      item.parm <- initials(Q, nstarts, randomseed = randomseed)  #a list with nstarts matrices of size J x 2^Kjmax
      ###### Multiple starting values
      if (nstarts > 1) {
        neg2LL <- NULL
        for (i in 1:nstarts) {
          LC.Prob <- uP(as.matrix(parloc), as.matrix(item.parm[[i]]))
          # calculate likelihood and posterior
          neg2LL <- c(neg2LL, -2 * Lik(mP = LC.Prob, mX = as.matrix(dat),
                                       vlogPrior = as.matrix(logprior), mIndmiss = as.matrix(missInd),vgroup = gr)$LL)
        }
        item.parm <- item.parm[[which.min(neg2LL)]]

      }
    } else {
      item.parm <- l2m(catprob.parm)
    }

    dif <- 10
    itr <- 0
    diag.itemprob <- diag.RN <- diag.likepost <- diag.lambda <- diag.opts <- vector("list", maxitr)
    diag.itrmaxchange <- NULL
    deltas <- calc_delta(item.parm, model, Kj)
    # print(deltas)
    lambda <- higher.order.struc.parm
    # cat('\nIteration maxchange -2LL\n') E-M algorithm
    dif.p <- dif.LL <- LL.1 <- LL.2 <- 0

    ############### START of while loop ##########
    designMlist <- vector("list",J)
    for(j in 1:J){ #for each item
      designMlist[[j]] <- designmatrix(Kj[j],model[j])
    }
    while (dif > conv.crit && itr < maxitr)
    {
      item.parm_copy <- item.parm
      # --------E step probability of getting each score for each latent
      # class
      # print(item.parm)
      LC.Prob <- uP(as.matrix(parloc), as.matrix(item.parm))
      LC.Prob[is.nan(LC.Prob)] <- 0
      # calculate likelihood and posterior
      # vlogPrior: col vector for single group; matrix with g columns for g groups;
      # g: group indicator from 1 to g;
      likepost <- Lik(mP = LC.Prob, mX = as.matrix(dat),
                      vlogPrior = as.matrix(logprior), mIndmiss = as.matrix(missInd),vgroup = gr)
      # expected # of examinees of answering items correctly in each latent
      # class expected # of examinees for each latent class
      RN <- NgRg(mlogPost = likepost$logpost, mX = as.matrix(dat), mloc = as.matrix(parloc), mIndmiss = as.matrix(missInd))
      # print(RN)
      LL.1 <- -2*likepost$LL
      dif.LL <- abs(LL.1-LL.2)
      if(tolower(conv.type)=="dev.change") {
        dif <- dif.LL
        if (dif.LL<=conv.crit) break
      }
      LL.2 <- LL.1
      if (higher.order)
      {

        Rl <- c(exp(likepost$logprior) * N)
        a <- b <- NULL
        if (!is.null(lambda)) a <- lambda[,1];b <- lambda[,2]
        HO.out <- HO.est(Rl, likepost$loglik, K = K, N = N, IRTmodel = higher.order.model,
                         method = higher.order.method, a = a, b = b)
        HOlogprior <- HO.out$logprior
        lambda <- HO.out$lambda
        colnames(lambda) <- c("slope", "intercept")
        rownames(lambda) <- paste("Attribute", 1:K, sep = " ")
        #print(HO.out$lambda)
        if (diagnosis >= 1) diag.lambda[[itr+1]] <- lambda
      }



      if (att.str) smallNcorrection <- c(-1,-1)
      optims <- Mstep(Kj=Kj, RN=RN, model=model, itmpar=item.parm, delta=deltas,
                      constr=mono.constraint,correction = smallNcorrection, lower.p = lower.p,
                      upper.p = upper.p, itr=itr+1,warning.immediate=Mstep.warning,
                      optimizer = optimizer, designmatrices = designMlist)
      item.parm <- optims$item.parm
      deltas <- optims$delta

      dif.p <- max(abs(item.parm - item.parm_copy),na.rm = TRUE)
      if(tolower(conv.type)=="max.p.change") {
        dif <- dif.p
      }
      if (diagnosis>=1) {
        diag.itemprob[[itr+1]] <- item.parm_copy
        diag.itrmaxchange <- rbind(diag.itrmaxchange,round(c(dif,-2*likepost$LL),6))
        if (diagnosis>=2) {
          diag.likepost[[itr+1]] <- likepost
          diag.RN[[itr+1]] <- RN
          diag.opts[[itr+1]] <- optims$optims
        }

      }

      itr <- itr + 1
      if(verbose==1) {
        cat('\rIteration =',itr,' Max change =',format(round(dif,4),nsmall=4),' Deviance =',format(round(-2*likepost$LL,2),nsmall=2))
      }else if (verbose==2) {
        cat('Iteration =',itr,' Max change =',format(round(dif,4),nsmall=4),' Deviance =',format(round(-2*likepost$LL,2),nsmall=2),"\n")
      }

      if (higher.order == FALSE & empirical == TRUE)
      {
        logprior <- likepost$logprior
      }else if (higher.order == FALSE & empirical == FALSE)
      {
        logprior <- log(att.prior)
      }else if(higher.order==TRUE){
        logprior <- HOlogprior
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
    # calculate likelihood and posterior

    if (higher.order == FALSE & empirical == TRUE) {
      likepost <- Lik(mP = LC.Prob, mX = as.matrix(dat),
                      vlogPrior = as.matrix(logprior), mIndmiss = as.matrix(missInd),vgroup = gr)
      logprior <- likepost$logprior
    } else if (higher.order == FALSE & empirical == FALSE) {
      logprior <- log(att.prior)
    } else if (higher.order == TRUE) {
      # logprior <- HOlogprior
      # if (att.str) stop("Joint attribute distribution cannot be structured higher-order.",call. = FALSE)
    }
    likepost <- Lik(mP = LC.Prob, mX = as.matrix(dat),
                    vlogPrior = as.matrix(logprior), mIndmiss = as.matrix(missInd),vgroup = gr)
    RN <- NgRg(mlogPost = likepost$logpost, mX = as.matrix(dat), mloc = as.matrix(parloc), mIndmiss = as.matrix(missInd))
    # -----------------Test Fit information----------------#

    neg2LL <- -2 * likepost$LL

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
    if (higher.order == FALSE) {
      if (!att.str) {
        npar <- npar + L - 1
      } else {
        npar <- npar + sum(is.finite(logprior)) - 1
      }
    } else if (higher.order == TRUE) {
      if (higher.order.model == "Rasch") {
        npar <- npar + K
      } else if (higher.order.model == "2PL") {
        npar <- npar + 2 * K
      } else if (higher.order.model == "1PL") {
        npar <- npar + 1 + K
      }

    }
    AIC <- 2 * npar + neg2LL

    BIC <- neg2LL + npar * log(N)


    item.prob <- vector("list",J)

    for (j in 1:J){
      item.prob[[j]] <- item.parm[j,1:Lj[j]]
      names(item.prob[[j]]) <- paste0("P(",apply(alpha(Kj[j]),1,paste0,collapse = ""),")")
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
    names(item.prob) <- rownames(LC.Prob) <- names(deltas) <- rownames(item.parm) <- item.names
    colnames(LC.Prob) <- colnames(postP) <- colnames(likepost$loglik) <-
      colnames(likepost$logpost) <- apply(alpha(K,T,Q),1,paste0,collapse = "")

    if(!is.null(group)) rownames(postP) <- paste("Group",gr.label)

    s2 <- Sys.time()

    output <-
      list(
        catprob.parm = item.prob, delta.parm = deltas, catprob.matrix = item.parm,
        higher.order.struc.parm = lambda, model = model.names,LC.prob = LC.Prob,
        testfit = list(npar = npar,Deviance = round(neg2LL,2),AIC = round(AIC,2),BIC = round(BIC,2)),
        posterior.prob = postP,
        technicals = list(logposterior.i = likepost$logpost, loglikelihood.i = likepost$loglik,
                          expectedCorrect = RN$Rg, expectedTotal = RN$Ng),
        options = list(timeused = s2 - s1,
                       start.time = s1, end.time = s2, dat = dat_d, Q = Qc, model = model,
                       itr = itr, dif.LL = dif.LL,dif.p=dif.p, npar = npar, item.npar= item.npar,
                       higher.order = higher.order, higher.order.model = higher.order.model,
                       mono.constraint = mono.constraint, item.names = item.names,
                       empirical = empirical, att.prior = att.prior, att.str=att.str,
                       nstarts = nstarts, conv.crit = conv.crit, maxitr = maxitr,
                       higher.order.method = higher.order.method, seq.dat = dat,
                       verbose = verbose, catprob.parm = catprob.parm,sequential = sequential,
                       higher.order.struc.parm = higher.order.struc.parm, digits = digits,
                       diagnosis = diagnosis,call=GDINAcall),
        diagnos = list(itemprob.matrix = diag.itemprob, RN = diag.RN, likepost=diag.likepost,
                       changelog = diag.itrmaxchange, HO.parm = diag.lambda)
      )
    class(output) <- "GDINA"
    invisible(output)
  }
