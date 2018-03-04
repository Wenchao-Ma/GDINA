#' @title Data simulation based on the G-DINA models
#'
#' @description
#'    Simulate responses based on the G-DINA model (de la Torre, 2011) and sequential G-DINA model
#'    (Ma & de la Torre, 2016), or CDMs subsumed by them, including the DINA model, DINO model, ACDM,
#'    LLM and R-RUM. Attributes can be simulated from uniform, higher-order or multivariate normal
#'    distributions, or be supplied by users. See \code{Examples} and \code{Details} for
#'    how item parameter specifications. See the help page of \code{\link{GDINA}}
#'    for model parameterizations.
#'
#' @details
#' Item parameter specifications in \code{simGDINA}:
#'
#' Item parameters can be specified in one of three different ways.
#'
#' The first and probably the easiest way is to specify the guessing and slip parameters for each item or nonzero category using
#' \code{gs.parm}, which is a matrix or data frame for \eqn{P(\bm{\alpha}_{lj}^*=0)} and \eqn{1-P(\bm{\alpha}_{lj}^*=1)}
#' for all items for dichotomous items and \eqn{S(\bm{\alpha}_{ljh}^*=0)} and \eqn{1-S(\bm{\alpha}_{ljh}^*=1)}
#' for all nonzero categories for polytomous items. Note that \eqn{1-P(\bm{\alpha}_{lj}^*=0)-P(\bm{\alpha}_{lj}^*=1)} or
#' \eqn{1-S(\bm{\alpha}_{lj}^*=0)-S(\bm{\alpha}_{lj}^*=1)} must be greater than 0.
#' For generating ACDM, LLM, and RRUM, delta parameters are generated randomly if \code{type="random"},
#' or in a way that each required attribute contributes equally, as in
#'  Ma, Iaconangelo, & de la Torre (2016) if \code{type="equal"}. For ACDM, LLM and RRUM, generated
#'  delta parameters are always positive, which implies that monotonicity constraints are always satisfied.
#'  If the generating model is the G-DINA model, \code{mono.constraint} can be used to specify whether monotonicity
#'  constraints should be satisfied.
#'
#' The second way of simulating responses is to specify success probabilities (i.e., \eqn{P(\bm{\alpha}_{lj}^*)}
#' or \eqn{S(\bm{\alpha}_{ljh}^*)}) for each nonzero category of each item directly
#' using the argument \code{catprob.parm}. If an item or category requires \eqn{K_j^*} attributes, \eqn{2^{K_j^*}} success probabilities
#' need to be provided. \code{catprob.parm} must be a list, where each element gives the success probabilities for nonzero category of each item.
#' Note that success probabilities cannot be negative or greater than one.
#'
#' The third way is to specify delta parameters for data simulation. For DINA and DINO model, each nonzero category requires two
#' delta parameters. For ACDM, LLM and RRUM, if a nonzero category requires \eqn{K_j^*} attributes, \eqn{K_j^*+1} delta parameters
#' need to be specified. For the G-DINA model, a nonzero category requiring \eqn{K_j^*} attributes has \eqn{2^{K_j^*}} delta parameters.
#' It should be noted that specifying delta parameters needs to ascertain the derived success probabilities are within the \eqn{[0,1]} interval.
#'
#' Please note that you need to specify item parameters in ONLY one of these three ways. If \code{gs.parm} is specified, it will be used regardless of
#' the inputs in \code{catprob.parm} and \code{delta.parm}. If \code{gs.parm} is not specified, \code{simGDINA} will check
#' if \code{delta.parm} is specified; if yes, it will be used for data generation. if both \code{gs.parm} and \code{delta.parm} are not specified,
#' \code{catprob.parm} is used for data generation.
#'
#'
#' @param N Sample size.
#' @param Q A required \eqn{J \times K} item/category and attribute association matrix, wher J represents the number of
#'    items/categories and K represents the number of attributes. For binary attributes,
#'    1 denotes attributes are measured by items and 0 means attributes are not
#'    necessary. For polytomous attributes, non-zero elements indicate which level
#'    of attributes are needed. Note that for polytomous items, the sequential G-DINA
#'    model is used and either restricted or unrestricted category-level Q-matrix is needed.
#'    The first column represents the item number and
#'    the second column indicates the category number. See \code{Examples}.
#' @param gs.parm A matrix or data frame for guessing and slip parameters. It
#'    must be of dimension \eqn{J \times 2}, where the first column represents the guessing parameters (or \eqn{P(0)}),
#'    and the second column represents slip parameters (or \eqn{1-P(1)}). This needs to be used in conjunction with
#'    the argument \code{type} if generating models include ACDM, LLM, or RRUM, and \code{model}.
#' @param type How are the delta parameters for ACDM, LLM, RRUM generated?
#'    It can be either \code{"random"} or \code{"equal"}. \code{"random"} means the delta parameters are simulated randomly,
#'    while \code{"equal"} means that each required attribute contributes equally to the probability of success (P), logit(P) or
#'    log(P) for ACDM, LLM and RRUM, respectively. See \code{Details} for more information.
#' @param model A vector for each item/category or a scalar which will be used for all
#'    items/categories to specify which model is fitted to each item/category. The possible options
#'    include \code{"GDINA"},\code{"DINA"},\code{"DINO"},\code{"ACDM"},\code{"LLM"}, and \code{"RRUM"}.
#'    If \code{model} is a scalar, the specified model is fitted to all items. Different
#'    models can be assigned to different items or categories.
#' @param sequential logical; \code{TRUE} if the sequential model is used for polytomous responses simulation, and \code{FALSE}
#'    if there is no polytomously scored items.
#' @param mono.constraint A vector for each item/category or a scalar which will be used for all
#'    items/categories to specify whether monotonicity constraints should be satisfied if the generating model is the G-DINA model. Note that
#'    this is applicable only for the G-DINA model when \code{gs.parm} is used. For ACDM, LLM and RRUM, monotonicity constraints
#'    are always satisfied and therefore this argument is ignored.
#' @param catprob.parm A list of success probabilities for each latent group for each non-zero category of each item. See \code{Examples} and
#'    \code{Details} for more information.
#' @param delta.parm A list of delta parameters for each latent group for each item or category.
#' @param item.names A vector giving the name of items or categories. If it is \code{NULL} (default), items are named as "Item 1", "Item 2", etc.
#' @param attribute optional user-specified person attributes. It is a \eqn{N\times K} matrix or data frame. If this is not supplied, attributes are simulated
#'    from a distribution specified in \code{att.dist}.
#' @param att.dist A string indicating the distribution for attribute simulation. It can be \code{"uniform"}, \code{"higher.order"},
#'    \code{"mvnorm"} or \code{"multinomial"} for uniform, higher-order, multivariate normal and multinomial distributions, respectively.
#'    The default is the uniform distribution. To specify structural parameters for the higher-order
#'    and multivariate normal distributions, see \code{higher.order.parm} and \code{mvnorm.parm}, respectively. To specify the probabilities
#'    for the multinomial distribution, use \code{att.prior} argument.
#'@param att.prior probability for each attribute pattern. Order is the same as that returned from \code{attributepattern(Q = Q)}. This is only
#'    applicable when \code{att.dist="fixed"}.
#' @param higher.order.parm A list specifying parameters for higher-order distribution for attributes
#'    if \code{att.dist=higher.order}. Particularly, \code{theta} is a
#'    vector of length \eqn{N} representing the higher-order ability
#'    for each examinee. and \code{lambda} is a \eqn{K \times 2} matrix. Column 1 gives the slopes for the higher-order
#'    model and column 2 gives the intercepts. See \code{\link{GDINA}} for the formulations of the higher-order
#'    models.
#' @param mvnorm.parm a list of parameters for multivariate normal attribute distribution. \code{mean} is a vector of length \eqn{K}
#'    specifying the mean of multivariate normal distribution; and \code{sigma} is a positive-definite
#'    symmetric matrix specifying the variance-covariance matrix. \code{cutoffs} is a vector giving the
#'    cutoff for each attribute. See \code{Examples}.
#' @param digits How many decimal places in each number? The default is 4.
#' @return an object of class \code{simGDINA}. Elements that can be extracted using method \code{extract}
#' include:
#' \describe{
#' \item{dat}{simulated item response matrix}
#' \item{Q}{Q-matrix}
#' \item{attribute}{A \eqn{N \times K} matrix for inviduals' attribute patterns}
#' \item{catprob.parm}{a list of non-zero category success probabilities for each latent group}
#' \item{delta.parm}{a list of delta parameters}
#' \item{higher.order.parm}{Higher-order parameters}
#' \item{mvnorm.parm}{multivariate normal distribution parameters}
#' \item{LCprob.parm}{A matrix of item/category success probabilities for each latent class}
#' }
#'
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @name simGDINA
#'
#' @export
#'
#' @references
#'
#' Chiu, C.-Y., Douglas, J. A., & Li, X. (2009). Cluster analysis for cognitive diagnosis: Theory and applications. \emph{Psychometrika, 74}, 633-665.
#'
#' de la Torre, J. (2011). The generalized DINA model framework. \emph{Psychometrika, 76}, 179-199.
#'
#' de la Torre, J., & Douglas, J. A. (2004). Higher-order latent trait models for cognitive diagnosis. \emph{Psychometrika, 69}, 333-353.
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
#' Ma, W., & de la Torre, J. (2016). A sequential cognitive diagnosis model for polytomous responses. \emph{British Journal of Mathematical and Statistical Psychology. 69,} 253-275.
#'
#' Ma, W., Iaconangelo, C., & de la Torre, J. (2016). Model similarity, model selection and attribute classification. \emph{Applied Psychological Measurement, 40}, 200-217.
#'
#' Maris, E. (1999). Estimating multiple classification latent class models. \emph{Psychometrika, 64}, 187-212.
#'
#' Templin, J. L., & Henson, R. A. (2006). Measurement of psychological disorders using cognitive diagnosis models.
#' \emph{Psychological Methods, 11}, 287-305.
#'
#' @examples
#'\dontrun{
#'####################################################
#'#                     Example 1                    #
#'#             Data simulation (DINA)               #
#'####################################################
#'N <- 500
#'Q <- sim30GDINA$simQ
#'J <- nrow(Q)
#'gs <- data.frame(guess=rep(0.1,J),slip=rep(0.1,J))
#'
#'# Simulated DINA model; to simulate G-DINA model
#'# and other CDMs, change model argument accordingly
#'
#'sim <- simGDINA(N,Q,gs.parm = gs,model = "DINA")
#'
#'# True item success probabilities
#'extract(sim,what = "catprob.parm")
#'
#'# True delta parameters
#'extract(sim,what = "delta.parm")
#'
#'# simulated data
#'extract(sim,what = "dat")
#'
#'# simulated attributes
#'extract(sim,what = "attribute")
#'
#'####################################################
#'#                     Example 2                    #
#'#             Data simulation (RRUM)               #
#'####################################################
#'N <- 500
#'Q <- sim30GDINA$simQ
#'J <- nrow(Q)
#'gs <- data.frame(guess=rep(0.2,J),slip=rep(0.2,J))
#'# Simulated RRUM
#'# deltas except delta0 for each item will be simulated
#'# randomly subject to the constraints of RRUM
#'sim <- simGDINA(N,Q,gs.parm = gs,model = "RRUM")
#'
#'# simulated data
#'extract(sim,what = "dat")
#'
#'# simulated attributes
#'extract(sim,what = "attribute")
#'
#'####################################################
#'#                     Example 3                    #
#'#             Data simulation (LLM)                #
#'####################################################
#'N <- 500
#'Q <- sim30GDINA$simQ
#'J <- nrow(Q)
#'gs <- data.frame(guess=rep(0.1,J),slip=rep(0.1,J))
#'# Simulated LLM
#'# By specifying type="equal", each required attribute is
#'# assumed to contribute to logit(P) equally
#'sim <- simGDINA(N,Q,gs.parm = gs,model = "LLM",type="equal")
#' #check below for what the equal contribution means
#'extract(sim,what = "delta.parm")
#'
#'# simulated data
#'extract(sim,what = "dat")
#'
#'# simulated attributes
#'extract(sim,what = "attribute")
#'
#'####################################################
#'#                   Example 4                      #
#'#          Data simulation (all CDMs)              #
#'####################################################
#'
#' set.seed(12345)
#'
#'N <- 500
#'Q <- sim10GDINA$simQ
#'J <- nrow(Q)
#'gs <- data.frame(guess=rep(0.1,J),slip=rep(0.1,J))
#'# Simulated different CDMs for different items
#'models <- c("GDINA","DINO","DINA","ACDM","LLM","RRUM","GDINA","LLM","RRUM","DINA")
#'sim <- simGDINA(N,Q,gs.parm = gs,model = models,type="random")
#'
#'# simulated data
#'extract(sim,what = "dat")
#'
#'# simulated attributes
#'extract(sim,what = "attribute")
#'
#'####################################################
#'#                   Example 5                      #
#'#          Data simulation (all CDMs)              #
#'#  using probability of success in list format     #
#'####################################################
#'
#' # success probabilities for each item need to be provided in list format as follows:
#' # if item j requires Kj attributes, 2^Kj success probabilities
#' # need to be specified
#' # e.g., item 1 only requires 1 attribute
#' # therefore P(0) and P(1) should be specified;
#' # similarly, item 10 requires 3 attributes,
#' # P(000),P(100),P(010)...,P(111) should be specified;
#' # the latent class represented by each element can be obtained
#' # by calling attributepattern(Kj)
#'itemparm.list <- list(item1=c(0.2,0.9),
#'                     item2=c(0.1,0.8),
#'                     item3=c(0.1,0.9),
#'                     item4=c(0.1,0.3,0.5,0.9),
#'                     item5=c(0.1,0.1,0.1,0.8),
#'                     item6=c(0.2,0.9,0.9,0.9),
#'                     item7=c(0.1,0.45,0.45,0.8),
#'                     item8=c(0.1,0.28,0.28,0.8),
#'                     item9=c(0.1,0.4,0.4,0.8),
#'                     item10=c(0.1,0.2,0.3,0.4,0.4,0.5,0.7,0.9))
#'set.seed(12345)
#' N <- 500
#' Q <- sim10GDINA$simQ
#' # When simulating data using catprob.parm argument,
#' # it is not necessary to specify model and type
#'sim <- simGDINA(N,Q,catprob.parm = itemparm.list)
#'
#'
#'####################################################
#'#                   Example 6                      #
#'#            Data simulation (all CDMs)            #
#'#      using delta parameters in list format       #
#'####################################################
#'
#' delta.list <- list(c(0.2,0.7),
#'                     c(0.1,0.7),
#'                     c(0.1,0.8),
#'                     c(0.1,0.7),
#'                     c(0.1,0.8),
#'                     c(0.2,0.3,0.2,0.1),
#'                     c(0.1,0.35,0.35),
#'                     c(-1.386294,0.9808293,1.791759),
#'                     c(-1.609438,0.6931472,0.6),
#'                     c(0.1,0.1,0.2,0.3,0.0,0.0,0.1,0.1))
#'
#' model <- c("GDINA","GDINA","GDINA","DINA","DINO","GDINA","ACDM","LLM","RRUM","GDINA")
#' N <- 500
#' Q <- sim10GDINA$simQ
#' # When simulating using delta.parm argument, model needs to be
#' # specified
#' sim <- simGDINA(N,Q,delta.parm = delta.list, model = model)
#'
#'
#'####################################################
#'#                   Example 7                      #
#'#      Data simulation (higher order DINA model)   #
#'####################################################
#'
#' Q <- sim30GDINA$simQ
#' gs <- matrix(0.1,nrow(Q),2)
#' N <- 500
#' set.seed(12345)
#' theta <- rnorm(N)
#' K <- ncol(Q)
#' lambda <- data.frame(a=rep(1,K),b=seq(-2,2,length.out=K))
#' sim <- simGDINA(N,Q,gs.parm = gs, model="DINA", att.dist = "higher.order",
#'                  higher.order.parm = list(theta = theta,lambda = lambda))
#'
#'####################################################
#'#                   Example 8                      #
#'#      Data simulation (higher-order CDMs)         #
#'####################################################
#'
#' Q <- sim30GDINA$simQ
#' gs <- matrix(0.1,nrow(Q),2)
#' models <- c(rep("GDINA",5),
#'             rep("DINO",5),
#'             rep("DINA",5),
#'             rep("ACDM",5),
#'             rep("LLM",5),
#'             rep("RRUM",5))
#' N <- 500
#' set.seed(12345)
#' theta <- rnorm(N)
#' K <- ncol(Q)
#' lambda <- data.frame(a=runif(K,0.7,1.3),b=seq(-2,2,length.out=K))
#' sim <- simGDINA(N,Q,gs.parm = gs, model=models, att.dist = "higher.order",
#'                  higher.order.parm = list(theta = theta,lambda = lambda))
#'
#'
#'####################################################
#'#                   Example 9                      #
#'#      Data simulation (higher-order model)        #
#'#  using the multivariate normal threshold model   #
#'####################################################
#'
#'
#' # See Chiu et al., (2009)
#'
#' N <- 500
#' Q <- sim10GDINA$simQ
#' K <- ncol(Q)
#' gs <- matrix(0.1,nrow(Q),2)
#' cutoffs <- qnorm(c(1:K)/(K+1))
#' m <- rep(0,K)
#' vcov <- matrix(0.5,K,K)
#' diag(vcov) <- 1
#' simMV <- simGDINA(N,Q,gs.parm = gs, att.dist = "mvnorm",
#'                  mvnorm.parm=list(mean = m, sigma = vcov,cutoffs = cutoffs))
#'
#'####################################
#'#          Example 10              #
#'#        Simulation using          #
#'#      user-specified att structure#
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
#' gs <- matrix(rep(0.1,2*nrow(Q)),ncol=2)
#' # data simulation
#' simD <- simGDINA(N,Q,gs.parm = gs,
#'                    model = "DINA",attribute = true.att)
#'
#'
#'####################################################
#'#                   Example 11                     #
#'#                Data simulation                   #
#'#  (GDINA with monotonicity constraints)           #
#'####################################################
#'
#' set.seed(12345)
#'
#'N <- 500
#'Q <- sim30GDINA$simQ
#'J <- nrow(Q)
#'gs <- data.frame(guess=rep(0.1,J),slip=rep(0.1,J))
#'# Simulated different CDMs for different items
#'sim <- simGDINA(N,Q,gs.parm = gs,model = "GDINA",mono.constraint=TRUE)
#'
#'# True item success probabilities
#'extract(sim,what = "catprob.parm")
#'
#'# True delta parameters
#'extract(sim,what = "delta.parm")
#'
#'# simulated data
#'extract(sim,what = "dat")
#'
#'# simulated attributes
#'extract(sim,what = "attribute")
#'
#'####################################################
#'#                   Example 12                     #
#'#                Data simulation                   #
#'# (Sequential G-DINA model - polytomous responses) #
#'####################################################
#'
#' set.seed(12345)
#'
#'N <- 2000
#'# restricted Qc matrix
#'Qc <- sim20seqGDINA$simQ
#'#total number of categories
#'J <- nrow(Qc)
#'gs <- data.frame(guess=rep(0.1,J),slip=rep(0.1,J))
#'# simulate sequential DINA model
#'simseq <- simGDINA(N, Qc, sequential = TRUE, gs.parm = gs, model = "DINA")
#'
#'# True item success probabilities
#'extract(simseq,what = "catprob.parm")
#'
#'# True delta parameters
#'extract(simseq,what = "delta.parm")
#'
#'# simulated data
#'extract(simseq,what = "dat")
#'
#'# simulated attributes
#'extract(simseq,what = "attribute")
#'}
#'
#'####################################################
#'#                   Example 13
#'#         DINA model Attribute generated using
#'#             multinomial distribution
#'####################################################
#'
#' Q <- sim10GDINA$simQ
#' gs <- matrix(0.1,nrow(Q),2)
#' N <- 5000
#' set.seed(12345)
#' prior <- c(0.1,0.2,0,0,0.2,0,0,0.5)
#' sim <- simGDINA(N,Q,gs.parm = gs, model="DINA", att.dist = "multinomial",att.prior = prior)
#' # check latent class sizes
#' table(sim$att.group)/N
#'
#'
simGDINA <- function(N, Q, gs.parm=NULL, model = "GDINA", sequential = FALSE, type = "random",
                      catprob.parm = NULL, delta.parm = NULL,mono.constraint = TRUE,
                      attribute = NULL, att.dist = "uniform", item.names = NULL,
                      higher.order.parm=list(theta = NULL, lambda = NULL),
                      mvnorm.parm=list(mean = NULL,sigma = NULL,cutoffs = NULL),
                     att.prior = NULL, digits=4)
  {
  simGDINAcall <- match.call()
  catprob.parm <- catprob.parm
  inputcheck.sim(N=N, Q=Q, sequential=sequential, gs.parm=gs.parm, model = model, type = type,
                             catprob.parm = catprob.parm, delta.parm = delta.parm)

  if(sequential){
    C <- table(Q[,1])
    S <- length(unique(Q[,1]))
    if (is.null(item.names)) item.names <- paste("Item",Q[,1],"Cat",Q[,2])
    Q <- Q[,-c(1:2)]
  }else{
    if (is.null(item.names)) item.names <- paste("Item",1:nrow(Q))
  }
  J <- nrow(Q)
  K <- ncol(Q)
  pattern <- attributepattern(Q = Q)
  pattern.t <- t(pattern)
  model <- model.transform(model,J)
  L <- nrow(pattern)  # the number of latent groups
  Kj <- apply(Q,1,function(x){sum(x>0)})  # The number of attributes for each item
  Kjmax <- max(Kj) # the maximum attributes required for each item
  catprob.matrix <- matrix(NA,J,2^Kjmax)
  par.loc <- eta(as.matrix(Q))
if (length(mono.constraint)==1)  mono.constraint <- rep(mono.constraint,J)

######################################################################################
  #
  #   Item parameter generation
  #
######################################################################################

######### based on gs.parm for easy item parameter simulation
if (!is.null(gs.parm)) {
  if(nrow(gs.parm)!=nrow(Q)) stop("The number of rows in gs is not equal to the number of non-zero categories.",call. = FALSE)
  if(any(1-rowSums(gs.parm)<0)) stop("Data cannot be simulated because 1-s-g<0 for some items - check your parameters or specify parameters using delta.parm or catprob.parm.",call. = FALSE)
  pd <- gs2p(Q=Q,gs=gs.parm,model=model,type=type,mono.constraint=mono.constraint,digits=8)
  delta.parm <- pd$delta.parm
  catprob.parm <- pd$itemprob.parm
  catprob.matrix <- pd$itemprob.matrix
}else if(!is.null(delta.parm))
  {
  catprob.parm <- vector("list",J)

    for (j in 1:J){

        Mj <- designmatrix(Kj[j],model[j])

        if (model[j]<=3){
          catprob.matrix[j,1:nrow(Mj)] <- catprob.parm[[j]] <- round(c(Mj%*%delta.parm[[j]]),digits)
        }else if (model[j]==4){
          catprob.matrix[j,1:nrow(Mj)] <- catprob.parm[[j]] <- round(plogis(c(Mj%*%delta.parm[[j]])),digits)
        }else if (model[j]==5){
          catprob.matrix[j,1:nrow(Mj)] <- catprob.parm[[j]] <-round(exp(c(Mj%*%delta.parm[[j]])),digits)
        }
        names(catprob.parm[[j]]) <- paste("P(",apply(attributepattern(Kj[j]),1,function(x){paste(x,collapse = "")}),")",sep = "")
    }
  delta.parm <- format_delta(delta.parm,model,Kj,digits=digits)
  }else if(!is.null(catprob.parm)){
    delta.parm <- vector("list",J)

      for (j in 1:J){
        catprob.matrix[j,1:length(catprob.parm[[j]])] <- round(c(catprob.parm[[j]]),digits)
        Mj <- designmatrix(Kj[j],model[j])
        delta.parm[[j]] <- round(c(solve(t(Mj)%*%Mj)%*%t(Mj)%*%c(catprob.parm[[j]])),digits)
        names(catprob.parm[[j]]) <- paste("P(",apply(attributepattern(Kj[j]),1,function(x){paste(x,collapse = "")}),")",sep = "")
      }
    delta.parm <- format_delta(delta.parm,model,Kj,digits=digits)
    }

  ######################################################################################
  #
  #   Person parameter generation
  #
  ######################################################################################

  if (is.null(attribute))
  {
    if (tolower(att.dist) == "uniform")
    {
      # uniform distribution
      att.group <- sample.int(L, N, replace = T)
    } else if (tolower(att.dist) == "higher.order")
    {
      if (max(Q) > 1) stop("Higher order structure is not allowed currently when attributes are polytomous.",call. = FALSE)

      if (is.null(higher.order.parm$theta)||is.null(higher.order.parm$lambda))
      {
        stop("Higher-order parameters must be provided.",call. = FALSE)
      }
      a <- higher.order.parm$lambda[, 1]
      b <- higher.order.parm$lambda[, 2]
      att <- (Pr_2PL_vec(theta = higher.order.parm$theta,
                 a = as.vector(a,mode = "numeric"),
                 b = as.vector(b,mode = "numeric"),
                 minvalue = 0, maxvalue = 1)> matrix(runif(N * K), nrow = N)) * 1
      # higher order IRT model (intercept - slope style)
      # z <- outer(c(higher.order.parm$theta),a) + matrix(rep(b, N), ncol = K, byrow = T)
      # att <- ((1/(1 + exp(-z))) > matrix(runif(N * K), nrow = N)) * 1
      att.group <- matchMatrix(pattern,att)
      # att.group <- apply(att, 1, function(x) which.max(colSums(x==pattern.t)))
    }else if (tolower(att.dist) == "mvnorm"){
      if (is.null(mvnorm.parm$mean)||is.null(mvnorm.parm$sigma||is.null(mvnorm.parm$cutoffs)))
      {
        stop("multivariate normal parameters must be provided.",call. = FALSE)
      }
      atts <- MASS::mvrnorm(N,mu = mvnorm.parm$mean, Sigma=mvnorm.parm$sigma, empirical = FALSE)
      if (max(Q)==1){# dichotomous Q matrix
        att <- 1*(atts>matrix(mvnorm.parm$cutoffs,nrow = N,ncol = length(mvnorm.parm$cutoffs),byrow = TRUE))
        # Calculate which latent group each examinee belongs to return a vector
        # of N elements ranging from 1 to 2^K
        att.group <- apply(att, 1, function(x) which.max(colSums(x==pattern.t)))
      }else{
        # if Q matrix is polytomous, cutoffs must be a list - each for one attribute
        stop("multivariate normal distribution is only available for dichotomous attributes.",call. = FALSE)
      }
    }else if(tolower(att.dist) == "multinomial"){
      att.group <- sample.int(L, N, replace = T, prob = att.prior)
    }

  } else
  {
    # users specified attributes
    att <- attribute
    # Calculate which latent group each examinee belongs to return a vector
    # of N elements ranging from 1 to 2^K
    att.group <- matchMatrix(pattern,att)
  }


if(!sequential){
  LC.Prob <- uP(as.matrix(par.loc), as.matrix(catprob.matrix))  #J x L
  LC.Prob <- t(LC.Prob)  #L x J
  if(any(LC.Prob>1)||any(LC.Prob<0)) stop("Some success probabilities are greater than 1 or less than 0.",call. = FALSE)
  # Probability for each examinee on each item
  att.Prob <- LC.Prob[att.group, ]  #N x J
   Y <- 1 * (att.Prob > matrix(runif(N * J), N, J))
}else{
  LC.Prob <- sequP(as.matrix(par.loc),as.matrix(catprob.matrix),C)
  # c.LC.Prob <- LC.Prob$cPr
  u.LC.Prob <- LC.Prob$uPr
  LC.Prob <- t(u.LC.Prob) #L x S0
  if(any(LC.Prob>1)||any(LC.Prob<0)) stop("Some success probabilities are greater than 1 or less than 0.",call. = FALSE)
  att.Prob <- LC.Prob[att.group, ]  #N x J
  C1 <- C + 1
  C0 <- NULL
  ###C0:Item number for expended Qmatrix where 0 category is included
  for(j in 1:S)  C0 <- c(C0,rep(j,C1[j]))
  Y <- matrix(0,N,S)
  for (j in 1:S){ #for each item
    #Pj <- P(0),P(1),...,P(h)
    Pj <- att.Prob[,which(C0==j)]
    Y[,j] <- apply(Pj,1,function(x){sample(c(0:C[j]),1,prob=x)})
  }
  LC.Prob <- LC.Prob[,-which(!duplicated(C0))]
}

  names(delta.parm) <- names(catprob.parm) <- item.names
  output <- list(dat = Y, Q = Q, attribute = pattern[att.group, ], att.group = att.group,
                 catprob.parm = catprob.parm, delta.parm = delta.parm,
                 LCprob.parm = LC.Prob,higher.order.parm = higher.order.parm,
       mvnorm.parm = mvnorm.parm,call=simGDINAcall)
  class(output) <- "simGDINA"
  invisible(output)
}
