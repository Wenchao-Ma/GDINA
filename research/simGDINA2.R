#' Simulate responses based on the G-DINA model, DINA, DINO, ACDM, LLM or RRUM
#'
#'
#'    This function can be used to simulate responses based on the G-DINA model
#'    or any of its subsumed models, including DINA, DINO, ACDM,
#'    LLM and R-RUM. Users need to specify the Q-matrix and item parameters.
#'    Attributes can be simulated from uniform, higher order or multivariate normal
#'    distributions or supplied by users. See \code{Examples} and \code{\link{GDINA}}.
#'
#'
#' @param N Sample size.
#' @param Q A required \eqn{J \times K} Q-matrix, wher J represents test
#'    length and K represents the number of attributes. For binary attributes,
#'    1 denotes attributes are measured by items and 0 means attributes are not
#'    necessary. For polytomous attributes, non-zero elements indicate which level
#'    of attributes are needed.
#' @param item.param A required item parameter matrix or list. If \code{param.type} is
#'    \code{'prob'}, \code{item.param} is user-supplied probabilities of success of each reduced latent
#'    class on each item.
#' @param param.type parameter type; It can be (1) \code{"prob"}, which represents the probabilities of success for each item -
#'    it should have the same format as \code{itemprob.matrix} or \code{itemprob.param} of the object from class \code{GDINA}
#'    if \code{param.format} is \code{"matrix"} or \code{"list"}, respectively.
#'    (2) \code{"gs"}, which represents guessing and slip parameters for DINA and DINO model only. It
#'    must have dimension of \eqn{J \times 2} where the first and second column represent guessing and
#'    slip parameters for all items, respectively; or (3) \code{"delta"}, which represents delta parameters
#'    in the GDINA model (de la Torre, 2011). It must have the same format as \code{delta.param} of the
#'    object from the class \code{GDINA}, which is a list of \eqn{J} elements.
#' @param param.format item parameter format; Only applicable when \code{param.type=prob}. It can be
#'    \code{"matrix"} or \code{"list"}.
#' @param model Only applicable if \code{param.type=gs} or \code{param.type=delta};
#'    A vector for each item or a scalar which will be used for all
#'    items to specify which model is simulated for each item. The possible options
#'    include \code{'GDINA'},\code{'DINA'},\code{'DINO'}, or \code{'ACDM'}.
#'    If \code{model} is a scalar, the specified model is used for all items. If \code{model} is a
#'    vector, it must have the same length as the test, where different
#'    models can be assigned to different items.
#' @param linkfunc link function used for data simulation. It is only applicable if \code{param.type=delta}.
#'    Possible options include \code{'identity'}, \code{'log'}, and \code{'logit'}. Note that
#'    RRUM is obtained if \code{model='ACDM'} and \code{linkfunc='log'} and LLM can be obtained
#'    if \code{model='ACDM'} and \code{linkfunc='logit'}.
#' @param attribute person attributes. If this is not supplied, it is simulated
#'    from a distribution specified in \code{att.dist}.
#' @param att.dist the attribute distribution. It can be \code{"uniform"}, \code{"higher.order"} or
#'    \code{"mvnorm"} for uniform, higher order and multivariate normal distribution, respectively.
#'    The default is the uniform distribution.
#' @param higher.order.par A list specifying parameters for higher order distribution for attributes
#'    if in \code{att.dist=higher.order}. Particularly, \code{theta} is a
#'    vector of length \eqn{N} representing the higher order ability
#'    for each examinee. and \code{lambda} is a \eqn{K \times 2} matrix. Column 1 gives slopes of higher-order
#'    model and column 2 gives the intercepts. See \code{\link{GDINA}} for the formulations of the higher-order
#'    models.
#' @param mvnorm.par a list of parameters for multivariate normal attribute distribution. \code{mean} is a vector of length \eqn{K}
#'    specifying the mean of multivariate normal distribution; and \code{sigma} is a positive-definite
#'    symmetric matrix specifying the variance-covariance matrix. \code{cutoffs} is a vector giving the
#'    cutoff for each attribute. See \code{Examples}.
#' @return a list with components:
#' \describe{
#' \item{dat}{A \eqn{N \times J} data matrix of \eqn{N} examinees to \eqn{J} items}
#' \item{att}{A \eqn{N \times K} inviduals' attribute patterns}
#' \item{itemprob.param}{item parameter matrix}
#' \item{higher.order.par}{Higher-order parameters}
#' \item{mvnorm.par}{multivariate normal distribution parameters}
#' \item{LC.prob}{Probability of success of each item for each latent class}
#' }
#'
#' @author Wenchao Ma, Jimmy de la Torre
#'
#' @export
#'
#' @references
#'
#' Chiu, C.-Y., Douglas, J. A., & Li, X. (2009). Cluster analysis for cognitive diagnosis: Theory and applications. \emph{Psychometrika, 74}, 633-665.
#'
#' de la Torre, J. (2011). The generalized DINA model framework. \emph{Psychometrika, 76}, 179-199.
#'
#' @examples
#'
#'####################################################
#'#                   Example 1.1                    #
#'#       Data simulation (DINA/DINO)                #
#'#  using probability of success in matrix format   #
#'####################################################
#' # Although there are multiple ways of simulating
#' # DINA/DINO data using this function, this is
#' # probably the easiest way
#' # -- guessing and slip parameters for each item
#' # need to be specified in a matrix or data frame
#' # of dimension J x 2 (column 1 is guessing and
#' # column 2 is slip)
#' # e.g.,
#' N <- 500
#' Q <- sim10GDINA$simQ
#' item.param <- matrix(c(0.1,0.2,
#'                        0.1,0.2,
#'                        0.1,0.2,
#'                        0.1,0.2,
#'                        0.1,0.2,
#'                        0.1,0.2,
#'                        0.1,0.2,
#'                        0.1,0.2,
#'                        0.1,0.2,
#'                        0.1,0.2),ncol = 2, byrow = TRUE)
#' # Item 1-5 are DINA; item 6-10 are DINO
#' model <- c(rep(c("DINA","DINO"),each=5))
#' # Simulated DINA and DINO model
#' simD <- GDINA.sim(N,Q,item.param = item.param,
#'                   param.type = "gs", param.format = "matrix", model = model)
#'
#' # check item parameters
#' simD$item.param
#'
#'####################################################
#'#                   Example 1.2                    #
#'#          Data simulation (all CDMs)              #
#'#  using probability of success in matrix format   #
#'####################################################
#' #     NOTE
#' # Any CDMs can be simulated in this way
#' # because the model does not matter if
#' # item parameter is success probability
#' # For example, different CDMs can be obtained
#' # using different probability of success
#' # e.g.     P(00) P(10) P(01) P(11)
#' # DINA      0.2   0.2   0.2   0.9
#' # DINO      0.2   0.8   0.8   0.8
#' # ACDM      0.2   0.5   0.5   0.8
#' # Therefore, specifying CDMs is NOT necessary.
#'
#' set.seed(12345)
#'
#' N <- 500
#' Q <- sim10GDINA$simQ
#' itempar.matrix <- sim10GDINA$simItempar
#' # the format of item param is as follows:
#' # rows represent items
#' # if item j requires Kj attributes, 2^Kj success probabilities
#' # need to be specified in the first 2^Kj cells
#' # e.g., item 1 only requires 1 attribute
#' # therefore P(0) and P(1) should be specified in the first
#' # and the second cells in the first row;
#' # similarly, item 10 requires 3 attributes,
#' # P(000),P(100),P(010)...,P(111) should be specified in
#' # the first 8 cells
#' # the latent class represented by each cell can be obtained
#' # by calling alpha(Kj)
#' # DINA and DINO can also be simulated in this way
#' # for example, item 5 is a DINA item and item 6 is a DINO item
#'
#' itempar.matrix
#' # 0.2	 0.9	  NA	NA	NA	NA	NA	NA
#' # 0.1	 0.8	  NA	NA	NA	NA	NA	NA
#' # 0.1	 0.9	  NA	NA	NA	NA	NA	NA
#' # 0.1	 0.3	 0.5	0.9	NA	NA	NA	NA
#' # 0.1 	 0.1	 0.1	0.8	NA	NA	NA	NA
#' # 0.2	 0.9	 0.9	0.9	NA	NA	NA	NA
#' # 0.1	0.45	0.45	0.8	NA	NA	NA	NA
#' # 0.1	0.28	0.28	0.8	NA	NA	NA	NA
#' # 0.1	 0.4	 0.4	0.8	NA	NA	NA	NA
#' # 0.1	 0.2	 0.3	0.4	0.4	0.5	0.7	0.9
#'
#'
#' sim <- GDINA.sim(N,Q,item.param = itempar.matrix, param.type = "prob", param.format = "matrix")
#'
#'####################################################
#'#                   Example 1.3                    #
#'#          Data simulation (all CDMs)              #
#'#  using probability of success in list format     #
#'####################################################
#'
#'# success probabilities for each item can also be provided in list format as follows:
#'
#'itempar.list <- list(item1=c(0.2,0.9),
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
#'sim2 <- GDINA.sim(N,Q,item.param = itempar.list, param.type = "prob", param.format = "list")
#'
#'# check results - identical
#'all(sim$dat==sim2$dat)
#'
#'####################################################
#'#                   Example 2                      #
#'#            Data simulation (GDINA)               #
#'#      using delta parameters in list format       #
#'####################################################
#'
#'# NOTE:
#'# CDMs matter if item parameter is delta
#'# There can be multiple parameterizations if considering different link functions
#'  delta.list <- list(c(0.2,0.7),
#'                     c(0.1,0.7),
#'                     c(0.1,0.8),
#'                     c(0.1,0.2,0.4,0.2),
#'                     c(0.1,0.0,0.0,0.7),
#'                     c(0.2,0.7,0.7,-0.7),
#'                     c(0.1,0.35,0.35,0.0),
#'                     c(0.1,0.18,0.18,0.34),
#'                     c(0.1,0.3,0.3,0.1),
#'                     c(0.1,0.1,0.2,0.3,0.0,0.0,0.1,0.1))
#'set.seed(12345)
#' N <- 500
#' Q <- sim10GDINA$simQ
#'sim3 <- GDINA.sim(N,Q,item.param = delta.list,
#'                  param.type = "delta", param.format = "list", model = "GDINA")
#'all(sim3$dat==sim$dat)
#'
#'
#'####################################################
#'#                   Example 3                      #
#'#            Data simulation (all CDMs)            #
#'#      using delta parameters in list format       #
#'####################################################
#'
#' delta2.list <- list(c(0.2,0.7),
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
#' model <- c("GDINA","GDINA","GDINA","DINA","DINO","GDINA","ACDM","ACDM","ACDM","GDINA")
#' linkf <- c("identity","identity","identity","identity","identity","identity",
#'            "identity","logit","log","identity")
#' N <- 500
#' Q <- sim10GDINA$simQ
#' sim5 <- GDINA.sim(N,Q,item.param = delta2.list,
#'                   param.type = "delta", param.format = "list", model = model,linkfunc = linkf)
#'
#'####################################################
#'#                   Example 4                      #
#'#      Data simulation (log logit GDINA)           #
#'#      using delta parameters in list format       #
#'####################################################
#'
#'# define logit function
#'logit <- function(p) log(p/(1-p))
#'
#'Q3 <- matrix(c(1,0,0,
#'               0,1,0,
#'               0,0,1,
#'               1,1,0,
#'               1,0,1,
#'               0,1,1),nrow=6,byrow=TRUE)
#'
#'# Note that when Kj*=2, f(P00)=delta0; f(P10)=delta0+delta1;
#'# f(P01)=delta0+delta2; f(P11)=delta0+delta1+delta2+delta12
#'# for GDINA model, where f can be log or logit
#'# in the following delta3.list, we need to specify
#'# these deltas
#'# Item 1,2,4 are logit link GDINA model
#'# item 3 and 5 are log like GDINA model
#'# item 6 is identity link GDINA model
#'# the data is simulated so that P(0)=P(00)=0.2
#'# P(10)=0.4; P(01)=0.5 and P(11)=0.8
#'# This can be verified by estimating the data using the GDINA function
#'delta3.list <- list(c(logit(0.2),logit(0.8)-logit(0.2)),
#'                     c(logit(0.2),logit(0.8)-logit(0.2)),
#'                     c(log(0.2),log(0.8)-log(0.2)),
#'                     c(logit(0.2),logit(0.4)-logit(0.2),logit(0.5)-logit(0.2),
#'                     logit(0.8)-logit(0.4)-logit(0.5)+logit(0.2)),
#'                     c(log(0.2),log(0.4)-log(0.2),log(0.5)-log(0.2),
#'                     log(0.8)-log(0.4)-log(0.5)+log(0.2)),
#'                     c(0.2,0.2,0.3,0.1))
#'
#' model <- "GDINA" #all items are GDINA models
#' linkf <- c("logit","logit","log","logit","log","identity")
#' N <- 20000 # large sample is needed for a short test
#' sim6 <- GDINA.sim(N,Q3,item.param = delta3.list,
#'                   param.type = "delta", param.format = "list", model = model,linkfunc = linkf)
#' mod6 <- GDINA(sim6$dat,Q3)
#' mod6$item.prob
#'
#'####################################################
#'#                   Example 5                      #
#'#      Data simulation (higher order model)        #
#'####################################################
#'
#' Q <- sim10GDINA$simQ
#' itempar.matrix <- sim10GDINA$simItempar
#'
#' set.seed(12345)
#' theta <- rnorm(N)
#' K <- ncol(Q)
#' lambda <- data.frame(a=runif(K,0.7,1.3),b=rnorm(K))
#' simHO <- GDINA.sim(N,Q,item.param = itempar.matrix, param.type = "prob", param.format = "matrix",
#'                  att.dist = "higher.order",
#'                  higher.order.par = list(theta = theta,lambda = lambda))
#'
#'####################################################
#'#                   Example 6                      #
#'#      Data simulation (higher order model)        #
#'#  using the multivariate normal threshold model   #
#'####################################################

#'
#' # See Chiu et al., (2009)
#'
#' N <- 500
#' Q <- sim10GDINA$simQ
#' cutoffs <- qnorm(c(1:K)/(K+1))
#' m <- rep(0,K)
#' vcov <- matrix(0.5,K,K)
#' diag(vcov) <- 1
#' simMV <- GDINA.sim(N,Q,item.param = itempar.matrix, param.type = "prob", param.format = "matrix",
#'                  att.dist = "mvnorm",
#'                  mvnorm.par=list(mean = m, sigma = vcov,cutoffs = cutoffs))
#'
#'

GDINA.sim <- function(N, Q, item.param, param.type = "prob",
                      param.format="matrix",model = "GDINA", linkfunc = "identity",
                      attribute = NULL, att.dist = "uniform",
                      higher.order.par=list(theta = NULL, lambda = NULL),
                      mvnorm.par=list(mean = rep(0,ncol(Q)),sigma = diag(rep(1,ncol(Q))),cutoffs = rep(0,ncol(Q))))
{


  item.param_copy <- item.param
  J <- nrow(Q)
  K <- ncol(Q)
  pattern <- alpha(K, T, Q)
  L <- nrow(pattern)  # the number of latent groups

  model <- model.transform(model,J)
  par.loc <- eta.loc(Q)
  ##### SOME INPUT CHECKS ARE NEEDED#########
  if (param.type=="gs"){ # DINA or DINO
    if(!all(model%in%c(1,2))) stop ("Only DINA or DINO is applicable if param.type is 'gs'.",call. = FALSE)
    if (nrow(item.param)!=J|ncol(item.param)!=2) stop ("item.param must have J rows and 2 columns if param.type is 'sg'.",call. = FALSE)
    Kj <- apply(Q,1,function(x){sum(x>0)})  # The number of attributes for each item
    Kjmax <- max(Kj) # the maximum attributes required for each item
    item.param <- matrix(-1,J,2^Kjmax)
    for (j in 1:J){
      if(model[j]==1){#DINA
        item.param[j,1:(2^Kj[j]-1)] <- item.param_copy[j,1] # guessing
        item.param[j,2^Kj[j]] <- 1-item.param_copy[j,2] # 1-slip
      }else if(model[j]==2){#DINO
        item.param[j,1] <- item.param_copy[j,1] # guessing
        item.param[j,2:2^Kj[j]] <- 1-item.param_copy[j,2] # 1-slip
      }
    }
  }
  else if(param.type=="delta")
  {
    if (!is.list(item.param)|!length(item.param)==J) stop ("item.param must be a list of J elements if param.type is 'delta'.",call. = FALSE)
    if (!is.list(item.param)|!length(item.param)==J) stop ("item.param must be a list of J elements if param.format is 'list'.",call. = FALSE)
    if (!all(linkfunc%in%c("identity","log","logit"))) stop ("linkfunc must be identity, log or logit.",call. = FALSE)
    if (length(linkfunc)==1) linkfunc <- rep(tolower(linkfunc),J)
    Kj <- apply(Q,1,function(x){sum(x>0)})  # The number of attributes for each item
    Kjmax <- max(Kj) # the maximum attributes required for each item
    item.param <- matrix(-1,J,2^Kjmax)
    for (j in 1:J){
      if (model[j]==0){
        Mj <- designM_GDINA(Kj[j])
      }else{
        Mj <- designM(alpha(Kj[j]),model[j])
      }
      if (linkfunc[j]=="identity"){
        tmp <- c(Mj%*%item.param_copy[[j]])
      }else if (linkfunc[j]=="logit"){
        tmp <- inv.logit(c(Mj%*%item.param_copy[[j]]))
      }else if (linkfunc[j]=="log"){
        tmp <- exp(c(Mj%*%item.param_copy[[j]]))
      }

      item.param[j,1:length(tmp)] <- tmp
    }

  }
  else if(param.type=="prob"){
    if(param.format=="list"){
      Kj <- apply(Q,1,function(x){sum(x>0)})  # The number of attributes for each item
      Kjmax <- max(Kj) # the maximum attributes required for each item
      item.param <- matrix(-1,J,2^Kjmax)
      for (j in 1:J){
        item.param[j,1:length(item.param_copy[[j]])] <- item.param_copy[[j]]
      }
    }
  }
  LC.Prob <- uP(as.matrix(par.loc), as.matrix(item.param))  #J x L
  LC.Prob <- t(LC.Prob)  #L x J

#higher.order.param <- list(theta=NULL,lambda)
  # specify true person parameters user do not supply attributes
  if (is.null(attribute))
  {
    if (tolower(att.dist) == "uniform")
    {
      # uniform distribution
      att.group <- sample(1:L, N, replace = T)  #uniform distribution
    } else if (tolower(att.dist) == "higher.order")
    {
      if (max(Q) > 1)
      {
        return(warning("Higher order structure is not allowed currently when attributes are polytomous."))
      }
      if (is.null(higher.order.par$lambda))
      {
        a <- runif(K, 0.7, 1.3)
        b <- rnorm(K, 0, 0.5)
      } else
      {
        a <- higher.order.par$lambda[, 1]
        b <- higher.order.par$lambda[, 2]
      }
      if (is.null(higher.order.par$theta))
      {
        theta <- higher.order.par$theta <- rnorm(N, 0, 1)
      }
      # higher order 2PL IRT model (intercept - slope style)
      z <- matrix(rep(a, N), ncol = K, byrow = T) * matrix(rep(theta,
                                                               K), ncol = K) + matrix(rep(b, N), ncol = K, byrow = T)
      att <- ((1/(1 + exp(-z))) > matrix(runif(N * K), nrow = N)) * 1

      lambda <- cbind(a, b)
      att.group <- apply(att, 1, function(x)
      {
        which.max(rowSums(matrix(rep(x,nrow(pattern)),ncol = length(x),byrow = T)==pattern))
      })
    }else if (tolower(att.dist) == "mvnorm"){
      atts <- MASS::mvrnorm(N,mu = mvnorm.par$mean, Sigma=mvnorm.par$sigma)
      if (max(Q)==1){# dichotomous Q matrix
        att <- 1*(atts>matrix(mvnorm.par$cutoffs,nrow = N,ncol = length(mvnorm.par$cutoffs),byrow = TRUE))
        # Calculate which latent group each examinee belongs to return a vector
        # of N elements ranging from 1 to 2^K
        att.group <- apply(att, 1, function(x)
        {
          which.max(rowSums(matrix(rep(x,nrow(pattern)),ncol = length(x),byrow = T)==pattern))
        })
      }else{
        # if Q matrix is polytomous, cutoffs must be a list - each for one attribute
        stop("multivariate normal distribution is only available for dichotomous attributes.",call. = FALSE)
      }
    }

  } else
  {
    # users specified attributes
    att <- attribute
    # Calculate which latent group each examinee belongs to return a vector
    # of N elements ranging from 1 to 2^K
    att.group <- apply(att, 1, function(x)
    {
      which.max(rowSums(matrix(rep(x,nrow(pattern)),ncol = length(x),byrow = T)==pattern))
    })
  }

  # Probability for each examinee
  att.Prob <- LC.Prob[att.group, ]  #N x S0
  Y <- 1 * (att.Prob > matrix(runif(N * J), N, J))

  return(list(dat = Y, Q = Q, att = pattern[att.group, ], att.group = att.group,
              item.param = item.param, LC.prob=LC.Prob,higher.order.par = higher.order.par,
              mvnorm.par = mvnorm.par))
}
