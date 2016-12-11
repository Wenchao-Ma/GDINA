#' Simulate responses based on the G-DINA model, DINA, DINO, ACDM, LLM or RRUM in an easy way
#'
#'
#'    This function can be used to simulate responses based on the G-DINA model
#'    or any of its subsumed models, including DINA, DINO, ACDM,
#'    LLM and R-RUM. Users need to specify the Q-matrix and item parameters (i.e., guessing and slip).
#'    Attributes can be simulated from uniform, higher order or multivariate normal
#'    distributions or supplied by users. See \code{Examples} and \code{\link{GDINA.sim}}.
#'
#'
#' @param N Sample size.
#' @param Q A required \eqn{J \times K} Q-matrix, wher J represents test
#'    length and K represents the number of attributes. For binary attributes,
#'    1 denotes attributes are measured by items and 0 means attributes are not
#'    necessary. For polytomous attributes, non-zero elements indicate which level
#'    of attributes are needed.
#' @param gs A required item parameter matrix or data frame. It
#'    must be dimension of \eqn{J \times 2}, where the first and second columns represent the guessing (or \eqn{P(0)}) and
#'    slip (or \eqn{1-P(1)}) parameters for all items, respectively.
#' @param type The type for simulating delta parameters for additive type CDMs and the G-DINA model.
#'    It can be either \code{'random'} or \code{'equal'}. The former means that the delta parameters are simulated randomly,
#'    while the latter means that each required attribute contributes equally to the probability of success (P), logit(P) or
#'    log(P) for ACDM, LLM and RRUM, respectively.
#' @param model A required vector for each item or a scalar which will be used for all
#'    items to specify which model is fitted to each item. The possible options
#'    include \code{'GDINA'},\code{'DINA'},\code{'DINO'},\code{'ACDM'},\code{'LLM'}, and \code{'RRUM'}.
#'    If \code{model} is a scalar, the specified model is fitted to all items. If \code{model} is a
#'    vector, it must have the same length as the test, where different
#'    models can be assigned to different items.
#'    It is also possible to specify models using numbers. Particularly, 0,1,2,3,4 and 5 represents
#'    \code{'GDINA'},\code{'DINA'},\code{'DINO'},\code{'ACDM'},\code{'LLM'}, and \code{'RRUM'}, respectively.
#'    The default is to fit the G-DINA model to all items.
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
#' @param digits How many decimal places in each number? The default is 4.
#' @return a list with components:
#' \describe{
#' \item{dat}{A \eqn{N \times J} data matrix of \eqn{N} examinees to \eqn{J} items}
#' \item{att}{A \eqn{N \times K} inviduals' attribute patterns}
#' \item{itemprob.param}{a list of item success probabilities for data generation}
#' \item{delta.param}{a list of delta parameters for data generation}
#' \item{higher.order.par}{Higher-order parameters}
#' \item{mvnorm.par}{multivariate normal distribution parameters}
#' \item{LC.prob}{probability of success of each item for each latent class}
#' }
#'
#' @author Wenchao Ma, Jimmy de la Torre
#'
#' @export
#'
#' @examples
#'
#'####################################################
#'#                   Example 1.1                    #
#'#            Data simulation (ACDM)                #
#'####################################################
#' set.seed(12345)
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
#' # By default, individuals are simulated from uniform distribution
#' # and deltas are simulated randomly
#' sim <- GDINA.sim.ez(N,Q,gs = item.param,model="ACDM")
#'
#' # True item probability success parameters
#' sim$itemprob.param
#'
#' # True delta parameters
#' sim$delta.param
#'
#'
#'####################################################
#'#                   Example 1.2                    #
#'#             Data simulation (LLM)                #
#'####################################################
#' set.seed(12345)
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
#' # Individuals are simulated from higher order distribution
#' # and deltas are simulated randomly
#' theta <- rnorm(N)
#' K <- ncol(Q)
#' lambda <- data.frame(a=runif(K,0.7,1.3),b=rnorm(K))
#' sim <- GDINA.sim.ez(N,Q,gs = item.param,model="LLM",att.dist = "higher.order",
#'                  higher.order.par = list(theta = theta,lambda = lambda))
#'
#' # True item probability success parameters
#' sim$itemprob.param
#'
#' # True delta parameters
#' sim$delta.param
#'
#'
#'####################################################
#'#                   Example 1.3                    #
#'#            Data simulation (RRUM)                #
#'####################################################
#' set.seed(12345)
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
#' sim <- GDINA.sim.ez(N,Q,gs = item.param,model="RRUM")
#'
#' # True item probability success parameters
#' sim$itemprob.param
#'
#' # True delta parameters
#' sim$delta.param
#'
#'
#'####################################################
#'#                   Example 1.4                    #
#'#            Data simulation (Different CDMs)      #
#'####################################################
#' set.seed(12345)
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
#' models <- c("DINA","DINO","GDINA","ACDM","LLM","RRUM","DINA","DINO","GDINA","LLM")
#' sim <- GDINA.sim.ez(N,Q,gs = item.param,model=models)
#'
#' # True item probability success parameters
#' sim$itemprob.param
#'
#' # True delta parameters
#' sim$delta.param
#'
#'
#'
GDINA.sim.ez <- function(N, Q, gs,model = "GDINA", type = "random",
                      attribute = NULL, att.dist = "uniform",
                      higher.order.par=list(theta = NULL, lambda = NULL),
                      mvnorm.par=list(mean = rep(0,ncol(Q)),
                                      sigma = diag(rep(1,ncol(Q))),cutoffs = rep(0,ncol(Q))),
                      digits=4){
  J <- nrow(Q)
  K <- ncol(Q)
  Kj <- apply(Q,1,function(x){sum(x>0)})  # The number of attributes for each item

  pattern <- alpha(K, T, Q)
  model <- model.transform(model,J)
  L <- nrow(pattern)  # the number of latent groups
  #if(length(model)==1) M <- rep(model,J)
  delta <- prob <- vector("list",J)
  for (j in 1:J){
    if (model[j]==1|model[j]==2){
      delta[[j]] <- c(gs[j,1],1 - gs[j,2]-gs[j,1])
    }else if (model[j]==3){
      p0 <- gs[j,1]
      p1 <- 1-gs[j,2]
      if(type=="equal"){
        d <- rep((p1-p0)/Kj[j],Kj[j])
        }else if(type=="random"){
        sumd <- p1-p0
        if(Kj[j]==1){d <- sumd}else{
        d <- rep(0,Kj[j])
        for(k in 1:(Kj[j]-1)){
          d[k] <- runif(1,0,sumd)
          sumd <- sumd - d[k]

        }
        d[Kj[j]] <- p1-p0-sum(d)
        }
      }
      delta[[j]] <- c(p0,d)

    }else if (model[j]==4){
      p0 <- plogis(gs[j,1])
      p1 <- plogis(1-gs[j,2])
      if(type=="equal"){
        d <- rep((p1-p0)/Kj[j],Kj[j])
      }else if(type=="random"){
        sumd <- p1-p0
        if(Kj[j]==1){d <- sumd}else{
          d <- rep(0,Kj[j])
          for(k in 1:(Kj[j]-1)){
            d[k] <- runif(1,0,sumd)
            sumd <- sumd - d[k]

          }
          d[Kj[j]] <- p1-p0-sum(d)
        }
      }
      delta[[j]] <- c(p0,d)
    }else if (model[j]==5){
      p0 <- log(gs[j,1])
      p1 <- log(1-gs[j,2])
      if(type=="equal"){
        d <- rep((p1-p0)/Kj[j],Kj[j])
      }else if(type=="random"){
        sumd <- p1-p0
        if(Kj[j]==1){d <- sumd}else{
          d <- rep(0,Kj[j])
          for(k in 1:(Kj[j]-1)){
            d[k] <- runif(1,0,sumd)
            sumd <- sumd - d[k]
            #print(sumd)
          }
          d[Kj[j]] <- p1-p0-sum(d)
        }
      }
      delta[[j]] <- c(p0,d)
    }else if (model[j]==0){
      p0 <- gs[j,1]
      p1 <- 1-gs[j,2]
       ps <- runif(2^Kj[j]-2,p0,p1)
       delta[[j]] <- c(solve(designM_GDINA(Kj[j]))%*%c(p0,ps,p1))

    }

  }
  delta <- format_delta(delta,model,Kj)


  item.param <- matrix(-1,J,2^max(Kj))
  for (j in 1:J){
    if (model[j]==0){
      Mj <- designM_GDINA(Kj[j])
    }else{
      Mj <- designM(alpha(Kj[j]),model[j])
    }
    if (model[j]<4){
      tmp <- c(Mj%*%delta[[j]])
    }else if (model[j]==4){
      tmp <- qlogis(c(Mj%*%delta[[j]]))
    }else if (model[j]==5){
      tmp <- exp(c(Mj%*%delta[[j]]))
    }

    prob[[j]] <- item.param[j,1:length(tmp)] <- round(tmp,digits)
    names(prob[[j]]) <- paste("P(",apply(alpha(Kj[j]),1,function(x){paste(x,collapse = "")}),")",sep = "")
  }

  sim <- GDINA.sim(N, Q, item.param, param.type = "prob",
                        param.format="matrix",
                        attribute = attribute, att.dist = att.dist,
                        higher.order.par=list(theta = higher.order.par$theta,
                                              lambda = higher.order.par$lambda),
                        mvnorm.par=list(mean = mvnorm.par$mean,
                                        sigma = mvnorm.par$sigma,
                                        cutoffs = mvnorm.par$cutoffs))
  names(delta) <- names(prob) <- paste("Item",1:J,sep = " ")
  return(list(itemprob.param=prob,delta.param = delta,dat = sim$dat, Q = sim$Q, att = sim$att,
              att.group = sim$att.group, LC.prob=sim$LC.Prob,higher.order.par = sim$higher.order.par,
              mvnorm.par = sim$mvnorm.par))

}
