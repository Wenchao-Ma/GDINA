#' Differential item functioning
#'
#' This function aims to detect differential item functioning based on the models estimated
#' in the \code{\link{GDINA}} function using the Wald test (Hou, de la Torre, & Nandakumar, 2014)
#'
#' @param dat A required \eqn{N \times J} \code{matrix} or \code{data.frame} consisting of binary
#' responses of \eqn{N} examinees to \eqn{J} items. Missing values need to be coded as \code{NA}.
#' @param Q A required \eqn{J \times K} item by attribute association matrix (Q-matrix; Tatsuoka, 1983),
#'    where \eqn{J} represents test
#'    length and \eqn{K} represents the number of attributes. For binary attributes,
#'    1 denotes attributes are measured by items and 0 means attributes are not
#'    necessary. For polytomous attributes, non-zero elements indicate which level
#'    of attributes are needed.
#' @param group a vector indicating the group each individual belongs to. Its length must be equal to
#'    the number of individuals. Only two groups can be dealt with for now.
#' @param model A required vector for each item or a scalar which will be used for all
#'    items to specify which model is fitted to each item. The possible options
#'    include \code{'GDINA'},\code{'DINA'},\code{'DINO'},\code{'ACDM'},\code{'LLM'}, and \code{'RRUM'}.
#'    If \code{model} is a scalar, the specified model is fitted to all items. If \code{model} is a
#'    vector, it must have the same length as the test, where different
#'    models can be assigned to different items.
#'    It is also possible to specify models using numbers. Particularly, 0,1,2,3,4 and 5 represents
#'    \code{'GDINA'},\code{'DINA'},\code{'DINO'},\code{'ACDM'},\code{'LLM'}, and \code{'RRUM'}, respectively.
#'    The default is to fit the G-DINA model to all items.
#' @param higher.order logical; \code{TRUE} indicates a higher order structure of attributes
#'    is assumed and higher order parameters will be estimated. Higher order model
#'    is specified in \code{higher.order.model}. The default is \code{FALSE}.
#' @param higher.order.model An IRT model for higher order attribute structure; It can be
#'    \code{"2PL"}, \code{"1PL"} or \code{"Rasch"}, representing two parameter logistic IRT model,
#'    one parameter logistic IRT model and Rasch model,
#'    respectively. Please note that in \code{"1PL"} model, a common slope parameter will be
#'    estimated (see \code{Details}). \code{"Rasch"} is the default.
#' @param higher.order.method The algorithm for estimation the higher order parameters; it can be either
#'    \code{"MMLE"} using marginal maximum likelihood estimation, or \code{"BL"} based on the Bock and
#'    Lieberman approach. \code{"BL"} is suitable when the number of attributes is few. It is not
#'    sensitive to sample size but can be very slow if the number of attributes is large.
#'    \code{"MMLE"}, which is the default, is suitable for most conditions but might be slow if
#'    sample size is extremely large.
#' @param higher.order.SE logical; whether the standard errors of higher order parameters are estimated?
#'    For now, the numerical method is adopted, and therefore it can be very slow. The default is \code{FALSE}.
#' @param mono.constraint logical; \code{TRUE} indicates that \eqn{P(\bm{\alpha}_1) <=P(\bm{\alpha}_2)} if
#'    for all \eqn{k}, \eqn{\alpha_{1k} <= \alpha_{2k}}. It can be a vector for each item or a scalar which will be used for all
#'    items to specify whether monotonicity constraint should be added for each item.
#' @param items The items that DIF detection will be implemented to. If \code{NULL}, all items will be conducted.
#' @param itemprob.matrix  initial item parameters; It must be a matrix giving probability of success of each reduced latent classes
#'    for each item.
#' @param SE.type Specifying different algorithms for standard errors calculation.
#' @param empirical Logical; whether empirical bayes is adopted or not? \code{TRUE} is
#'    the default when higher order attribute structure is not assumed. If estimating
#'    higher order structure, it will be \code{FALSE}.
#' @param att.prior attribute prior distribution for \eqn{2^K} latent classes. Only available for dichotomous attributes.
#'    It can be used to specify the hierarchical structure of attributes.
#'    Its length must be equal to \eqn{2^K}. Element 0 specifies which latent class does not exist.
#'    The sum of all elements does not have to be equal to 1; however, it will be standardized so that the sum is equal to 1
#'    before model calibration. When any latent class is given prior 0, standard errors for item parameters are not available.
#'    (1) If \code{empirical=F} and \code{higher.order=F}, the attribute prior distribution is fixed during model
#'    calibration; (2) if \code{empirical=T} and \code{higher.order=F}, the distribution for all latent classes
#'    with non-zero priors is updated using the empirical bayes method;
#'    (3) if \code{empirical=F} and \code{higher.order=T}, the distribution for all latent classes
#'    with non-zero priors is updated using the higher order model.
#'    The label for each latent class can be obtained by calling \code{alpha(K)}. See \code{examples} for more info.
#' @param att.str logical; whether attributes have any structure?
#' @param nstarts how many sets of starting values? The default is 3.
#' @param conv.crit The convergence criterion for max absolute change in item parameters.
#' @param maxitr The maximum iterations of EM cycles allowed.
#' @param higher.order.param A matrix or data frame providing higher order parameters. If supplied, it must be of dimension \eqn{K\times 2}.
#'    The first column is the slope parameters and the second column is the intercept.
#' @param digits How many decimal places in each number? The default is 4.
#'
#'  @return a data frame giving the Wald statistics and associated p-values.
#'
#' @seealso \code{\link{GDINA}}
#'
#' @examples
#' set.seed(12345)
#' N <- 1000
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
#' sim1 <- GDINA.sim.ez(N,Q,gs = item.param,model="ACDM",type="equal")
#' sim2 <- GDINA.sim.ez(N,Q,gs = item.param,model=c(rep("ACDM",9),"DINA"),type="equal")
#' dat <- rbind(sim1$dat,sim2$dat)
#' gr <- c(rep(1,N),rep(2,N))
#' dif.wald(dat,Q,group=gr)
#'
#'
#' @references
#' Hou, L., de la Torre, J., & Nandakumar, R. (2014). Differential item functioning assessment in cognitive diagnostic modeling: Application of the Wald test to
#' investigate DIF in the DINA model. \emph{Journal of Educational Measurement, 51}, 98-125.
#'


dif.wald <- function(dat, Q, group, model = "GDINA", items=NULL, higher.order = FALSE, higher.order.model =
                        "Rasch", higher.order.method = "MMLE", higher.order.SE = FALSE,
                      itemprob.matrix = NULL, higher.order.param = NULL,
                      mono.constraint = FALSE, SE.type = 1,
                      empirical = !higher.order, att.prior = NULL, att.str = FALSE,
                      nstarts = 3, conv.crit = 0.001, maxitr = 1000,digits = 4){
  if (nrow(dat)!=length(group))stop("The length of group variable must be equal to the number of individuals.",call. = FALSE)
  if (length(unique(group))!=2)stop("Only two group DIF can be examined.",call. = FALSE)
  out <- vector("list",length(unique(group)))
  for (g in 1:length(unique(group))){
    out[[g]] <- GDINA(dat[group==unique(group)[g],], Q, model, higher.order = higher.order, higher.order.model =
            higher.order.model, higher.order.method = higher.order.method,
          higher.order.SE = higher.order.SE, verbose = FALSE,
          itemprob.matrix = itemprob.matrix, higher.order.param = higher.order.param,
          mono.constraint = mono.constraint, person.est = FALSE, SE = TRUE, SE.type = SE.type,
          empirical = empirical, att.prior = att.prior, att.str = att.str,
          nstarts = nstarts, conv.crit = conv.crit, maxitr = maxitr)
  }
  J <- nrow(Q)
  if (is.null(items)) items <- c(1:J)
  Wp <- matrix(0,length(items),2)
  for (j in items){
    x <- c(out[[1]]$itemprob.param[[j]],out[[2]]$itemprob.param[[j]])
    npar <- length(out[[1]]$itemprob.param[[j]])
    m0 <- matrix(0,npar,npar)
    R <- cbind(diag(npar),-1*diag(npar))
    itmloc <- which(out[[1]]$technicals$itemprob.covIndex[,2]==j,arr.ind = T)
    vcov <- cbind(rbind(out[[1]]$technicals$itemprob.cov[itmloc,itmloc],m0),
                  rbind(m0,out[[2]]$technicals$itemprob.cov[itmloc,itmloc]))
    Wp[j,1] <- t(R%*%x)%*%solve(R%*%vcov%*%t(R))%*%(R%*%x)
    Wp[j,2] <- pchisq(Wp[j,1],nrow(R),lower.tail = FALSE)
  }
Wp <- round(data.frame(Wp),digits)
colnames(Wp) <- c("Wald stat.","p-value")
return(Wp)

}
