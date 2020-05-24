#' Q-matrix validation
#'
#' Q-matrix validation for the (sequential) G-DINA model based on PVAF (de la Torre & Chiu, 2016; Najera, Sorrel, & Abad, 2019), stepwise Wald test (Ma & de la Torre, 2020) or mesa plot (de la Torre & Ma, 2016).
#' All these methods are suitable for dichotomous and ordinal response data. If too many modifications are suggested based on the default PVAF method, you are suggested to try the stepwise Wald test method or predicted cutoffs.
#' You should always check the mesa plots for further examination.
#'
#' @param GDINA.obj An estimated model object of class \code{GDINA}
#' @param eps cutoff value for PVAF from 0 to 1. Default = 0.95. Note that it can also be -1, indicating the predicted cutoff based on Najera, P., Sorrel, M., and Abad, P. (2019).
#' @param method which Q-matrix validation method is used? Can be either \code{"PVAF"} or \code{"wald"}.
#' @param wald.args a list of arguments for the stepwise Wald test method.
#' \describe{
#' \item{SE.type}{Type of covariance matrix for the Wald test}
#' \item{alpha.level}{alpha level for the wald test}
#' \item{GDI}{It can be 0, 1 or 2; 0 means GDI is not used to choose the attribute -
#' when more than one attributes are significant, the one with the largest p-value will be selected;
#' GDI=1 means the attribute with the largest GDI will be selected; GDI=2 means the q-vector with
#' the largest GDI will be selected.}
#' \item{verbose}{Print detailed information or not?}
#' \item{stepwise}{\code{TRUE} for stepwise approach and \code{FALSE} for forward approach}
#' }
#' @param digits How many decimal places in each number? The default is 4.
#' @return An object of class \code{Qval}. Elements that can be
#' extracted using \code{extract} method include:
#' \describe{
#' \item{sug.Q}{suggested Q-matrix}
#' \item{Q}{original Q-matrix}
#' \item{varsigma}{varsigma index}
#' \item{PVAF}{PVAF}
#' }
#'
#' @include GDINA.R
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @references
#' de la Torre, J. & Chiu, C-Y. (2016). A General Method of Empirical Q-matrix Validation. \emph{Psychometrika, 81}, 253-273.
#'
#' de la Torre, J., & Ma, W. (2016, August). Cognitive diagnosis modeling: A general framework approach and its implementation in R. A Short Course at the Fourth Conference on Statistical Methods in Psychometrics, Columbia University, New York.
#'
#' Ma, W., & de la Torre, J. (2020). An empirical Q-matrix validation method for the sequential G-DINA model. \emph{British Journal of Mathematical and Statistical Psychology, 73}, 142-163.
#'
#' Najera, P., Sorrel, M., & Abad, P. (2019). Reconsidering Cutoff Points in the General Method of Empirical Q-Matrix Validation. \emph{Educational and Psychological Measurement}.
#'
#' @seealso \code{\link{GDINA}}
#' @export
#' @examples
#'\dontrun{
#' ################################
#' #
#' # Binary response
#' #
#' ################################
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' Q[10,] <- c(0,1,0)
#'
#' # Fit the G-DINA model
#' mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
#'
#' # Q-validation using de la Torre and Chiu's method
#' pvaf <- Qval(mod1,method = "PVAF",eps = 0.95)
#' pvaf
#' extract(pvaf,what = "PVAF")
#' #See also:
#' extract(pvaf,what = "varsigma")
#' extract(pvaf,what = "sug.Q")
#'
#' # Draw mesa plots using the function plot
#'
#' plot(pvaf,item=10)
#'
#' #The stepwise Wald test
#' stepwise <- Qval(mod1,method = "wald")
#' stepwise
#' extract(stepwise,what = "PVAF")
#' #See also:
#' extract(stepwise,what = "varsigma")
#' extract(stepwise,what = "sug.Q")
#'
#' #Set eps = -1 to determine the cutoff empirically
#' pvaf2 <- Qval(mod1,method = "PVAF",eps = -1)
#' pvaf2
#' ################################
#' #
#' # Ordinal response
#' #
#' ################################
#' seq.est <- GDINA(sim20seqGDINA$simdat,sim20seqGDINA$simQ,sequential = TRUE)
#' stepwise <- Qval(seq.est, method = "wald", eps = -1)
#'}

Qval <- function(GDINA.obj, method = "PVAF", eps = 0.95, digits = 4, wald.args = list()){

  if (class(GDINA.obj) != "GDINA")
    stop("GDINA.obj must be a GDINA object from GDINA function.", call. = FALSE)


  if (!is.null(extract(GDINA.obj, "att.str")))
    stop("Q-matrix validation is not available if attributes are structured.",
         call. = FALSE)

  if(eps == -1){
    gs <- coef(GDINA.obj, what = "gs") # item parameters
    eps <- plogis(-0.4045140782147791 +
                       4.8404850955032684E-4*extract(GDINA.obj,"nobs") +
                       2.8667570118638275*(1-sum(colMeans(gs))) +
                       -0.003315555999671906*extract(GDINA.obj,"nitem"))
  }

  if (eps > 1 || eps < 0)
    stop("eps must be between 0 and 1.", call. = FALSE)

  if(extract(GDINA.obj,"ngroup")>1) stop("Only available for single-group models.",call. = FALSE)

  # if(any(extract(GDINA.obj,"models_numeric")<0)||any(extract(GDINA.obj,"models_numeric")>5))
  #   stop("Models must be GDINA, DINA, DINO, ACDM, LLM or RRUM",call. = FALSE)

  if (max(extract(GDINA.obj, "Q")) > 1)
    stop("Q-validation can only be used for dichotomous attribute G-DINA model.",
         call. = FALSE)



  updated.wald.args <- NULL
  if(toupper(method)=="PVAF"){
    # if(any(extract(GDINA.obj,"models")!="GDINA"))
    #   warning("Saturated G-DINA model may be used to calibrate all items for better performance.",call. = FALSE)
    ret <- Qval_PVAF(GDINA.obj,eps = eps, digits = digits)
  }else if (toupper(method)=="WALD"){
    if(any(extract(GDINA.obj,"models")!="GDINA"))
      stop("Saturated G-DINA model needs to be fitted to all items.",call. = FALSE)
    args.default <- list(GDINA.obj = GDINA.obj, SE.type = 2,
                         alpha.level = 0.05, GDI = 2, PVAF = eps,
                         verbose = FALSE, stepwise = TRUE,digits=digits)
    updated.wald.args <- modifyList(args.default,wald.args)
    ret <- do.call(Qval_wald,updated.wald.args)
  }
  ret$method <- method
  ret$wald.args <- updated.wald.args
  ret$sequential <- extract(GDINA.obj,"sequential")
  class(ret) <- "Qval"
  return(ret)
}


Qval_wald <- function(GDINA.obj, SE.type = 2,
                      alpha.level = 0.05, GDI = 2, PVAF = 0.95,
                      verbose = FALSE, stepwise = TRUE,digits=4,...){
  seqent <- extract(GDINA.obj,"sequential")

  if(seqent){
    dat <- extract(GDINA.obj,"seq.dat")
  }else{
    dat <- extract(GDINA.obj,"dat")
  }

  Qc <- extract(GDINA.obj,"Qc")

  Qr <- Q <- extract(GDINA.obj,"Q")

  N <- extract(GDINA.obj,"nobs")

  # number of categories
  J <- extract(GDINA.obj,"ncat")

  K <- extract(GDINA.obj,"natt")

  L <- 2^K

  Kj <- rowSums(attributepattern(K)[-1,])
  w <- extract(GDINA.obj,"posterior.prob") #1 x L

  scofun <- score(GDINA.obj,parm="prob") # a list with # of category elements

  # RN <- NgRg(GDINA.obj$technicals$logposterior.i,seqdat,eta.loc(matrix(1,J,K)),1-is.na(seqdat))
  expectedR <- extract(GDINA.obj,"expectedCorrect.LC")
  expectedN <- extract(GDINA.obj,"expectedTotal.LC")

  est.p <- t((expectedR + 1e-10)/(expectedN + 2*1e-10))
  # est.p[is.nan(est.p)] <- 0.5
  patt <- attributepattern(K)[-1,]
  loc <- LC2LG(patt) #2^K-1 x 2^K
  vsg <- varsigma(as.matrix(t(loc)),as.matrix(est.p),c(w)) # versigma
  vsg0 <- vsg/vsg[,ncol(vsg)] #pvaf
  vsgK <- vsg0[,1:K] # pvaf for single att
  if(verbose){
    cat("\nGDI\n")
    print(vsg[,1:K])
    cat("\nPVAF\n")
    print(vsgK)
    cat("First selected attribute based on GDI is",apply(vsgK,1,which.max.randomtie))
  }
  #first attribute is the one with the largest GDI
  first.att <- apply(vsgK,1,which.max.randomtie)
  att.monitor <- Qrr <- vector("list",J)
  inichoose <- numeric(J)
  fullset <- c(1:K)
  # iteras <- NULL
  item <- c(1:nrow(Q))
  for (j in item) {

    item.no <- Qc[j,1]

    inichoose[j] <- currentset <- first.att[j] # initial att. ---largest GDI


    att.monitor[[j]] <- c(att.monitor[[j]],currentset)
    loop <- ifelse(vsgK[j,currentset]>=PVAF,FALSE,TRUE)
    it <- 1
    while(loop&&it<K){
      difset <- setdiff(fullset,currentset)
      add.a <- NULL
      #********************************************************************Second round eval. Wald forward
      Rm <- Rmatrix.att(length(currentset)+1)
      Wp.a <- NULL
      for (k in difset){ # 2nd att.
        #
        Qr <- extract(GDINA.obj,"Q")
        Qr[j,seq_len(ncol(Qr))] <- 0
        Qr[j,c(currentset,k)] <- 1

        etas <- LC2LG(as.matrix(Qr))
        itemparj <- extract(GDINA.obj,"catprob.parm")

        itemparj[[j]] <- aggregate(expectedR[j,],by=list(etas[j,]),sum)$x/
          aggregate(expectedN[j,],by=list(etas[j,]),sum)$x
        index <- data.frame(Cat=rep(1:length(rowSums(Qr) ),2^rowSums(Qr)) )
        index$Column <- seq_len(length(index$Cat))
        sco <- scofun

        sco[[j]] <- score_pj(Xj = dat[,j],                   # a vector of item responses to item j
                             parloc.j=etas[j,,drop=FALSE],         # parameter locations for item j - H by 2^K matrix
                             catprob.j=itemparj[j],        # a list with H elements giving the reduced catprob.parm for each nonzero category
                             logpost=indlogPost(GDINA.obj))[[1]]
        if(SE.type==2){
          if(extract(GDINA.obj,"att.dist")!="saturated"){
            v <- inverse_crossprod(do.call(cbind,sco))
          }else{
            v <- inverse_crossprod(do.call(cbind,sco[-length(sco)]))
          }

        }else if(SE.type==3){
          if(extract(GDINA.obj,"att.dist")!="saturated")
            warning("structural parameters are not considered in the calculation of covariance matrix.",call. = FALSE)
          v <- inverse_crossprod(do.call(cbind,sco))
        }
        # print(j)
        # print(dim(v))
        # print(index)
        Varj <- v[index$Column[which(index$Cat==j)],index$Column[which(index$Cat==j)]]
        # print(Varj)
        Cparj <- itemparj[[j]]
        # print(Cparj)
        # p values for attribute in target attribute + current set
        Wp.a <- rbind(Wp.a,sapply(Rm,function(x) {
          # wald statistic & p values
          W <-t(x %*% Cparj) %*% MASS::ginv(x %*% Varj %*% t(x)) %*% (x %*% Cparj)
          p <- pchisq(W,nrow(x),lower.tail = FALSE)
          return (p)
        }))
      }
      #********************************************************************end 2nd round eval. att. Forward
      remove.new <- add.new <- NULL
      for (d in 1:length(difset)){
        newq <- vector("numeric",K)
        newq[c(currentset,difset[d])] <- 1
        newq.loc <- rowMatch(patt,newq)$row.no
        # p values, GDI for target att.
        add.new <- rbind(add.new,c(difset[d],Wp.a[d,1+sum(difset[d]>currentset)],vsgK[j,difset[d]],vsg0[j,newq.loc]))
        remove.new <- rbind(remove.new,Wp.a[d,-(1+sum(difset[d]>currentset)),drop=FALSE]) #p values for att in current set
      }
      # iteras <- rbind(iteras,cbind(j,it,add.new))
      if(verbose){
        cat("\nItem",j,"Att",k,"\n")
        info <- cbind(add.new,remove.new)
        colnames(info) <- c("att","p(att)","GDI-largest att","GDI-largest vec",paste0("p-A",currentset))
        print(info)
      }


      #*********************************************************If additional elements should be added
      if (any(add.new[,2]<alpha.level)){ # Yes - sig some attributes should be required
        add.new <- add.new[which(add.new[,2]<alpha.level),,drop=FALSE] # which should be added
        remove.new <- remove.new[which(add.new[,2]<alpha.level),,drop=FALSE]
        loc.a <- which.max(add.new[,2+GDI])
        add.a <- add.new[loc.a,1]
        remove.a <- remove.new[loc.a,] # pvalues for all current attributes in the q-vector
        if(any(remove.a>alpha.level)) {
          remove.a <- currentset[which(remove.a>alpha.level)] # the att that needs to be removed
        }else{
          remove.a <- NULL
        }
        if(!stepwise) remove.a <- NULL
        currentset <- sort(setdiff(c(currentset,add.a),remove.a)) # all att. chosen
        att.monitor[[j]] <- c(att.monitor[[j]],add.a)

        if(is.null(remove.a)){ #if no attributes will be removed
          current_vec_PVAF <- add.new[which(add.new[,1,drop=FALSE]==add.a),4]
        }else{ # if some need to be removed - vector PVAF need to be recalculated
          newq <- vector("numeric",K)
          newq[currentset] <- 1
          current_vec_PVAF <- vsg0[j,rowMatch(patt,newq)$row.no]
        }

        if(current_vec_PVAF>=PVAF){
          loop <- FALSE
        }else{
          loop <- TRUE
          it <- it + 1
        }

      }else{   # No
        loop <- FALSE
      }
    }

    Q[j,] <- 0
    Q[j,currentset] <- 1

j <- j+1
  }

  if(seqent){
    ret.Q <- Qc
    ret.sugQ <- cbind(Qc[,1:2],Q)
  }else{
    ret.Q <- extract(GDINA.obj,"Q")
    ret.sugQ <- Q
  }
  out.vsg <- round(t(vsg), digits)
  out.PVAF <- round(t(vsg0), digits)
  rownames(out.vsg) <-
    rownames(out.PVAF) <-
    apply(patt, 1, paste, collapse = "")

  qvalid <-
    list(
      sug.Q = ret.sugQ,
      Q = ret.Q,
      varsigma = out.vsg,
      PVAF = out.PVAF,
      eps = PVAF,
      initialAtt=inichoose,
      est.p = est.p
    )

  return(qvalid)

}



Qval_PVAF <- function(GDINA.obj,
                 eps = 0.95,
                 digits = 4) {


  seqent <- extract(GDINA.obj,"sequential")

  if(seqent){
    dat <- extract(GDINA.obj,"seq.dat")
  }else{
    dat <- extract(GDINA.obj,"dat")
  }
  Y <- dat

  Qr <- Qc <- extract(GDINA.obj,"Qc")

  Q <- extract(GDINA.obj,"Q")

  N <- extract(GDINA.obj, "nobs")

  J <- extract(GDINA.obj, "nitem")

  K <- extract(GDINA.obj, "natt")

  L <- no_LC(Q)

  Kj <- rowSums(attributepattern(K)[-1, ])
  w <- extract(GDINA.obj, "posterior.prob") #1 x L
  YY <- Y
  YY[is.na(YY)] <- 0
  rc <- apply(YY, 2, function(x) {
    colSums(x * exp(extract(GDINA.obj, "logposterior.i")))
  })
  rn <- apply(1 * (!is.na(Y)), 2, function(x) {
    colSums(x * exp(extract(GDINA.obj, "logposterior.i")))
  })
  # est.p <- rc/c(w*N)
  est.p <- (rc + 1e-10) / (rn + 2* 1e-10)
  patt <- attributepattern(K)[-1, ]
  loc <- eta(patt) #2^K-1 x 2^K
  vsg <- varsigma(as.matrix(t(loc)), as.matrix(est.p), c(w))
  PVAF <- vsg / vsg[, L - 1]

      val_q <- NULL
    for (k in sort(unique(Kj))) {
      tmp <- PVAF[, which(Kj == k)]
      if (length(which(Kj == k)) == 1) {
        tmp[which(tmp > eps)] <- 1
        tmp[which(tmp <= eps)] <- 0
        val_q <- cbind(val_q, tmp)
      } else{
        val_q <- cbind(val_q, apply(tmp, 1, function(x) {
          ifelse (max(x) > eps, which.max(x), 0)
        }))
      }
    }
    if (ncol(val_q) > 1) {
      for (k in 2:ncol(val_q)) {
        val_q[which(val_q[, k] > 0), k] <-
          val_q[which(val_q[, k] > 0), k] + sum(Kj < k)
      }
    }
    #### modified Q-matrix and associated PVAF
    loc_q <- apply(val_q, 1, function(x) {
      x[which.max(x > 0)]
    })
    val_q <- attributepattern(K)[-1, ][loc_q, ]

  out.vsg <- round(t(vsg), digits)
  out.PVAF <- round(t(PVAF), digits)
  rownames(out.vsg) <-
    rownames(out.PVAF) <-
    apply(patt, 1, paste, collapse = "")
  # Q <- data.frame(Q, row.names = extract(GDINA.obj, "item.names"))
  if(seqent){
    ret.Q <- Qc
    ret.sugQ <- cbind(Qc[,1:2],val_q)
  }else{
    ret.Q <- Q
    ret.sugQ <- val_q
  }

  qvalid <-
    list(
      sug.Q = ret.sugQ,
      Q = ret.Q,
      varsigma = out.vsg,
      PVAF = out.PVAF,
      eps = eps,
      est.p = est.p
    )

  return(qvalid)
}
