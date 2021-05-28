#' Q-matrix validation
#'
#' Q-matrix validation for the (sequential) G-DINA model based on PVAF (de la Torre & Chiu, 2016; Najera, Sorrel, & Abad, 2019; Najera et al., 2020), stepwise Wald test (Ma & de la Torre, 2020) or mesa plot (de la Torre & Ma, 2016).
#' All these methods are suitable for dichotomous and ordinal response data. If too many modifications are suggested based on the default PVAF method, you are suggested to try the stepwise Wald test method, iterative procedures or predicted cutoffs.
#' You should always check the mesa plots for further examination.
#'
#' @param GDINA.obj an estimated model object of class \code{GDINA}
#' @param method which Q-matrix validation method is used? Can be either \code{"PVAF"} or \code{"wald"}.
#' @param iter implement the method iteratively? Can be \code{"none"} for non-iterative validation (by default), \code{"test"}, \code{"test.att"}, or \code{"item"} (Najera et al., 2020).
#' @param eps cutoff value for PVAF from 0 to 1. Default = 0.95. Note that it can also be -1, indicating the predicted cutoff based on Najera, Sorrel, and Abad (2019).
#' @param digits how many decimal places in each number? The default is 4.
#' @param wald.args a list of arguments for the stepwise Wald test method.
#' \describe{
#' \item{SE.type}{type of covariance matrix for the Wald test}
#' \item{alpha.level}{alpha level for the wald test}
#' \item{GDI}{it can be 0, 1 or 2; 0 means GDI is not used to choose the attribute -
#' when more than one attributes are significant, the one with the largest p-value will be selected;
#' GDI=1 means the attribute with the largest GDI will be selected; GDI=2 means the q-vector with
#' the largest GDI will be selected.}
#' \item{verbose}{print detailed information or not?}
#' \item{stepwise}{\code{TRUE} for stepwise approach and \code{FALSE} for forward approach}
#' }
#' @param iter.args a list of arguments for the iterative implementation.
#' \describe{
#' \item{empty.att}{can a Q-matrix with an empty attribute (i.e., measured by no items) be provided? Default is FALSE}
#' \item{max.iter}{maximum number of iterations. Default is 150}
#' \item{verbose}{print information after each iteration? Default is FALSE}
#' }
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
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu}, \cr Miguel A. Sorrel, Universidad Autónoma de Madrid, \cr Jimmy de la Torre, The University of Hong Kong}
#' @references
#' de la Torre, J. & Chiu, C-Y. (2016). A General Method of Empirical Q-matrix Validation. \emph{Psychometrika, 81}, 253-273.
#'
#' de la Torre, J., & Ma, W. (2016, August). Cognitive diagnosis modeling: A general framework approach and its implementation in R. A Short Course at the Fourth Conference on Statistical Methods in Psychometrics, Columbia University, New York.
#'
#' Ma, W., & de la Torre, J. (2020). An empirical Q-matrix validation method for the sequential G-DINA model. \emph{British Journal of Mathematical and Statistical Psychology, 73}, 142-163.
#'
#' Najera, P., Sorrel, M. A., & Abad, F.J. (2019). Reconsidering cutoff points in the general method of empirical Q-matrix validation. \emph{Educational and Psychological Measurement, 79}, 727-753.
#'
#' Najera, P., Sorrel, M. A., de la Torre, J., & Abad, F. J. (2020). Improving robustness in Q-matrix validation using an iterative and dynamic procedure. \emph{Applied Psychological Measurement}.
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
#'
#' #Iterative procedure (test-attribute level)
#' pvaf3 <- Qval(mod1, method = "PVAF", eps = -1,
#'               iter = "test.att", iter.args = list(verbose = 1))
#' pvaf3
#'
#' ################################
#' #
#' # Ordinal response
#' #
#' ################################
#' seq.est <- GDINA(sim20seqGDINA$simdat,sim20seqGDINA$simQ, sequential = TRUE)
#' stepwise <- Qval(seq.est, method = "wald")
#'}

Qval <- function(GDINA.obj, method = "PVAF", iter = "none", eps = 0.95,
                    digits = 4, wald.args = list(),
                    iter.args = list(empty.att = FALSE, max.iter = 150, verbose = FALSE)){

  if (class(GDINA.obj) != "GDINA")
    stop("GDINA.obj must be a GDINA object from GDINA function.", call. = FALSE)

  if (!is.null(extract(GDINA.obj, "att.str")))
    stop("Q-matrix validation is not available if attributes are structured.",
         call. = FALSE)

  if(method == "wald" & iter != "none"){
    warning("Iterative implementation is not available for the Wald method.")
  }

  if(eps != -1){
    if (eps > 1 || eps < 0) stop("eps must be between 0 and 1, or equal to -1.", call. = FALSE)
  }

  if(!(iter %in% c("none", "test", "test.att", "item"))){stop("iter must be 'none', 'test', 'test.att', or 'item'.")}

  if(extract(GDINA.obj,"ngroup")>1) stop("Only available for single-group models.",call. = FALSE)

  # if(any(extract(GDINA.obj,"models_numeric")<0)||any(extract(GDINA.obj,"models_numeric")>5))
  #   stop("Models must be GDINA, DINA, DINO, ACDM, LLM or RRUM",call. = FALSE)

  if (max(extract(GDINA.obj, "Q")) > 1)
    stop("Q-validation can only be used for dichotomous attribute G-DINA model.",
         call. = FALSE)

  updated.wald.args <- NULL
  if(toupper(method)=="PVAF"){

    if(is.null(iter.args$empty.att)){iter.args$empty.att <- FALSE}
    if(is.null(iter.args$max.iter)){iter.args$max.iter <- 150}
    if(is.null(iter.args$verbose)){iter.args$verbose <- FALSE}
    iter.args <- iter.args[c("empty.att", "max.iter", "verbose")]
    # if(any(extract(GDINA.obj,"models")!="GDINA"))
    #   warning("Saturated G-DINA model may be used to calibrate all items for better performance.",call. = FALSE)
    ret <- Qval_PVAF(GDINA.obj,eps = eps, digits = digits,
                     iter = iter, iter.args = iter.args)

  }else if (toupper(method)=="WALD"){
    if(any(extract(GDINA.obj,"models")!="GDINA"))
      stop("Saturated G-DINA model needs to be fitted to all items.",call. = FALSE)
    if(eps == -1){
      gs <- coef(GDINA.obj, what = "gs") # item parameters
      eps <- plogis(-0.4045140782147791 +
                      4.8404850955032684E-4*extract(GDINA.obj,"nobs") +
                      2.8667570118638275*(1-sum(colMeans(gs))) +
                      -0.003315555999671906*extract(GDINA.obj,"nitem"))
    }
    args.default <- list(GDINA.obj = GDINA.obj, SE.type = 2,
                         alpha.level = 0.05, GDI = 2, PVAF = eps,
                         verbose = FALSE, stepwise = TRUE,digits=digits)
    updated.wald.args <- modifyList(args.default,wald.args)
    ret <- do.call(Qval_wald,updated.wald.args)
  }
  ret$method <- method
  ret$iter <- iter
  ret$iter.args <- iter.args
  ret$wald.args <- updated.wald.args
  ret$sequential <- extract(GDINA.obj,"sequential")
  class(ret) <- "Qval"
  return(ret)
}

Qval_wald <- function(GDINA.obj, SE.type = 2,
                      alpha.level = 0.05, GDI = 2, PVAF = 0.95,
                      verbose = FALSE, stepwise = TRUE, digits = 4,...){
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
    cat("First selected attribute based on GDI is",apply(vsgK, 1, which.max.randomtie))
  }
  #first attribute is the one with the largest GDI
  first.att <- apply(vsgK, 1, which.max.randomtie)
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
                      digits = 4,
                      iter = "none",
                      iter.args = list(empty.att = FALSE, max.iter = 150, verbose = FALSE)){

  seqent <- extract(GDINA.obj,"sequential")
  eps.base <- eps

  fit <- GDINA.obj

  Qiter <- list()
  if(seqent){
    Y <- extract(GDINA.obj,"seq.dat")
    dat <- extract(GDINA.obj,"dat")
    Qiter[["0"]] <- Qtest <- Qprov <- Q <- Qr <- Qc <- extract(GDINA.obj,"Qc")
  }else{
    Y <- dat <- extract(GDINA.obj,"dat")
    Qiter[["0"]] <- Qtest <- Qprov <- Q <- extract(GDINA.obj,"Q")
  }

  model <- GDINA.obj$options$model
  mono.constraint <- GDINA.obj$options$mono.constraint
  if(!is.null(GDINA.obj$options$higher.order)){
    ho.model <- GDINA.obj$options$higher.order$model
    ho.nquad <- GDINA.obj$options$higher.order$nquad
    ho.sRange <- GDINA.obj$options$higher.order$SlopeRange
    ho.iRange <- GDINA.obj$options$higher.order$InterceptRange
    ho.sPrior <- GDINA.obj$options$higher.order$SlopePrior
    ho.iPrior <- GDINA.obj$options$higher.order$InterceptPrior
    ho.prior <- GDINA.obj$options$higher.order$Prior
  }

  N <- extract(GDINA.obj, "nobs")
  J <- extract(GDINA.obj, "nitem")
  K <- extract(GDINA.obj, "natt")
  if(seqent){
    L <- no_LC(Q[,-c(1:2)])
    Js <- nrow(Q)
  } else {
    L <- no_LC(Q)
  }

  M <- attributepattern(K)[-1,]
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

  i <- 0

  if(iter == "none"){
    if(eps.base == -1){
      gs <- coef(fit, what = "gs") # item parameters
      eps <- plogis(-0.4045140782147791 +
                      4.8404850955032684E-4*extract(fit,"nobs") +
                      2.8667570118638275*(1-sum(colMeans(gs))) +
                      -0.003315555999671906*extract(fit,"nitem"))
    }
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

    i <- n.iter <- convergence <- 1

    qvalid <-
      list(
        sug.Q = ret.sugQ,
        Q = ret.Q,
        varsigma = out.vsg,
        PVAF = out.PVAF,
        eps = eps,
        est.p = est.p,
        n.iter = n.iter,
        convergence = convergence
      )

  } else {
    while(i < iter.args$max.iter){
      i <- i + 1

      if(eps.base == -1){
        gs <- coef(fit, what = "gs") # item parameters
        eps <- plogis(-0.4045140782147791 +
                        4.8404850955032684E-4*extract(fit,"nobs") +
                        2.8667570118638275*(1-sum(colMeans(gs))) +
                        -0.003315555999671906*extract(fit,"nitem"))
      }

      w <- extract(fit, "posterior.prob") #1 x L

      rc <- apply(YY, 2, function(x) {
        colSums(x * exp(extract(fit, "logposterior.i")))
      })
      rn <- apply(1 * (!is.na(Y)), 2, function(x) {
        colSums(x * exp(extract(fit, "logposterior.i")))
      })
      # est.p <- rc/c(w*N)
      est.p <- (rc + 1e-10) / (rn + 2* 1e-10)

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

      Qsug.obj <-
        list(
          sug.Q = ret.sugQ,
          Q = ret.Q,
          varsigma = out.vsg,
          PVAF = out.PVAF,
          eps = eps,
          est.p = est.p)

      Qsug <- Qsug.obj$sug.Q
      res.eps <- Qsug.obj$eps

      if(seqent){
        maxloc <- matrix(NA, nrow = K, ncol = Js, dimnames = list(paste0("K", 1:K), paste0("Js", 1:Js)))
      } else {
        maxloc <- matrix(NA, nrow = K, ncol = J, dimnames = list(paste0("K", 1:K), paste0("J", 1:J)))
      }
      for(kj in 1:K){maxloc[kj,] <- apply(Qsug.obj$PVAF[Kj == kj,, drop = FALSE], 2, which.max)}
      best.pos <- cumsum(c(0, table(Kj)[-K])) + maxloc
      if(seqent){
        best.PVAF <- sapply(1:Js, function(j) Qsug.obj$PVAF[best.pos[,j],j])
        rownames(best.PVAF) <- paste0("K", 1:K); colnames(best.PVAF) <- paste0("Js", 1:Js)
      } else {
        best.PVAF <- sapply(1:J, function(j) Qsug.obj$PVAF[best.pos[,j],j])
        rownames(best.PVAF) <- paste0("K", 1:K); colnames(best.PVAF) <- paste0("J", 1:J)
      }
      res.PVAF <- Qsug.obj$PVAF
      res.varsigma <- Qsug.obj$varsigma
      if(seqent){
        hit <- hit.cand <- which(sapply(1:Js, function(j) any(Qprov[j,] != Qsug[j,])))
      } else {
        hit <- hit.cand <- which(sapply(1:J, function(j) any(Qprov[j,] != Qsug[j,])))
      }

      if(length(hit.cand) == 0){
        i <- i - 1
        convergence <- 1 # YES convergence
        break
      } else {
        if(seqent){
          Kj.prov <- rowSums(Qprov[hit, -c(1:2), drop = F])
          Kj.sug <- rowSums(Qsug[hit, -c(1:2), drop = F])
        } else {
          Kj.prov <- rowSums(Qprov[hit,, drop = F])
          Kj.sug <- rowSums(Qsug[hit,, drop = F])
        }

        if(iter == "test.att"){
          att.change <- ifelse(Kj.prov < Kj.sug, 1, ifelse(Kj.prov == Kj.sug, 0, -1))
          new.att <- Kj.prov + att.change
          for(r in 1:length(hit)){
            att.r <- new.att[r]
            if(seqent){
              Qsug[hit[r], -c(1:2)] <- M[best.pos[att.r, hit[r]],]
            } else {
              Qsug[hit[r],] <- M[best.pos[att.r, hit[r]],]
            }
          }
        }

        if(iter == "item"){
          if(seqent){
            prov.q <- match(apply(Qprov[hit.cand, -c(1:2), drop = F], 1, paste, collapse = ""), apply(M, 1, paste, collapse = ""))
          } else {
            prov.q <- match(apply(Qprov[hit.cand,, drop = F], 1, paste, collapse = ""), apply(M, 1, paste, collapse = ""))
          }
          PVAF.diff <- c()
          for(j in 1:length(hit.cand)){
            PVAF.diff <- c(PVAF.diff, max(abs(best.PVAF[Kj.sug[j], hit.cand[j]] - best.PVAF[Kj.prov[j], hit.cand[j]]),
                                          abs(best.PVAF[Kj.sug[j], hit.cand[j]] - res.PVAF[prov.q[j], hit.cand[j]])))
          }
          hit <- hit.cand[which.max(PVAF.diff)]
        }
      }

      Qtest[hit,] <- Qsug[hit,]
      if(any(sapply(1:i, function(x) mean(Qtest == Qiter[[x]])) == 1)){convergence <- 2; break} # NO convergence (loop)
      if(any(colSums(Qtest) == 0)){convergence <- 3; break} # NO convergence (null attribute)
      Qprov <- Qiter[[as.character(i)]] <- Qtest
      if(iter.args$verbose){cat("Iteration =", sprintf("%03d", i), "| Item/s modified =", paste(hit, collapse = ", "), "\n")}
      if(is.null(GDINA.obj$options$higher.order)){
        fit <- GDINA(dat, Qprov, model, mono.constraint = mono.constraint, sequential = seqent, verbose = 0)
      } else {
        fit <- GDINA(dat, Qprov, model, mono.constraint = mono.constraint, att.dist = "higher.order", sequential = seqent, verbose = 0,
                     higher.order = list(model = ho.model, nquad = ho.nquad, SlopeRange = ho.sRange, InterceptRange = ho.iRange,
                                         SlopePrior = ho.sPrior, InterceptPrior = ho.iPrior, Prior = ho.prior))
      }
    }

    n.iter <- i
    if(n.iter == iter.args$max.iter){convergence <- 4} # NO convergence (max.iter reached)
    if(iter.args$empty.att & convergence == 3){
      Qsug <- Qtest
    } else {
      Qsug <- Qprov
    }

    qvalid <-
      list(
        # sug.Q = ret.sugQ,
        sug.Q = Qsug,
        Q = ret.Q,
        varsigma = out.vsg,
        PVAF = out.PVAF,
        eps = eps,
        est.p = est.p,
        n.iter = n.iter,
        convergence = convergence
      )
  }

  return(qvalid)
}
