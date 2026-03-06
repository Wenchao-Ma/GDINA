#' @include GDINA.R

SG.Est <- function(dat, Q, weight=NULL, model, sequential,att.dist, att.prior, saturated,
                att.str, mono.constraint, no.bugs, verbose,
                catprob.parm,loglinear,item.names,solnp_args,item.prior,
                linkfunc,higher.order, solver,auglag_args,nloptr_args,
                DesignMatrices,ConstrPairs,control){



  #################################
  #
  # data, Q-matrix, item.names
  #
  ##################################
  Qcm <- originalQ <- Q
  if(is.null(weight))
    weight <- rep(1, nrow(dat))
  del.ind <- NULL
  if (any(is.na(dat))) {
    # some missings individuals with one or fewer valid response are
    # removed
    del.ind <- which(rowSums(1L - is.na(dat)) < 2L, arr.ind = TRUE)
    if (length(del.ind) > 0L) {
      warning(length(del.ind)," individuals with one or fewer valid responses are removed.", call. = FALSE)
      dat <- dat[-del.ind,]
      weight <- weight[-del.ind]
    }

  }
  originalData <- dat

  dat <- unique(originalData)
  raw2unique <- match(apply(originalData,1,paste0,collapse=""),apply(dat,1,paste0,collapse=""))
  freq <- aggregate(weight,by = list(raw2unique),sum)$x



  N <- nrow(dat)
  nitems <- ncol(dat)
  ncat <- nrow(Q)


  model <- model2numeric(model, ncat)

  #     model.char model.num linkf.num linkf.char rule
  # 1    LOGGDINA        -3         3        log    0
  # 2  LOGITGDINA        -2         2      logit    0
  # 3         UDF        -1        -1        UDF   -1
  # 4       GDINA         0         1   identity    0
  # 5        DINA         1         1   identity    1
  # 6        DINO         2         1   identity    2
  # 7        ACDM         3         1   identity    3
  # 8         LLM         4         2      logit    3
  # 9        RRUM         5         3        log    3
  # 10     MSDINA         6         1   identity    4
  # 11    BUGDINO         7         1   identity    5
  # 12       SISM         8         1   identity    6

  if (sequential) {
    if(any(model >= 6))
      stop("The model speicified and sequential model cannot be estimated together.",call. = FALSE)
    dat <- seq_coding(dat, Q)
    Q <- Q[,-c(1, 2)]
    if (is.null(item.names))
      item.names <- paste("Item", originalQ[, 1], "Cat", originalQ[, 2])

  } else if (any(model == 6)) { #MS-DINA

    msQ <- unrestrQ(Q[which(model == 6), ])
    for (j in unique(msQ[, 1])) {
      Q[which(Q[, 1] == j & Q[, 2] == 1), ] <- msQ[which(msQ[, 1] == j & msQ[, 2] == 1), ]
      loc <- which(Q[, 1] == j & Q[, 2] != 1)
      if(length(loc)>0){
        Q <- Q[-loc, ]
        model <- model[-loc]
      }

    }
    Qcm <- Q
    Q <- Q[,-c(1, 2)]

    ncat <- nitems
  }else if(any(model == 7)){ # BUG-DINO
    no.bugs <- ncol(Q)
  }


  ##############################################
  #
  # model, condensation rule and link functions
  #
  ##############################################
  # rule: 0 -> saturated model; 1 ->DINA; 2 ->DINO; 3 ->additive model; 4 ->MS-DINA; -1 -> UDF
  rule <- model2rule(model)
  # identitiy link -> 1
  # logit link -> 2
  # log link -> 3
  LF.numeric <- linkf.numeric(linkfunc, model)


  if(length(mono.constraint)==1)
    mono.constraint <- rep(mono.constraint, ncat)


  if (is.null(item.names))
    item.names <- paste("Item", seq_len(nitems))

  if(!is.matrix(Q))
    Q <- as.matrix(Q)

  if(!is.matrix(dat))
    dat <- as.matrix(dat)



  #################################
  #
  # control arguments
  #
  ##################################

  control <- init_control(control, ncat, model, is_mg = FALSE)


  # input check
  inputcheck(dat = dat, Q = Q, model = model, sequential = sequential, att.dist = att.dist,loglinear=loglinear,
             no.bugs = no.bugs, verbose = verbose, catprob.parm = catprob.parm, mono.constraint = mono.constraint,
             att.prior = att.prior, lower.p = control$lower.p,upper.p = control$upper.p, att.str = att.str,
             nstarts = control$nstarts, conv.crit = control$conv.crit, maxitr = control$maxitr)



  #################################
  #
  # solver arguments
  #
  ##################################

  solver_init <- init_solver_args(auglag_args, solnp_args, nloptr_args)
  auglag_args <- solver_init$auglag_args
  solnp_args <- solver_init$solnp_args
  nloptr_args <- solver_init$nloptr_args
  
  solver <- init_solver(solver, ncat)


  #########################################
  #
  # att. space and reduced att. space
  #
  #########################################

  K <- ncol(Q)
  Kj <- rowSums(Q > 0)

  if(is.null(att.str)){ # no structure
    AlphaPattern <- as.matrix(att.structure(hierarchy.list = att.str,K = K,Q = Q,att.prob="uniform")$`att.str`)
    parloc <- eta(Q)  #J x L
    reduced.LG <- item_latent_group(Q)
  }else if(is.matrix(att.str)){
    AlphaPattern <- att.str
    parloc <- eta(Q, AlphaPattern)  #J x L
    reduced.LG <- item_latent_group(Q, AlphaPattern)
  }else{
    AlphaPattern <- as.matrix(att.structure(hierarchy.list = att.str,K = K,Q = Q,att.prob="uniform")$`att.str`)
    parloc <- eta(Q, AlphaPattern)  #J x L
    reduced.LG <- item_latent_group(Q, AlphaPattern)
  }
  L <- nrow(AlphaPattern)  # The number of latent classes


  Lj <- sapply(reduced.LG,nrow)


  #################################
  #
  # joint att. distribution
  #
  ##################################
  no.mg <- 1L
  gr <- rep(1L, N)
  gr.label <- "all data"

  att.dist <- tolower(att.dist)
  lambda <- NULL
  if(att.dist == "higher.order"){

    myHO <- list(model = "Rasch",nquad = ifelse(is.null(higher.order) | is.null(higher.order$nquad),25L,higher.order$nquad),
                 SlopeRange = c(0.1,5), InterceptRange = c(-4,4), Prior = FALSE,
                 SlopePrior = c(0,0.25), InterceptPrior = c(0L,1L), anchor = "all")
    myHO$QuadNodes = matrix(seq(-6,6,length.out = myHO$nquad),nrow = myHO$nquad,ncol = no.mg)
    myHO$QuadWghts = ColNormalize(dnorm(myHO$QuadNodes))

    higher.order <- modifyList(myHO,higher.order)

    lambda <- matrix(c(rep(1,K),rnorm(K,0,0.5)),ncol = 2)

  }else if(att.dist == "independent"){
    lambda <- matrix(c(rep(0,K),rnorm(K,0,0.5)),ncol = 2)
  }else if(att.dist == "saturated"){
    if(!is.null(att.prior))
      lambda <- att.prior
  }

  # Generate the log prior - a L x no.mg matrix
  if (is.null(att.prior)) {
    att.prior <- matrix(1/L, nrow = L, ncol = 1)
    logprior <- matrix(log(att.prior), nrow = L, ncol = 1)
  } else if (is.vector(att.prior)) {
    att.prior <- matrix(att.prior, ncol = 1) # vector -> matrix
  }else if (is.matrix(att.prior)) {
    if (nrow(att.prior) != L || ncol(att.prior) != 1)
      stop(
        "Joint attribute distribution priors must be a matrix of dimension 2^K by number of groups if specified.",
        call. = FALSE
      )
  }else if(is.data.frame(att.prior)){
    if (nrow(att.prior) != L || ncol(att.prior) != 1){
      stop(
        "Joint attribute distribution priors must be a matrix of dimension 2^K by number of groups if specified.",
        call. = FALSE
      )
    }else{
      att.prior <- as.matrix(att.prior)
    }
  }

  logprior <- log(ColNormalize(att.prior))



  #########################################
  #
  # constr and design matrices
  #
  #########################################



  ConstrType <- rep(1, ncat)
  ConstrMatrix <- vector("list", ncat)
  if(is.null(ConstrPairs)){
    ConstrPairs <- vector("list", ncat)
    for(j in seq_len(ncat)){
      if(mono.constraint[[j]]){
        ConstrType[j] <- 3
        ConstrPairs[[j]] <- partial_order2(Kj[j],reduced.LG[[j]])
        nctj <- nrow(ConstrPairs[[j]])
        tmp <- matrix(0,nctj,Lj[j])
        tmp[matrix(c(seq_len(nctj),ConstrPairs[[j]][,1]),ncol = 2)] <- 1
        tmp[matrix(c(seq_len(nctj),ConstrPairs[[j]][,2]),ncol = 2)] <- -1
        ConstrMatrix[[j]] <- tmp
      }
    }
  }

  if(is.null(DesignMatrices)){
    if(any(model == -1))
      stop("design.matrix must be provided for user-defined models.",call. = FALSE)
    DesignMatrices <-  vector("list", ncat)
    for(j in seq_len(ncat)) {
      if(model[j] == 6){
        DesignMatrices[[j]] <- designmatrix(model = model[j],Qj = originalQ[which(originalQ[,1]==j),-c(1:2),drop=FALSE])
      }else if(model[j] %in% c(7,8)){
        DesignMatrices[[j]] <- designmatrix(model = model[j],Qj = Q[j,],no.bugs=no.bugs)

      }else if(rule[j] >= 0 & rule[j]<= 3){

        DesignMatrices[[j]] <- designM(Kj[j], rule[j], reduced.LG[[j]])
      }
    }
  }



  #########################################
  #
  # initial catprob.parm
  #
  #########################################

  #print(DesignMatrices)
  if (is.null(catprob.parm)) {
    item.parm <- initials(Q, control$nstarts, DesignMatrices = DesignMatrices,att.str = att.str)  #a list with nstarts matrices of size J x 2^Kjmax
    ###### Multiple starting values
    if (control$nstarts > 1L) {
      neg2LL <- vector("numeric",control$nstarts)
      for (i in seq_len(control$nstarts)) {
        neg2LL[i] <- -2*ObsLogLik(mpar = as.matrix(item.parm[[i]]),
                                  mX = dat,
                                  vlogPrior = as.matrix(logprior),
                                  vgroup = rep(1,N),
                                  mloc = as.matrix(parloc),
                                  weights = freq)
      }
      item.parm <- item.parm[[which.min.randomtie(neg2LL)]]
    }
  } else {
    item.parm <- l2m(catprob.parm)
    neg2LL <- NA
  }
  initial.parm <- item.parm



  ##############################
  #
  # variable declarations
  #
  #############################
  npar.items <- sapply(DesignMatrices,ncol)
  stru.npar <- 0
  itr <- 0L

  delta <- calc_delta(item.parm, DesignMatrices = DesignMatrices, linkfunc = LF.numeric)

  parm0 <- list(ip = c(item.parm), prior = c(exp(logprior)), neg2LL = 0, delt = unlist(delta))

  dif.parm <- list(ip = 0,
                   prior = 0,
                   neg2LL = 0,
                   delt = 0)
  ##############################
  #
  #           E-M
  #
  #############################

  initial.logprior <- logprior

  if(control$Cpp){
     stopifnot(all(model %in% c(0,1,2)),all(mono.constraint==FALSE),att.dist == "saturated")

    ret <- fast_GDINA_EM(parloc, item.parm, dat, logprior, model,control$vmaxitr,control$lower.p,
                  control$upper.p,control$smallNcorrection,item.prior$beta,item.prior$on,control$conv.crit)
    item.parm <- ret$ip
    logprior <- ret$logprior
    itr <- ret$itr
    delta <- calc_delta(item.parm, DesignMatrices, LF.numeric)
    higher.order <- NULL
    lambda <- exp(logprior)
    stru.npar <- length(lambda) - 1
  }else{
    while(itr < control$maxitr)
    {

      estep <- LikNR(mpar = as.matrix(item.parm),
                     mX = dat,
                     vlogPrior = as.matrix(logprior),
                     vgroup = rep(1,N),
                     mloc = as.matrix(parloc),
                     weights = freq,
                     simplify = TRUE)



      # length of J indicating whether a nonzero cat should  be est. (TRUE) or not (FALSE)
      est.bin <- control$vmaxitr>itr
      # print(estep$Rg)
      optims <- Mstep(J=ncat, Kj = Kj, Lj = Lj, Rg = estep$Rg, Ng = estep$Ng, model = model, item.parm = item.parm, delta = delta,
                      ConstrType = ConstrType, correction = control$smallNcorrection, lower.p = control$lower.p, itr = itr + 1,
                      upper.p = control$upper.p, ConstrMatrix=ConstrMatrix,linkfunc=LF.numeric,MstepMessage = control$MstepMessage,
                      solver = solver, DesignMatrices = DesignMatrices, est.bin = est.bin, item.prior = item.prior,
                      auglag_args = auglag_args,solnp_args = solnp_args,nloptr_args = nloptr_args)

      item.parm <- optims$item.parm
      delta <- optims$delta


      struc.parm <- structural.parm.sg(AlphaPattern = AlphaPattern, logprior=estep$logprior,
                                       att.dist=att.dist,att.str=att.str, saturated = saturated,initial.logprior = initial.logprior,
                                       lower.prior = control$lower.prior,loglinear=loglinear,
                                       K=K,N=sum(freq),higher.order=higher.order,lambda = lambda)

      higher.order <- struc.parm$higher.order
      lambda <- struc.parm$lambda
      logprior <- struc.parm$logprior
      stru.npar <- struc.parm$npar
      parm1 <- list(ip = c(item.parm),
                    prior = c(exp(estep$logprior)),
                    neg2LL = -2 * estep$LL,
                    delt = unlist(delta))

      dif.parm <- list(ip = max(abs(parm1$ip-parm0$ip),na.rm = TRUE),
                       prior = max(abs(parm1$prior-parm0$prior),na.rm = TRUE),
                       neg2LL = parm0$neg2LL-parm1$neg2LL,
                       delt = max(abs(parm1$delt-parm0$delt),na.rm = TRUE))

      parm0 <- parm1
      itr <- itr + 1
      conv_check <- check_em_convergence(dif_parm = dif.parm, conv_type = control$conv.type,
                                         conv_crit = control$conv.crit, neg2LL_current = parm0$neg2LL)
      maxchg <- conv_check$maxchg

      print_em_progress(itr = itr, maxchg = maxchg, neg2LL = -2 * estep$LL, verbose = verbose)

      if(conv_check$converged) break
    }
  }



  logprior0 <- logprior # old log priors

  estep <- LikNR(mpar = as.matrix(item.parm),
                 mX = dat,
                 vlogPrior = as.matrix(logprior),
                 vgroup = rep(1,N),
                 mloc = as.matrix(parloc),
                 weights = freq,
                 simplify = FALSE)


  total.item.npar <- sum(sapply(DesignMatrices,ncol))

  free.item.npar <- 0

  if(any(control$vmaxitr>0)){
    est.item <- which(control$vmaxitr>0)
    free.item.npar <- sum(sapply(DesignMatrices[est.item],ncol))
  }




  npar <- free.item.npar + stru.npar

  neg2LL <- -2 * estep$LL

  item.prob <- vector("list",ncat)
  initial.parm <- m2l(initial.parm)
  IRF.labels <- lapply(reduced.LG,function(x)apply(x,1,paste0,collapse=""))
  for (j in seq_len(ncat)){
    item.prob[[j]] <- item.parm[j,1:Lj[j]]
    names(initial.parm[[j]]) <- names(item.prob[[j]]) <- paste0("P(",IRF.labels[[j]],")")
  }
  postP <- exp(t(estep$logprior))
  pf <- LC.Prob <- uP(parloc,item.parm)
  if(sequential){ # transform LC.prob
    p <- LC.Prob
    for (j in seq_len(nitems)){
      locj <- which(originalQ[,1]==j)
      Qj <- originalQ[locj,-c(1:2),drop=FALSE]
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

  LC.labels <- apply(AlphaPattern,1,paste0,collapse = "")
  names(item.prob) <- names(initial.parm) <- rownames(LC.Prob) <- names(delta) <- rownames(item.parm) <- item.names
  colnames(LC.Prob) <- colnames(pf) <- colnames(postP) <- LC.labels

  att.prior = c(exp(logprior0))

  list(catprob.parm = item.prob, delta.parm = delta, catprob.matrix = item.parm,
       struc.parm = lambda, model = model2character(model), LC.prob = LC.Prob,
       posterior.prob = postP, pf = pf, attributepattern = AlphaPattern,
       testfit = list(Deviance=neg2LL,npar = npar,item.npar = free.item.npar,
                      AIC=2 * npar + neg2LL,
                      BIC=neg2LL + npar * log(length(raw2unique)),
                      CAIC=neg2LL + npar * (log(length(raw2unique))+1),
                      AICc=neg2LL + 2*npar*(npar+1) / (length(raw2unique)-npar-1),
                      SABIC=neg2LL + npar * log((length(raw2unique)+2)/24)),
       technicals = list(logposterior.i = estep$logpost[raw2unique, ], loglikelihood.i = estep$loglik[raw2unique, ],
                         free.item.npar = free.item.npar,
                         total.item.npar = total.item.npar, stru.npar = stru.npar, total.npar = npar,
                         expectedCorrect = estep$Rg, expectedTotal = estep$Ng,initial.parm = initial.parm,
                         LC.labels = LC.labels,reduced.LG=reduced.LG,eta = parloc,del.ind=del.ind),
       options = list(dat = originalData, Q = originalQ, Qm = Q, Qcm = Qcm, model = model,
                      itr = itr, dif.LL = dif.parm$neg2LL,dif.p=dif.parm$ip,dif.prior=dif.parm$prior,
                      att.dist=att.dist, higher.order=higher.order,att.prior = att.prior, no.bugs = no.bugs,
                      mono.constraint = mono.constraint, item.names = item.names, group = rep(1,N), gr = gr,
                      att.str= att.str,  seq.dat = dat[raw2unique, ], no.group = 1, group.label = "all",
                      verbose = verbose, catprob.parm = catprob.parm,sequential = sequential,
                      nloptr_args = nloptr_args,auglag_args=auglag_args,solnp_args = solnp_args,
                      linkfunc = LF.numeric,higher.order = higher.order, loglinear = loglinear, solver = solver,
                      DesignMatrices = DesignMatrices,ConstrPairs = ConstrPairs),
       control = control)

}
