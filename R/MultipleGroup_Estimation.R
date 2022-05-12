#' @include GDINA.R
MG.Est <- function(dat, Q, model, sequential,att.dist, att.prior,saturated,
                att.str, mono.constraint, group, no.bugs, verbose,
                catprob.parm,loglinear,item.names,solnp_args,item.prior,
                linkfunc,higher.order, solver,auglag_args,nloptr_args,
                DesignMatrices,ConstrPairs,control){

  #######################
  #
  # Groups
  #
  #######################

  if (length(group)==1){
    if(is.positiveInteger(group) & group<ncol(dat))
      group <- dat[,group]
    dat <- dat[,-group]
  }else{
    if (nrow(dat) != length(group))
      stop("Group indicator variable is not correctly specified.", call. = FALSE)
  }

  #original group indicator
  ori.group <- group

  if(is.factor(group)){
    ori.gr.label <- levels(group)
    group <- as.vector.factor(group) # factor -> vector
  }else if(is.vector(group)){
    ori.gr.label <- unique(group)
  }

  group <- vector("numeric",length = length(ori.group))
  no.mg <- length(ori.gr.label) # the number of groups
  gr.label <- seq_len(no.mg)
  # group -> numeric vector
  if(!is.numeric(ori.group)){
    for(g in seq_len(no.mg)){
      group[ori.group==ori.gr.label[g]] <- g
    }
  }else{
    group <- ori.group
  }

  #################################
  #
  # data, Q-matrix, item.names
  #
  ##################################
  Qcm <- originalQ <- Q


  del.ind <- NULL
  if (any(is.na(dat))) {
    # some missings individuals with one or fewer valid response are
    # removed
    del.ind <- which(rowSums(1L - is.na(dat)) < 2L, arr.ind = TRUE)
    if (length(del.ind) > 0L) {
      warning(
        length(del.ind),
        " individuals with one or fewer valid responses are removed.",
        call. = FALSE
      )
      dat <- dat[-del.ind,]
      group <- group[-del.ind]
      ori.group <- ori.group[-del.ind]
    }
  }

  originalData <- dat
  N <- nrow(dat)
  nitems <- ncol(dat)
  ncat <- nrow(Q)
  Ng <- table(group)

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
    if(any(model == 6))
      stop("MSDINA model and sequential model cannot be estimated together.",call. = FALSE)
    dat <- seq_coding(dat, Q)
    Q <- Q[,-c(1, 2)]
    if (is.null(item.names))
      item.names <- paste("Item", originalQ[, 1], "Cat", originalQ[, 2])

  } else if (any(model == 6)) {

    msQ <- unrestrQ(Q[which(model == 6), ])
    for (j in unique(msQ[, 1])) {
      Q[which(Q[, 1] == j &
                Q[, 2] == 1), ] <- msQ[which(msQ[, 1] == j & msQ[, 2] == 1), ]
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


  if(!is.matrix(Q))
    Q <- as.matrix(Q)

  if(!is.matrix(dat))
    dat <- as.matrix(dat)


  #################################
  #
  # control arguments
  #
  ##################################


  myControl <- list(
    maxitr = 2000,
    conv.crit = 1e-4,
    conv.type = c("ip","mp"),
    nstarts = 3L,
    lower.p = 1e-4,
    upper.p = 1 - 1e-4,
    lower.prior = .Machine$double.eps,
    randomseed = 123456,
    smallNcorrection = c(.0005, .001),
    MstepMessage = FALSE,
    Cpp = FALSE,
    countitemparm = 0 # if an item parameter is fixed, it will not count as a parameter
  )


  control <- utils::modifyList(myControl,control)


  if (length(control$lower.p) == 1) {
    control$lower.p <- rep(control$lower.p, ncat)
  } else {
    if (length(control$lower.p) != ncat)
      stop("lower.p must have length of 1 or number of nonzero categories", call. = FALSE)
  }
  if (length(control$upper.p) == 1) {
    control$upper.p <- rep(control$upper.p, ncat)
  } else {
    if (length(control$upper.p) != ncat)
      stop("upper.p must have length of 1 or number of nonzero categories", call. = FALSE)
  }

  if(length(control$maxitr)==1L) {
    control$vmaxitr <- rep(control$maxitr,ncat)
  }else if(length(control$maxitr)!=ncat){
    warning("Length of maxitr must be equal to 1 or the number of nonzero categories.",call. = FALSE)
  }else{
    control$vmaxitr <- control$maxitr
    control$maxitr <- max(control$maxitr)
  }

  if(control$Cpp)
    warning("No fast implementation is available for multiple-group analysis.",call. = FALSE)

  control$nstarts <- ifelse(any(model>6),1L,3L) #if model is SISM or BUGDINO, only one set of starting values is used

  set.seed(control$randomseed)


  # input check
  inputcheck(dat = dat, Q = Q, model = model, sequential = sequential, att.dist = att.dist,
             no.bugs = no.bugs, verbose = verbose, catprob.parm = catprob.parm, mono.constraint = mono.constraint,
             att.prior = att.prior, lower.p = control$lower.p,upper.p = control$upper.p, att.str = att.str,
             nstarts = control$nstarts, conv.crit = control$conv.crit, maxitr = control$maxitr)



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



  #################################
  #
  # solver arguments
  #
  ##################################


  myAuglag_args <-
    list(control.outer = list(
      trace = FALSE,
      method = "nlminb",
      kkt2.check = FALSE,
      eps = 1e-6
    ))
  auglag_args <- modifyList(myAuglag_args, auglag_args)
  mySolnp_args <- list(trace = 0)
  solnp_args <- modifyList(mySolnp_args, solnp_args)
  Mynloptr_args <- list(xtol_rel = 1e-4)
  nloptr_args <- modifyList(Mynloptr_args, nloptr_args)

  if(is.null(solver)){
    solver <- rep("auto", ncat)
  }else if (length(solver) == 1) {
    solver <- rep(solver, ncat)
  } else {
    if (length(solver) != ncat)
      stop("solver must have length of 1 or number of items.", call. = FALSE)
  }


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


  att.dist <- tolower(att.dist)

  if (length(att.dist) == 1)
    att.dist <- rep(att.dist, no.mg)
  if (length(control$lower.prior) == 1)
    control$lower.prior <- rep(control$lower.prior, no.mg)
  if (length(loglinear) == 1)
    loglinear <- rep(loglinear, no.mg)
  if (any(att.dist == "higher.order")){
    if(any(att.dist != "higher.order"))
      stop("Higher-order model must be used for all groups.",call. = FALSE)
    att.dist <- rep("higher.order", no.mg)
  }

  lambda <- vector("list",no.mg)

  if(any(att.dist=="higher.order")){

    myHO <- list(model = "Rasch",nquad = 25L, SlopeRange = c(0.1,5), InterceptRange = c(-4,4), Prior = FALSE,
                 SlopePrior = c(0,0.25), InterceptPrior = c(0L,1L), anchor = "all")
    # GH <- gaussHermiteData(myHO$nquad)
    # myHO$QuadNodes = matrix(GH$x*sqrt(2),nrow = myHO$nquad,ncol = no.mg)
    # myHO$QuadWghts = matrix(GH$w/sqrt(pi),nrow = myHO$nquad,ncol = no.mg)
    myHO$QuadNodes = matrix(seq(-4,4,length.out = myHO$nquad),nrow = myHO$nquad,ncol = no.mg)
    myHO$QuadWghts = ColNormalize(dnorm(myHO$QuadNodes))

    higher.order <- modifyList(myHO,higher.order)

    for(g in seq_len(no.mg))
      lambda[[g]] <- matrix(c(rep(1,K),rnorm(K,0,0.5)),ncol = 2)

  }else{
    higher.order <- NULL
    if(any(att.dist=="independent")){
      for(g in seq_len(no.mg)){
        if(att.dist[g]=="independent")
          lambda[[g]] <- matrix(c(rep(0,K),rnorm(K,0,0.5)),ncol = 2)
      }
    }

  }

  # Generate the log prior - a L x no.mg matrix
  if (is.null(att.prior)) {
    att.prior <- matrix(1/L, nrow = L, ncol = no.mg)
    logprior <- matrix(log(att.prior), nrow = L, ncol = no.mg)
  } else if (is.vector(att.prior)) {
    att.prior <- matrix(att.prior, ncol = no.mg) # vector -> matrix
  }else if (is.matrix(att.prior)) {
    if (nrow(att.prior) != L || ncol(att.prior) != no.mg)
      stop(
        "Joint attribute distribution priors must be a matrix of dimension 2^K by number of groups if specified.",
        call. = FALSE
      )
  }else if(is.data.frame(att.prior)){
    if (nrow(att.prior) != L || ncol(att.prior) != no.mg){
      stop(
        "Joint attribute distribution priors must be a matrix of dimension 2^K by number of groups if specified.",
        call. = FALSE
      )
    }else{
      att.prior <- as.matrix(att.prior)
    }

  }
  if(!is.matrix(saturated$prior)|!is.data.frame(saturated$prior))
    saturated$prior <- matrix(saturated$prior,ncol = no.mg)


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

  if (is.null(catprob.parm)) {
    item.parm <- initials(Q, control$nstarts, DesignMatrices = DesignMatrices)  #a list with nstarts matrices of size J x 2^Kjmax
    ###### Multiple starting values
    if (control$nstarts > 1L) {
      neg2LL <- vector("numeric",control$nstarts)
      for (i in seq_len(control$nstarts)) {
        neg2LL[i] <- -2*ObsLogLik(mpar = as.matrix(item.parm[[i]]),
                                  mX = as.matrix(dat),
                                  vlogPrior = as.matrix(logprior),
                                  vgroup = as.vector(group),
                                  mloc = as.matrix(parloc),
                                  weights = rep(1,N))
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
  stru.npar <- 0
  while(itr < control$maxitr)
  {
    estep <- LikNR(mpar = as.matrix(item.parm),
                   mX = as.matrix(dat),
                   vlogPrior = as.matrix(logprior),
                   vgroup = group,
                   mloc = as.matrix(parloc),
                   weights = rep(1,N),
                   simplify = TRUE)



    # length of J indicating whether a nonzero cat should  be est. (TRUE) or not (FALSE)
    est.bin <- control$vmaxitr>itr
    # print(estep$Rg)
    optims <- Mstep(J=ncat, Kj = Kj, Lj = Lj,Rg = estep$Rg, Ng = estep$Ng, model = model, item.parm = item.parm, delta = delta,
                    ConstrType = ConstrType, correction = control$smallNcorrection, lower.p = control$lower.p, itr = itr + 1,
                    upper.p = control$upper.p, ConstrMatrix=ConstrMatrix,linkfunc=LF.numeric,MstepMessage = control$MstepMessage,
                    solver = solver, DesignMatrices = DesignMatrices, est.bin = est.bin, item.prior = item.prior,
                    auglag_args = auglag_args,solnp_args = solnp_args,nloptr_args = nloptr_args)

    item.parm <- optims$item.parm
    delta <- optims$delta


    struc.parm <- structural.parm.mg(AlphaPattern = AlphaPattern, no.mg = no.mg, logprior=estep$logprior,
                                  att.dist=att.dist,att.str=att.str, saturated = saturated,initial.logprior = initial.logprior,
                                  lower.prior = control$lower.prior,loglinear=loglinear,
                                  Ng=Ng,K=K,N=N,higher.order=higher.order,lambda = lambda)

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
    maxchg <- 0
    if(any(tolower(control$conv.type)=="ip"))  maxchg <- max(maxchg,dif.parm$ip)
    if(any(tolower(control$conv.type)=="delta"))  maxchg <- max(maxchg,dif.parm$delt)
    if(any(tolower(control$conv.type)=="mp"))  maxchg <- max(maxchg,dif.parm$prior)
    if(any(tolower(control$conv.type)=="neg2ll"))  maxchg <- max(maxchg,abs(dif.parm$neg2LL))
    if(any(tolower(control$conv.type)=="relneg2ll"))  maxchg <- max(maxchg,abs(dif.parm$neg2LL)/parm0$neg2LL)

    if(verbose==1L) {
      cat('\rIter =',itr,' Max. abs. change =',formatC(maxchg,digits = 5, format = "f"),
          ' Deviance  =',formatC(-2 * estep$LL,digits = 2, format = "f"),'                                                                                 ')
    }else if (verbose==2L) {
      cat('Iter =',itr,' Max. abs. change =',formatC(maxchg,digits = 5, format = "f"),
          ' Deviance  =',formatC(-2 * estep$LL,digits = 2, format = "f"),'                                                                                \n')
    }

    if(maxchg < control$conv.crit) break
  }

  logprior0 <- logprior # old log priors

  estep <- LikNR(mpar = as.matrix(item.parm),
                 mX = as.matrix(dat),
                 vlogPrior = as.matrix(logprior),
                 vgroup = group,
                 mloc = as.matrix(parloc),
                 weights = rep(1,N),
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
    for (j in seq_len(ncat)){
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

  rownames(postP) <- paste("Group",ori.gr.label)

  att.prior = exp(logprior0)


  list(catprob.parm = item.prob, delta.parm = delta, catprob.matrix = item.parm,
       struc.parm = lambda, model = model2character(model), LC.prob = LC.Prob,
       posterior.prob = postP, pf = pf,attributepattern = AlphaPattern,
       testfit = list(Deviance=neg2LL,npar = npar,item.npar = free.item.npar,
                      AIC=2 * npar + neg2LL, BIC=neg2LL + npar * log(N),
                      CAIC=neg2LL + npar * (log(N)+1),
                      SABIC=neg2LL + npar * log((N+2)/24)),
       technicals = list(logposterior.i = estep$logpost, loglikelihood.i = estep$loglik, free.item.npar = free.item.npar,
                         total.item.npar = total.item.npar, stru.npar = stru.npar, total.npar = npar,
                         expectedCorrect = estep$Rg, expectedTotal = estep$Ng,initial.parm = initial.parm,
                         LC.labels = LC.labels,reduced.LG=reduced.LG,eta = parloc,del.ind=del.ind),
       options = list(dat = originalData, Q = originalQ, Qm = Q, Qcm = Qcm, model = model,
                      itr = itr, dif.LL = dif.parm$neg2LL,dif.p=dif.parm$ip,dif.prior=dif.parm$prior,
                      att.dist=att.dist, higher.order=higher.order,att.prior = att.prior, no.bugs = no.bugs,
                      mono.constraint = mono.constraint, item.names = item.names,group = ori.group, gr = group,
                      att.str= att.str,  seq.dat = dat, no.group = no.mg, group.label = ori.gr.label,
                      verbose = verbose, catprob.parm = catprob.parm,sequential = sequential,
                      nloptr_args = nloptr_args,auglag_args=auglag_args,solnp_args = solnp_args,
                      linkfunc = LF.numeric,higher.order = higher.order, loglinear = loglinear, solver = solver,
                      DesignMatrices = DesignMatrices,ConstrPairs = ConstrPairs),
       control = control)

}
