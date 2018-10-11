Est <- function(dat, Q, model, sequential,att.dist, att.prior,saturated,
                att.str, mono.constraint, group, latent.var, verbose,
                catprob.parm,loglinear,item.names,solnp_args,item.prior,
                linkfunc,higher.order, solver,auglag_args,nloptr_args,
                DesignMatrices,ConstrPairs,control){

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
  if(att.str){
    if(!is.null(control$lower.prior)){
      if(control$lower.prior!=0){
        warning("lower.prior must be 0 when att.str = TRUE.",call. = FALSE)
        lower.prior <- 0
      }else{
        lower.prior <- 0
      }
    }else{
      lower.prior <- 0
    }

    if(mono.constraint){
      warning("Monotonic constraint cannot be imposed when att.str = TRUE.",call. = FALSE)
      mono.constraint <- FALSE
    }
  }else{
      lower.prior <- .Machine$double.eps
      }
  myControl <- list(
    maxitr = 2000,
    conv.crit = 1e-4,
    conv.type = c("ip","mp"),
    nstarts = 3L,
    lower.p = 1e-4,
    upper.p = 1 - 1e-4,
    lower.prior = lower.prior,
    randomseed = 123456,
    smallNcorrection = c(.0005, .001),
    MstepMessage = FALSE,
    countitemparm = 0 # if an item parameter is fixed, it will not count as a parameter
  )

  control <- utils::modifyList(myControl,control)

  set.seed(control$randomseed)


  if (is.null(group)){
    no.mg <- 1L
    gr <- rep(1L,nrow(dat))
    gr.label <- "all data"

  }else{
    if (nrow(dat) != length(group) | min(group)!=1 | !all(is.positiveInteger(group))) stop("Group indicator variable is not correctly specified.", call. = FALSE)
      gr <- group # group indicator variable
    }

    gr.label <- unique(gr) # group labels
    no.mg <- length(gr.label) # the number of groups

    if (no.mg > 1) {
      if (length(att.dist) == 1)
        att.dist <- rep(att.dist, no.mg)
      # if (length(att.str) == 1)
      #   att.str <- rep(att.str, no.mg)
      if (length(control$lower.prior) == 1)
        control$lower.prior <- rep(control$lower.prior, no.mg)
      if (length(loglinear) == 1)
        loglinear <- rep(loglinear, no.mg)
      if (any(att.dist == "higher.order")){
        if(any(att.dist != "higher.order")) stop("Higher-order model must be used for all groups.",call. = FALSE)
        att.dist <- rep("higher.order", no.mg)
      }

    }

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
        gr <- gr[-del.ind]
      }
    }
    originalQ <- Q
    originalData <- dat
    if (sequential)
      dat <- seq_coding(dat, Q)

    Ng <- N <- nrow(dat)

    if (no.mg > 1)
      Ng <- table(gr)

    J <- ncol(dat)

    M <- c("logGDINA","logitGDINA","UDF", "GDINA", "DINA", "DINO", "ACDM", "LLM", "RRUM", "MSDINA")

    model <- model.transform(model, nrow(Q))

    if (any(model == 6)) {
      #MSDINA
      msQ <- unrestrQ(Q[which(model == 6), ])
      for (j in unique(msQ[, 1])) {
        Q[which(Q[, 1] == j &
                  Q[, 2] == 1), ] <- msQ[which(msQ[, 1] == j & msQ[, 2] == 1), ]
        loc <- which(Q[, 1] == j & Q[, 2] != 1)
        Q <- Q[-loc, ]
        model <- model[-loc]
      }
    }

    Qcm <- Q

  model.names <- M[model + 4]
  names(model.names) <- item.names

  inputcheck(
    dat = dat,
    Q = Q,
    model = model,
    sequential = sequential,
    att.dist = att.dist,
    latent.var = latent.var,
    verbose = verbose,
    catprob.parm = catprob.parm,
    mono.constraint = mono.constraint,
    att.prior = att.prior,
    lower.p = control$lower.p,
    upper.p = control$upper.p,
    att.str = att.str,
    nstarts = control$nstarts,
    conv.crit = control$conv.crit,
    maxitr = control$maxitr
  )

  if(att.str&(any(model<0)|any(model>2)))stop("Attribute structure can only be imposed for the DINA, DINO and G-DINA models.",call. = FALSE)

  if (sequential) {
    Q <- Q[,-c(1, 2)]
    if (is.null(item.names))
      item.names <- paste("Item", originalQ[, 1], "Cat", originalQ[, 2])

  } else {
    if(any(model==6)) Q <- Q[,-c(1, 2)]
    if (max(dat, na.rm = TRUE) > 1)
      stop("Maximum response is greater than 1 - set sequential = TRUE to fit a sequential model.", call. = FALSE)
    if (is.null(item.names)) {
      item.names <- paste("Item", 1:ncol(dat))
    }

  }

  K <- ncol(Q)

  Kj <- rowSums(Q > 0)  # vector with length of J

  Lj <- 2 ^ Kj

  L <- no_LC(Q)  # The number of latent groups

  parloc <- eta(as.matrix(Q))  #J x L

  AlphaPattern <- attributepattern(Q = Q)



  if(length(mono.constraint)==1) mono.constraint <- rep(mono.constraint,J)

  ConstrType <- rep(1,J)
  ConstrMatrix <- vector("list",J)
  if(is.null(ConstrPairs)) ConstrPairs <- vector("list",J)

  if(is.null(DesignMatrices)){
    if(any(model==-1)) stop("design.matrix must be provided for user-defined models.",call. = FALSE)
    DesignMatrices <-  vector("list",J)
  }
  for(j in seq_len(J)) {
    if(model[j]>=0&model[j]<=5){
      DesignMatrices[[j]] <- designmatrix(Kj[j],model[j])
    }else if(model[j]==6){
      DesignMatrices[[j]] <- designmatrix(model = model[j],Qj = originalQ[which(originalQ[,1]==j),-c(1:2),drop=FALSE])
    }else if(model[j]==-2||model[j]==-3){
      DesignMatrices[[j]] <- designM(Kj[j],0)
    }
    if(mono.constraint[[j]]){
      ConstrType[j] <- 3
      ConstrPairs[[j]] <- partial_order2(Kj[j])
      nctj <- nrow(ConstrPairs[[j]])
      tmp <- matrix(0,nctj,2^Kj[j])
      tmp[matrix(c(seq_len(nctj),ConstrPairs[[j]][,1]),ncol = 2)] <- 1
      tmp[matrix(c(seq_len(nctj),ConstrPairs[[j]][,2]),ncol = 2)] <- -1
      ConstrMatrix[[j]] <- tmp
    }
  }


  if(is.null(linkfunc)){
    linkfunc <- rep(1,J)
    for(j in seq_len(J)){
      if(model[j] == 4 || model[j] == -2) linkfunc[j] <- 2 else
        if(model[j] == 5 || model[j] == -3) linkfunc[j] <- 3
    }
  }else if (length(linkfunc) == 1){
    linkfunc <- which(tolower(linkfunc)==c("identity","logit","log"))
    linkfunc <- rep(linkfunc, J)
  } else {
    if (length(linkfunc) != J) stop("linkfunc must have length of 1 or J.", call. = FALSE)
    tmp <- linkfunc
    linkfunc <- rep(1,J)
    linkfunc[which(tolower(tmp)=="logit")] <- 2
    linkfunc[which(tolower(tmp)=="log")] <- 3
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

    for(g in seq_len(no.mg)) lambda[[g]] <- matrix(c(rep(1,K),rnorm(K,0,0.5)),ncol = 2)

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
  if(!is.matrix(saturated$prior)|!is.data.frame(saturated$prior)) saturated$prior <- matrix(saturated$prior,ncol = no.mg)

  if (any(att.prior < 0)){
    stop(
      "Joint attribute distribution prior can only contain numerical values between 0 and 1.",
      call. = FALSE
    )
  }

  logprior <- log(ColNormalize(att.prior))

  if (any(att.str)) control$smallNcorrection <- c(-1,-1)


  if(is.null(solver)){
    solver <- rep("auto",J)
  }else if (length(solver) == 1) {
    solver <- rep(solver, J)
  } else {
    if (length(solver) != J)
      stop("solver must have length of 1 or J.", call. = FALSE)
  }


  if (length(control$lower.p) == 1) {
    control$lower.p <- rep(control$lower.p, J)
  } else {
    if (length(control$lower.p) != J)
      stop("lower.p must have length of 1 or J.", call. = FALSE)
  }
  if (length(control$upper.p) == 1) {
    control$upper.p <- rep(control$upper.p, J)
  } else {
    if (length(control$upper.p) != J)
      stop("upper.p must have length of 1 or J.", call. = FALSE)
  }

  if(length(control$maxitr)==1L) {
    control$vmaxitr <- rep(control$maxitr,J)
  }else if(length(control$maxitr)!=J){
    warning("Length of maxitr must be equal to 1 or the number of nonzero categories.",call. = FALSE)
  }else{
    control$vmaxitr <- control$maxitr
    control$maxitr <- max(control$maxitr)
  }

  # if(is.null(weights)) weights <- rep(1,N)

  if (is.null(catprob.parm)) {
    item.parm <- initials(Q, control$nstarts, randomseed = control$randomseed,latent.var=latent.var)  #a list with nstarts matrices of size J x 2^Kjmax
    ###### Multiple starting values
    if (control$nstarts > 1L) {
      neg2LL <- vector("numeric",control$nstarts)
      for (i in seq_len(control$nstarts)) {
        neg2LL[i] <- -2*ObsLogLik(mpar = as.matrix(item.parm[[i]]),
                                  mX = as.matrix(dat),
                                  vlogPrior = as.matrix(logprior),
                                  vgroup = gr,
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
  delta <- calc_delta(item.parm, DesignMatrices = DesignMatrices, linkfunc = linkfunc)

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

  while(itr < control$maxitr)
  {
    pseudo.parm <- item.parm
    if(any(att.str)) pseudo.parm[is.nan(item.parm)] <- 0.5
    estep <- LikNR(mpar = as.matrix(pseudo.parm),
                   mX = as.matrix(dat),
                   vlogPrior = as.matrix(logprior),
                   vgroup = gr,
                   mloc = as.matrix(parloc),
                   weights = rep(1,N),
                   simplify = 1)



    # length of J indicating whether a nonzero cat should  be est. (TRUE) or not (FALSE)
    est.bin <- control$vmaxitr>itr
    # print(estep$Rg)
    optims <- Mstep(J=J, Kj = Kj, Rg = estep$Rg, Ng = estep$Ng, model = model, item.parm = item.parm, delta = delta,
                    ConstrType = ConstrType, correction = control$smallNcorrection, lower.p = control$lower.p, itr = itr + 1,
                    upper.p = control$upper.p, ConstrMatrix=ConstrMatrix,linkfunc=linkfunc,MstepMessage = control$MstepMessage,
                    solver = solver, DesignMatrices = DesignMatrices, est.bin = est.bin, item.prior = item.prior,
                    auglag_args = auglag_args,solnp_args = solnp_args,nloptr_args = nloptr_args)

    item.parm <- optims$item.parm
    delta <- optims$delta


    struc.parm <- structural.parm(AlphaPattern = AlphaPattern, no.mg = no.mg, logprior=estep$logprior,
                                  att.dist=att.dist,att.str=att.str, saturated = saturated,initial.logprior = initial.logprior,
                                  lower.prior = control$lower.prior,loglinear=loglinear,
                                  Ng=Ng,K=K,N=N,higher.order=higher.order,lambda = lambda)

    if(!is.null(higher.order)){
      higher.order <- struc.parm$higher.order
      ho.npar <- struc.parm$higher.order.npar
      }



    lambda <- struc.parm$lambda
    logprior <- struc.parm$logprior

    # logprior <- estep$logprior


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
  pseudo.parm <- item.parm
  if(any(att.str)) pseudo.parm[is.nan(item.parm)] <- 0.5

  estep <- LikNR(mpar = as.matrix(pseudo.parm),
                 mX = as.matrix(dat),
                 vlogPrior = as.matrix(logprior),
                 vgroup = gr,
                 mloc = as.matrix(parloc),
                 weights = rep(1,N),
                 simplify = 0)

  total.item.npar <- 0

  if(!att.str){
    total.item.npar <- sum(sapply(DesignMatrices,ncol))
  }else{
    total.item.npar <- sum(is.finite(c(item.parm)))
  }

  free.item.npar <- 0

  if(any(control$vmaxitr>0)){
    est.item <- which(control$vmaxitr>0)
    if(!att.str){
      free.item.npar <- sum(sapply(DesignMatrices[est.item],ncol))
    }else{
      free.item.npar <- sum(is.finite(c(item.parm[est.item,,drop=FALSE])))
    }
  }

  stru.npar <- 0

  if (any(att.dist=="higher.order")) {
    stru.npar <- ho.npar
  }
  for(g in 1:no.mg){
    if (att.dist[g]=="saturated") {
      if (!att.str) {
        stru.npar <- stru.npar + L - 1
      } else {
        stru.npar <- stru.npar + sum(is.finite(logprior[,g])) - 1
      }
    }else if(att.dist[g]=="loglinear"){
      stru.npar <- stru.npar + sum(sapply(seq_len(loglinear[g]),choose,n=K)) + 1
    }else if(att.dist[g]=="independent"){
      stru.npar <- stru.npar + K
    }
  }

  npar <- free.item.npar + stru.npar

  neg2LL <- -2 * estep$LL

  item.prob <- vector("list",J)
  initial.parm <- m2l(initial.parm)
  for (j in seq_len(J)){
    item.prob[[j]] <- item.parm[j,1:Lj[j]]
    names(initial.parm[[j]]) <- names(item.prob[[j]]) <- paste0("P(",apply(attributepattern(Kj[j]),1,paste0,collapse = ""),")")
  }
  postP <- exp(t(estep$logprior))
  pf <- LC.Prob <- uP(parloc,item.parm)
  if(sequential){ # transform LC.prob
    p <- LC.Prob
    for (j in 1:J){
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

  if(!is.null(group)) rownames(postP) <- paste("Group",gr.label)

  if(no.mg==1) {
    att.prior = c(exp(logprior0))
  }else{
    att.prior = exp(logprior0)
  }

  list(catprob.parm = item.prob, delta.parm = delta, catprob.matrix = item.parm,
       struc.parm = lambda, model = model.names, LC.prob = LC.Prob,
       posterior.prob = postP, pf = pf,
       testfit = list(Deviance=neg2LL,npar = npar,item.npar = free.item.npar, AIC=2 * npar + neg2LL, BIC=neg2LL + npar * log(N)),
       technicals = list(logposterior.i = estep$logpost, loglikelihood.i = estep$loglik, free.item.npar = free.item.npar,
                         total.item.npar = total.item.npar, stru.npar = stru.npar, total.npar = npar,
                         expectedCorrect = estep$Rg, expectedTotal = estep$Ng,initial.parm = initial.parm,
                         LC.labels = LC.labels),
       options = list(dat = originalData, Q = originalQ, Qm = Q, Qcm = Qcm, model = model,
                      itr = itr, dif.LL = dif.parm$neg2LL,dif.p=dif.parm$ip,dif.prior=dif.parm$prior,
                      att.dist=att.dist, higher.order=higher.order,att.prior = att.prior, latent.var = latent.var,
                      mono.constraint = mono.constraint, item.names = item.names,group = group, gr = gr,
                      att.str= att.str,  seq.dat = dat, no.group = no.mg, group.label = gr.label,
                      verbose = verbose, catprob.parm = catprob.parm,sequential = sequential,
                      nloptr_args = nloptr_args,auglag_args=auglag_args,solnp_args = solnp_args,
                      linkfunc = linkfunc,higher.order = higher.order, loglinear = loglinear, solver = solver,
                      DesignMatrices = DesignMatrices,ConstrPairs = ConstrPairs),
       control = control)

}
