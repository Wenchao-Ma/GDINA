

internalextract <- function(object, what, ...) {
  out <- switch(
    what,
    catprob.parm = object$catprob.parm,
    catprob.matrix = object$catprob.matrix,
    itemprob.parm = {
      if(internalextract(object,"sequential")){
        itemparm <- internalextract(object,"catprob.matrix")
        p <- vector("list",internalextract(object,"nitem"))
        Q <- internalextract(object,"Qc")
        for (j in 1:internalextract(object,"nitem")){
          Qj <- Q[which(Q[,1]==j),-c(1:2),drop=FALSE]
          itemparj <- itemparm[which(Q[,1]==j),,drop=FALSE]
          colj <- which(colSums(Qj)>0)
          Kj <- length(colj)
          if(nrow(Qj)>1){ # polytomous items
            if(length(colj)>1){
              redQj <- Qj[,colj]
              pj <- sj <- rbind(uP(as.matrix(eta.loc(redQj)), as.matrix(itemparj)),0)
            }else{
              pj <- sj <- rbind(itemparj[,1:2],0)
            }
            for (s in 1:(nrow(sj)-1)){
              pj[s,] <- apply(sj[1:s,,drop=FALSE],2,prod)*(1-sj[s+1,])
            }
            pj <- pj[-nrow(pj),]
          }else{ #dichotomous items
            pj <- matrix(itemparj[!is.na(itemparj)],nrow = 1)
          }
          colnames(pj) <- paste0("P(", apply(alpha(Kj), 1, paste0, collapse = ""), ")")
          rownames(pj) <- paste("Cat",1:nrow(pj))
          p[[j]] <- pj
        }
        names(p) <- paste("Item",1:internalextract(object,"nitem"))
        p
      }else{
        object$itemprob.parm
      }

    },
    catprob.se = {
      if (length(list(...)) == 0)
        stop("Please specify variance-covariance estimation method using type = 1 or 2.",
             call. = FALSE)

      Kj <- internalextract(object, "Kj")
      m <- internalextract(object, "models_numeric")
      if (all(m < 3)) {
        # G-DINA, DINO or DINA only
        se <- itemprob_se(object, ...)$se
        for (j in 1:length(se))
          names(se[[j]]) <-
            paste0("SE[P(", apply(alpha(Kj[j]), 1, paste0, collapse = ""), ")]")
      } else{
        v <- internalextract(object, "catprob.cov", ...)
        sej <- sqrt(diag(v$cov))
        ind <- v$index
        se <- vector("list", length(unique(ind$item)))

        for (j in sort(unique(ind$item))) {
          se[[j]] <- sej[which(ind$item == j)]
          names(se[[j]]) <-
            paste0("SE[P(", apply(alpha(Kj[j]), 1, paste0, collapse = ""), ")]")
        }
      }
      se
    },
    itemprob.se = {
      if(internalextract(object,"sequential")){
        NULL
      }else{
        internalextract(object,"catprob.se",...)
      }
    },
    delta.parm = object$delta.parm,
    delta.se = delta_se(object, ...)$se,
    higher.order.struc.parm = object$higher.order.struc.parm,
    logLik = -0.5 * object$testfit$Deviance,
    deviance = object$testfit$Deviance,
    npar = object$options$npar,
    npar.item = object$options$item.npar,
    npar.att = object$options$npar - object$options$item.npar,
    nitr = object$options$itr,
    nobs = nrow(object$options$dat),
    nitem = ncol(object$options$dat),
    ncat = nrow(object$options$Q),
    natt = ifelse(
      object$options$sequential,
      ncol(object$options$Q) - 2,
      ncol(object$options$Q)
    ),
    AIC = object$testfit$AIC,
    BIC = object$testfit$BIC,
    models = object$model,
    models_numeric = object$options$model,
    Kj = rowSums(internalextract(object, "Q")),
    LC.prob = object$LC.prob,
    discrim = {

      wp <- t(object$LC.prob) * c(object$posterior.prob)  #L x J w*p matrix
      gdi <-
        colSums(t((object$LC.prob - colSums(wp)) ^ 2) * c(object$posterior.prob)) # vector of length J
      dj <-
        sapply(object$catprob.parm, function(x)
          x[length(x)] - x[1])
      Discrim <- data.frame(dj = dj, GDI = gdi)
      rownames(Discrim) <- object$options$item.names
      colnames(Discrim) <- c("P(1)-P(0)", "GDI")
      Discrim
    },
    time = object$options$timeused,
    start.time = object$options$start.time,
    end.time = object$options$end.time,
    dat = object$options$dat,
    seq.dat = object$options$seq.dat,
    Qc = {if (object$options$sequential) {
      data.frame(object$options$Q)
    } else{
      out <- data.frame(
        Item = 1:nrow(object$options$Q),
        Cat = rep(1, nrow(object$options$Q)),
        object$options$Q
      )
      colnames(out)[-c(1:2)] <-
        paste0("A", 1:ncol(object$options$Q))
      out
    }},
    Q = {
      if (object$options$sequential) {
        out <- data.frame(object$options$Q[, -c(1:2)])
      } else{
        out <- data.frame(object$options$Q)
      }
      colnames(out) <- paste0("A", 1:ncol(out))
      out
    },
    posterior.prob = object$posterior.prob,
    prevalence = {
      if (!object$options$sequential) {
        tmp <- object$options$Q
      } else{
        tmp <- object$options$Q[, -c(1:2)]
      }
      preva <- NULL
      pattern <- alpha(ncol(tmp), T, tmp)
      for (i in c(0:max(tmp))) {
        preva <-
          cbind(preva, t(object$posterior.prob %*% ((pattern == i) * 1)))
      }
      colnames(preva) <-
        paste("Level", 0:max(tmp), sep = "")
      rownames(preva) <-
        paste("A", seq(1, ncol(tmp)), sep = "")
      preva
    },
    catprob.cov = {
      if (length(list(...)) == 0)
        stop("Please specify variance-covariance estimation method using type = 1 or 2.",
             call. = FALSE)

      m <- internalextract(object, "models_numeric")
      if (all(m < 3)) {
        # G-DINA, DINO or DINA only
        var <- itemprob_se(object, ...)
      } else{
        J <- internalextract(object, "ncat")
        Kj <- internalextract(object, "Kj")
        var <- vector("list", J)
        delta = delta_se(object, ...)
        for (j in 1:J) {
          if (m[j] <= 3) {
            var[[j]] <-
              designmatrix(Kj[j], m[j]) %*% delta$cov[delta$ind[delta$ind$item == j, "loc"],
                                                      delta$ind[delta$ind$item ==
                                                                  j, "loc"]] %*% t(designmatrix(Kj[j], m[j]))
          } else if (m[j] == 4) {
            grad <-
              diag(internalextract(object, "catprob.parm")[[j]] * (1 - internalextract(object, "catprob.parm")[[j]])) %*%
              designmatrix(Kj[j], m[j])
            var[[j]] <-
              grad %*% delta$cov[delta$ind[delta$ind$item == j, "loc"],
                                 delta$ind[delta$ind$item ==
                                             j, "loc"]] %*% t(grad)

          } else if (m[j] == 5) {
            grad <-
              diag(internalextract(object, "catprob.parm")[[j]]) %*% designmatrix(Kj[j], m[j])
            var[[j]] <-
              grad %*% delta$cov[delta$ind[delta$ind$item == j, "loc"],
                                 delta$ind[delta$ind$item ==
                                             j, "loc"]] %*% t(grad)

          }

        }
        var <- list(cov = bdiag(var),
                    ind = data.frame(item = rep(1:J, unlist(
                      lapply(var, nrow)
                    )),
                    loc = seq(1, sum(
                      unlist(lapply(var, nrow))
                    ))))
      }


      list(cov = var$cov, index = var$ind)
    },
    delta.cov = {
      if (length(list(...)) == 0)
        stop("Please specify variance-covariance estimation method using type = 1 or 2.",
             call. = FALSE)
      var <- delta_se(object, ...)
      list(cov = var$cov, index = var$ind)
    },
    logposterior.i = object$technicals$logposterior.i,
    loglikelihood.i = object$technicals$loglikelihood.i,
    expectedCorrect = object$technicals$expectedCorrect,
    expectedTotal = object$technicals$expectedTotal,
    higher.order = object$options$higher.order,
    higher.order.model = object$options$higher.order.model,
    mono.constraint = object$options$mono.constraint ,
    SE = object$options$SE,
    SE.type = object$options$SE.type,
    empirical = object$options$empirical,
    att.prior = object$options$att.prior,
    att.str = object$options$att.str,
    nstarts = object$options$nstarts,
    conv.crit = object$options$conv.crit,
    maxitr = object$options$maxitr,
    higher.order.method = object$options$higher.order.method,
    higher.order.SE = object$options$higher.order.SE,
    verbose = object$options$verbose,
    sequential = object$options$sequential,
    higher.order.parm = object$options$higher.order.parm,
    digits = object$options$digits,
    dif.LL = object$options$dif.LL,
    dif.p = object$options$dif.p,
    item.names = object$options$item.names,
    itemprob.history = object$diagnos$itemprob.matrix,
    RN.history = object$diagnos$RN,
    likepost.history = object$diagnos$likepost,
    iter.history = object$diagnos$changelog,
    HO.parm.history = object$diagnos$HO.parm,
    call = object$options$call,
    expectedCorrect.LC = {
      if (object$options$sequential) {
        out <- apply(seq_coding(object$options$dat, object$options$Q), 2,
                     function(x)
                       colSums(x * exp(object$technicals$logposterior.i), na.rm = TRUE))
      } else{
        out <- apply(object$options$dat, 2,
                     function(x)
                       colSums(x * exp(object$technicals$logposterior.i), na.rm = TRUE))
      }
      out <- t(out)
      row.names(out) <- object$options$item.names
      out
    },
    expectedTotal.LC = {
      if (object$options$sequential) {
        out <- apply(seq_coding(object$options$dat, object$options$Q), 2,
                     function(x)
                       colSums(
                         as.numeric(!is.na(x)) * exp(object$technicals$logposterior.i),
                         na.rm = TRUE
                       ))
      } else{
        out <-
          matrix(
            colSums(exp(object$technicals$logposterior.i), na.rm = TRUE),
            ncol = nrow(object$options$Q),
            nrow = ncol(object$technicals$logposterior.i)
          )
      }
      out <- t(out)
      row.names(out) <- object$options$item.names
      colnames(out) <-
        colnames(object$technicals$logposterior.i)
      out
    },
    stop(sprintf("Can not extract element \'%s\'", what), call. =
           FALSE)
  )
  return(out)
}
