#' Model fit statistics
#'
#' Calculate various absolute model-data fit statistics
#'
#' Various model-data fit statistics including M2 statistic for G-DINA model with dichotmous responses (Liu, Tian, & Xin, 2016; Hansen, Cai, Monroe, & Li, 2016) and for sequential G-DINA model with graded responses (Ma, 2020).
#' It also calculates SRMSR and RMSEA2.
#'
#' @param GDINA.obj An estimated model object of class \code{GDINA}
#' @param CI numeric value from 0 to 1 indicating the range of the confidence interval for RMSEA. Default returns the 90\% interval.
#' @param ItemOnly should joint attribute distribution parameters be considered? Default = FALSE. See Ma (2019).
#'
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu}}
#' @export
#' @references
#'
#' Hansen, M., Cai, L.,  Monroe, S., & Li, Z. (2016). Limited-information goodness-of-fit testing of diagnostic classification item response models. \emph{British Journal of Mathematical and Statistical Psychology. 69,} 225--252.
#'
#' Liu, Y., Tian, W., & Xin, T. (2016). An Application of M2 Statistic to Evaluate the Fit of Cognitive Diagnostic Models. \emph{Journal of Educational and Behavioral Statistics, 41}, 3-26.
#'
#' Ma, W. (2020). Evaluating the fit of sequential G-DINA model using limited-information measures. \emph{Applied Psychological Measurement, 44}, 167-181.
#'
#' Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive Diagnosis Modeling. \emph{Journal of Statistical Software, 93(14)}, 1-26.
#'
#' Maydeu-Olivares, A. (2013). Goodness-of-Fit Assessment of Item Response Theory Models. \emph{Measurement, 11}, 71-101.
#'
#' @examples
#' \dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' mod1 <- GDINA(dat = dat, Q = Q, model = "DINA")
#' modelfit(mod1)
#'}


modelfit <- function(GDINA.obj, CI = 0.90, ItemOnly = FALSE)
{

  if(CI>=1||CI<=0) stop("CI must be between 0 and 1.",call. = FALSE)

  if(extract(GDINA.obj, "ngroup")!=1) {
    stop("modelfit is only applicable to single group analysis.", call. = FALSE)
  }


  delta <- extract(GDINA.obj, "delta.parm")
  Q <- extract(GDINA.obj, "Q")
  if (max(Q) > 1) {
    stop("modelfit is only available for dichotomous attribute models.",
         call. = FALSE)
  }
  Qc <- extract(GDINA.obj, "Qc")
  item.no <- c(Qc[, 1])
  dat <- as.matrix(extract(GDINA.obj, "dat"))
  N <- nrow(dat)
  nitems <- ncol(dat)   # number of items
  ncat <- extract(GDINA.obj, "ncat")
  K <- extract(GDINA.obj, "natt")
  att <- as.matrix(extract(GDINA.obj,"attributepattern"))
  L <- extract(GDINA.obj, "nLC")
  Kj <- extract(GDINA.obj, "Kj")
  nparJ <- npar(GDINA.obj)$`No. of total item parameters`
  models <- extract(GDINA.obj, "models")
  post <- c(extract(GDINA.obj, "posterior.prob"))
  pj <- extract(GDINA.obj, "LCprob.parm") # P(Xi=s|alpha) S x L without category 0
  pf <- extract(GDINA.obj, "LCpf.parm") # S(Xi=s|alpha) S x L without category 0

  # Observed p
  crossp <- crossprod.na(dat, dat, val = 0) / crossprod(!is.na(dat), !is.na(dat))
  # univariate and bivariate observed
  # length of nitems + nitems*(nitems-1)/2
  p <- c(colMeans(dat, na.rm = TRUE), crossp[lower.tri(crossp)])

  Xi <- Mord(item.no, as.matrix(pj), post)
  Xi2 <- cbind(rbind(Xi$Xi11, Xi$Xi21), rbind(t(Xi$Xi21), Xi$Xi22))

  e <- c(Xi$uni, Xi$bi[lower.tri(Xi$bi)])   # length of nitems + nitems*(nitems-1)/2
  se <- sqrt(diag(Xi$bi) - c(Xi$uni) ^ 2)
  difr <- cor(dat, use = "pairwise.complete.obs") - (Xi$bi - Xi$uni %*% t(Xi$uni)) /
    (se %*% t(se))
  SRMSR <- sqrt(sum((difr[lower.tri(difr)]) ^ 2 / (nitems * (nitems - 1) / 2)))
  M2 <- sig <- df <- rmsea <- ci <- NULL
  if (ItemOnly &
      nitems * (nitems + 1) / 2 - extract(GDINA.obj, "npar.item") < 0) {
    warning("M2 statistic cannot be calculated - Degrees of freedom are too low.",
            call. = FALSE)
  } else if (!ItemOnly &
             nitems * (nitems + 1) / 2 - extract(GDINA.obj, "npar") < 0) {
    warning("M2 statistic cannot be calculated - Degrees of freedom are too low.",
            call. = FALSE)
  } else{
    # parameter locations

    patt <- extract(GDINA.obj,"eta")
    Mj <- extract(GDINA.obj,"designmatrix")
    linkf <- extract(GDINA.obj,"linkfunc")
    parloc <- matrix(c(cumsum(c(
      1, unlist(lapply(Mj, ncol))[-ncat]
    )),
    cumsum(unlist(lapply(
      Mj, ncol
    )))), ncol = 2)
    # partial s/ partial d
    extMj <- vector("list", ncat)
    for (j in 1:ncat) {
      if(linkf[j]=="identity"){
        extMj[[j]] <- Mj[[j]][patt[j, ],]
      }else if (linkf[j]=="logit") {
        extMj[[j]] <- pj[j, ] * (1 - pj[j, ]) * Mj[[j]][patt[j, ],]
      } else if (linkf[j]=="log") {
        extMj[[j]] <- pj[j, ] * Mj[[j]][patt[j, ],]
      }

    }
    # seq_component is partial E/partial s
    seq_component <- matrix(0, ncat, L)
    expected.score <- matrix(0, nitems, L)
    for (j in 1:nitems) {
      locj <- which(item.no == j)
      if (length(locj) > 1) {
        # sequential model
        pjs <- pj[locj, , drop = FALSE]
        pfs <- pf[locj, , drop = FALSE]
        expected.score[j, ] <- colSums(pjs * seq_len(nrow(pjs)))

        cumpf <- apply(pfs, 2, cumprod)
        cumpf <- rbind(0, cumpf[-nrow(cumpf), , drop = FALSE]) # sj x L
        cumpf <- apply(cumpf, 2, cumsum)
        seq_component[locj, ] <-
          (rep(1, nrow(pfs)) %o% expected.score[j, ] - cumpf) / pfs
      } else{
        # dichotomous item
        seq_component[locj, ] <- rep(1, L)
        expected.score[j, ] <- pj[locj, ]
      }
    }
    delta11 <- matrix(0, nitems, nparJ)
    delta21 <- matrix(0, length(p) - nitems, nparJ)

    delta22E <- matrix(0, nitems * (nitems - 1) / 2, L)
    loc <- 1
    for (j in 1:nitems) {
      locj <- which(item.no == j)
      for (sj in locj) {
        delta11[j, parloc[sj, 1]:parloc[sj, 2]] <-
          colSums(extMj[[sj]] * post * seq_component[sj, ])
      }

      for (i in 1:nitems) {
        if (j < i) {
          loci <- which(item.no == i)
          for (sj in locj) {
            delta21[loc, parloc[sj, 1]:parloc[sj, 2]] <-
              colSums(extMj[[sj]] * post * seq_component[sj, ] * expected.score[i, ])
          }
          for (si in loci) {
            delta21[loc, parloc[si, 1]:parloc[si, 2]] <-
              colSums(extMj[[si]] * post * seq_component[si, ] * expected.score[j, ])
          }
          delta22E[loc, ] <- expected.score[j, ] * expected.score[i, ]
          loc <- loc + 1
        }
      }
    }

    if (ItemOnly || extract(GDINA.obj, "att.dist")=="fixed") {
      delt <- rbind(delta11, delta21)
    } else{
      if (extract(GDINA.obj, "att.dist") == "saturated") {
        pcpl <- rbind(diag(L - 1), -1)
        delta12 <- expected.score %*% pcpl
        delta22 <- delta22E %*% pcpl
      } else if (extract(GDINA.obj, "att.dist") == "loglinear"){
        Z <- designM(K, 0)
        loglinear <- extract(GDINA.obj, "loglinear")
        Z <- Z[, seq_len(1 + sum(sapply(seq_len(extract(GDINA.obj,"loglinear")),choose,n=K)))]

        delta12 <- expected.score %*% (post * Z)
        delta22 <- delta22E %*% (post * Z)
      }else if(extract(GDINA.obj, "att.dist") == "higher.order"){

        higher.order <- extract(GDINA.obj,"higher.order")
        HOpar <- extract(GDINA.obj,"struc.parm")

        P.att.theta <- exp(logLikPattern(att,
                      higher.order$QuadNodes,
                      HOpar[,1],HOpar[,2]) +
              matrix(1,nrow(att),1)%*%t(log(higher.order$QuadWghts))) # 2^K x nnodes

        Pk_theta <- Pr_2PL_vec(higher.order$QuadNodes, HOpar[,1],HOpar[,2]) #P(\alpha_k|theta) nnodes x K
        if(higher.order$model=="Rasch"){
          delta12 <- matrix(0,nitems,K)
          delta22 <- matrix(0,nrow(delta22E),K)
          for(k in seq_len(nrow(HOpar))){

            # partial E/partial intercept
            delta12[,k] <- rowSums(expected.score%*%(P.att.theta * outer(att[,k],Pk_theta[,k],"-")))
            delta22[,k] <- rowSums(delta22E%*%(P.att.theta * outer(att[,k],Pk_theta[,k],"-")))

          }
        }else if(higher.order$model=="2PL"){
          delta12 <- matrix(0,nitems,2*K)
          delta22 <- matrix(0,nrow(delta22E),2*K)
          for(k in seq_len(nrow(HOpar))){

            # partial E/partial intercept
            delta12[,k] <- rowSums(expected.score%*%(P.att.theta * outer(att[,k],Pk_theta[,k],"-")))
            delta22[,k] <- rowSums(delta22E%*%(P.att.theta * outer(att[,k],Pk_theta[,k],"-")))
            # partial E/partial slope <= 2PL
            delta12[,k+K] <- c(expected.score%*%(P.att.theta * outer(att[,k],Pk_theta[,k],"-"))%*%higher.order$QuadNodes)
            delta22[,k+K] <- c(delta22E%*%(P.att.theta * outer(att[,k],Pk_theta[,k],"-"))%*%higher.order$QuadNodes)

          }
        }else if(higher.order$model=="1PL"){
          delta12 <- matrix(0,nitems,1+K)
          delta22 <- matrix(0,nrow(delta22E),1+K)
          for(k in seq_len(nrow(HOpar))){

            # partial E/partial intercept
            delta12[,k] <- rowSums(expected.score%*%(P.att.theta * outer(att[,k],Pk_theta[,k],"-")))
            delta22[,k] <- rowSums(delta22E%*%(P.att.theta * outer(att[,k],Pk_theta[,k],"-")))

          }
          # partial E/partial slope <= 2PL
          delta12[,1+K] <- c(expected.score%*%(P.att.theta * outer(rowSums(att),rowSums(Pk_theta),"-"))%*%higher.order$QuadNodes)
          delta22[,1+K] <- c(delta22E%*%(P.att.theta * outer(rowSums(att),rowSums(Pk_theta),"-"))%*%higher.order$QuadNodes)


        }
      }else if(extract(GDINA.obj, "att.dist") == "independent"){
        pr <- extract(GDINA.obj,"struc.parm")
        pr[pr < .Machine$double.eps] <- .Machine$double.eps
        pr[pr > 1 - .Machine$double.eps] <- 1 - .Machine$double.eps
        mp <- matrix(pr,nrow(att),ncol(att),byrow = TRUE)
        delta12 <- expected.score %*% diag(post) %*% ( (att - mp) / (mp * (1-mp)) )
        delta22 <- delta22E %*% diag(post) %*% ( (att - mp) / (mp * (1-mp)) )
      }
      delt <- cbind(rbind(delta11, delta21), rbind(delta12, delta22))
    }

    tmp <- qr.Q(qr(delt), complete = TRUE)
    deltac <- tmp[, (ncol(delt) + 1L):ncol(tmp), drop = FALSE]
    C2 <- deltac %*% solve(t(deltac) %*% Xi2 %*% deltac) %*% t(deltac)
    M2 <- N * t(p - e) %*% C2 %*% (p - e)
    df <- nrow(delt) - ncol(delt)
    sig <- 1 - pchisq(M2, df)

    rmsea <- sqrt(max((M2 - df) / (N * df), 0))

    ci <- RMSEAfun(M2, df, N, CI)
  }


  out <- c(list(M2=c(M2),M2.pvalue=c(sig),M2.df=df,SRMSR = SRMSR,RMSEA2=rmsea,RMSEA2.CI = ci,CI=CI),
           GDINA.obj$testfit,sequential=GDINA.obj$options$sequential)


  class(out) <- "modelfit"
 return(out)
}

func.lambda <- function(lambda, X2, df, b) pchisq(X2, df=df, ncp=lambda) - (1 - b)

RMSEAfun <- function(X2, df, N, CI) {

  bb <- (1 - CI)/2

  lower <- 0
  l.upper <- X2
  u.upper <- max(N, X2*5)
  if(func.lambda(lower,X2,df,b = bb)*func.lambda(l.upper,X2,df,b=bb)>0){
    l.lambda <- 0
  }else{
    l.lambda <- uniroot(f=func.lambda, lower=lower, upper=l.upper, X2=X2, b=bb, df=df)$root
  }


  if(func.lambda(lower,X2,df,b=CI + bb)*func.lambda(u.upper,X2,df,b = CI + bb)>0){
    u.lambda <- 0
  }else{
    u.lambda <- uniroot(f=func.lambda, lower=lower, upper=u.upper, X2=X2, b=CI + bb, df=df)$root
  }


  return(c(sqrt(l.lambda/(N*df)), sqrt(u.lambda/(N*df))))

}
