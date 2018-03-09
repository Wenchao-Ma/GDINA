#' Model fit statistics
#'
#' Calculate various model-data fit statistics
#'
#' Various model-data fit statistics including M2 statistic for dichotmous CDMs and Mord or M2* statistic for graded responses.
#' It also calculates SRMSR and RMSEA2.
#'
#' @param GDINA.obj An estimated model object of class \code{GDINA}
#' @param CI numeric value from 0 to 1 indicating the range of the confidence interval for RMSEA. Default returns the 90\% interval.
#' @param ... arguments passed to the function
#'
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @export
#' @references
#'
#' Maydeu-Olivares, A. (2013). Goodness-of-Fit Assessment of Item Response Theory Models. \emph{Measurement, 11}, 71-101.
#'
#' Liu, Y., Tian, W., & Xin, T. (2016). An Application of M2 Statistic to Evaluate the Fit of Cognitive Diagnostic Models. \emph{Journal of Educational and Behavioral Statistics, 41}, 3-26.
#' @examples
#' \dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' mod1 <- GDINA(dat = dat, Q = Q, model = "DINA")
#' modelfit(mod1)
#'}


modelfit <- function(GDINA.obj,CI=0.90,...)
{

  ItemOnly <- dots("ItemOnly",FALSE,...)
  if(extract(GDINA.obj,"sequential")) stop("modelfit is not available for sequential models.",call. = FALSE)

  if(extract(GDINA.obj, "ngroup")!=1) {
    stop("modelfit is only applicable to single group analysis.", call. = FALSE)
  }
  if (any(extract(GDINA.obj, "models_numeric") < 0) ||
      any(extract(GDINA.obj, "models_numeric") > 6))
    stop("modelfit is only applicable to GDINA, DINA, DINO, ACDM, LLM and RRUM.",
         call. = FALSE)
  if (extract(GDINA.obj, "att.dist") %in% c("higher.order","independent","fixed")){
    stop(paste("modelfit is not available for ",extract(GDINA.obj, "att.dist"),"joint attribute distribution."),call. = FALSE)
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
  L <- 2^K
  Kj <- extract(GDINA.obj, "Kj")
  nparJ <- npar(GDINA.obj)$`No. of item parameters`
  models <- extract(GDINA.obj, "models")
  post <- c(extract(GDINA.obj, "posterior.prob"))
  pj <- extract(GDINA.obj, "LCprob.parm") # P(Xi=s|alpha) S x L without category 0
  pf <- extract(GDINA.obj, "LCpf.parm") # S(Xi=s|alpha) S x L without category 0

  # Observed p
  crossp <-
    crossprod.na(dat, dat, val = 0) / crossprod(!is.na(dat), !is.na(dat))
  # univariate and bivariate observed
  p <-
    c(colMeans(dat, na.rm = TRUE), crossp[lower.tri(crossp)]) # length of nitems + nitems*(nitems-1)/2

  Xi <- Mord(item.no, as.matrix(pj), post)
  Xi2 <- cbind(rbind(Xi$Xi11, Xi$Xi21), rbind(t(Xi$Xi21), Xi$Xi22))

  e <-
    c(Xi$uni, Xi$bi[lower.tri(Xi$bi)])   # length of nitems + nitems*(nitems-1)/2
  se <- sqrt(diag(Xi$bi) - c(Xi$uni) ^ 2)
  difr <-
    cor(dat, use = "pairwise.complete.obs") - (Xi$bi - Xi$uni %*% t(Xi$uni)) /
    (se %*% t(se))
  SRMSR <-
    sqrt(sum((difr[lower.tri(difr)]) ^ 2 / (nitems * (nitems - 1) / 2)))
  ##### SRMSR
  #  Maydeu-Olivares, A. (2013). Goodness-of-Fit Assessment of Item Response Theory Models. Measurement, 11(3), 71–101.
  #  https://doi.org/10.1080/15366367.2013.831680
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
    patt <- eta(as.matrix(Q))
    Mj <- list()
    for (s in 1:ncat)
      Mj[[s]] <- designmatrix(Kj[s], models[s])

    parloc <- matrix(c(cumsum(c(
      1, unlist(lapply(Mj, ncol))[-ncat]
    )),
    cumsum(unlist(lapply(
      Mj, ncol
    )))), ncol = 2)
    # partial s/ partial d
    extMj <- vector("list", ncat)
    for (j in 1:ncat) {
      extMj[[j]] <- Mj[[j]][patt[j, ],]
      if (models[j] == "LLM") {
        extMj[[j]] <- pj[j, ] * (1 - pj[j, ]) * extMj[[j]]
      } else if (models[j] == "RRUM") {
        extMj[[j]] <- pj[j, ] * extMj[[j]]
      }

    }
    # seq_component is partial E/partial s
    seq_component <- matrix(0, nrow(Qc), L)
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
    pcpl <- rbind(diag(L - 1), -1)
    ##=============== need to double check
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

    if (ItemOnly) {
      delt <- rbind(delta11, delta21)
    } else{
      if (extract(GDINA.obj, "att.dist") == "saturated") {
        delta12 <- expected.score %*% rbind(diag(L - 1), -1)
        delta22 <- delta22E %*% rbind(diag(L - 1), -1)
      } else if (extract(GDINA.obj, "att.dist") == "loglinear") {
        Z <- designM(K, 0)
        loglinear <- extract(GDINA.obj, "loglinear")
        if (loglinear == 1) {
          Z <- Z[, seq_len(1 + K)]
        } else if (loglinear == 2) {
          Z <- Z[, seq_len(1 + K * (K + 1) / 2)]
        } else if (loglinear == 3) {
          Z <- Z[, seq_len(1 + choose(K, 1) + choose(K, 2) + choose(K, 3))]
        } else{
          stop("Argument 'loglinear' must be 1, 2 or 3.", call. = FALSE)
        }
        delta12 <- expected.score %*% (post * Z)
        delta22 <- delta22E %*% (post * Z)
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

    ci <- RMSEA.CI(M2, df, N)
  }


  out <- c(list(M2=c(M2),M2.pvalue=c(sig),M2.df=df,SRMSR = SRMSR,RMSEA=rmsea,RMSEA.CI = ci,CI=CI),GDINA.obj$testfit)


  class(out) <- "modelfit"
 return(out)
}


RMSEA.CI <- function(X2, df, N, ci.lower=.05, ci.upper=.95) {

  ##########################################################################################################################
  # RMSEA.CI
  # This function is copied from mirt package with minor modifications
  # R. Philip Chalmers (2012). mirt: A Multidimensional Item Response Theory Package for the R Environment. Journal of Statistical
  # Software, 48(6), 1-29. doi:10.18637/jss.v048.i06
  ##########################################################################################################################


  lower.lambda <- function(lambda) pchisq(X2, df=df, ncp=lambda) - ci.upper

  upper.lambda <- function(lambda) pchisq(X2, df=df, ncp=lambda) - ci.lower



  lambda.l <- try(uniroot(f=lower.lambda, lower=0, upper=X2)$root, silent=TRUE)

  lambda.u <- try(uniroot(f=upper.lambda, lower=0, upper=max(N, X2*5))$root, silent=TRUE)

  if(!inherits(lambda.l, 'try-error')){

    RMSEA.lower <- sqrt(lambda.l/(N*df))

  } else {

    RMSEA.lower <- 0

  }

  if(!inherits(lambda.u, 'try-error')){

    RMSEA.upper <- sqrt(lambda.u/(N*df))

  } else {

    RMSEA.upper <- 0

  }



  return(c(RMSEA.lower, RMSEA.upper))

}
