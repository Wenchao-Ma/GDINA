#' Item fit statistics from the power-divergence family
#'
#' Calculate item fit statistics from the power-divergence family
#'
#' @param GDINA.obj Object containing a model fitted with the \code{GDINA::GDINA} function.
#' @param bootstrap Logical; whether parametric bootstrap should be used.
#' @param R Integer; number of replicates in the bootstrap procedure if used.
#' @param Stone Logical; whether Stone indices should be computed (only available if \code{bootstrap = TRUE}).
#' @param init.parm Logical; whether the estimated item parameters are used in the estimation of the bootstrap replications.
#' @param lambda Numeric; parameter for the power-divergence fit statistic.
#' @param p.adjust.method p-values can be adjusted for multiple comparisons at item level. This is conducted using \code{p.adjust} function in \pkg{stats},
#'  and therefore all adjustment methods supported by \code{p.adjust} can be used, including \code{"holm"},
#'  \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"} and \code{"BY"}. See \code{p.adjust}
#'  for more details. \code{"holm"} is the default.
#' @param person.sim Character; how to simulate attribute profiles in the bootstrap replications.
#' @param cores Integer; number of cores for parallelization during bootstrap.
#' @param digits Integer; number of decimal digits to report.
#' @param bound Numeric; minimum possible value for probabilities.
#' @param seed random seed.
#' @return an object of class \code{itemfitPD} consisting of several elements including:
#' \describe{
#' \item{p}{the proportion correct statistics, adjusted and unadjusted p values for each item}
#' \item{r}{the transformed correlations, adjusted and unadjusted p values for each item pair}
#' \item{logOR}{the log odds ratios, adjusted and unadjusted p values for each item pair}
#' \item{maxitemfit}{the maximum proportion correct, transformed correlation, and log-odds ratio for each item with associated item-level adjusted p-values}
#' }
#'#' @importFrom foreach %dopar% foreach
#' @author Pablo Najera
#'   Universidad Pontificia Comillas
#'   \email{pnajera@comillas.edu}
#'
#' @author Wenchao Ma
#'   University of Minnesota
#'   \email{wma@umn.edu}
#' @export
#' @references
#' Najera, P., Ma, W., Sorrel, M. A. and Abad, F. J. (2025). Assessing Item-Level Fit for the Sequential G-DINA Model.\emph{Behaviormetrika}.
#'
#' @examples
#' \dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#'
#' mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' mod1
#' PDfit <- itemfitPD(mod1)
#' PDfit
#'
#' dat <- sim21seqDINA$simdat
#' Q <- sim21seqDINA$simQ
#' sDINA <- GDINA(dat,Q,model="DINA",sequential = TRUE)
#' PDfit <- itemfitPD(sDINA)
#' PDfit
#' PDfit <- itemfitPD(sDINA, bootstrap = TRUE, Stone = TRUE, cores = 10)
#' PDfit
#'}

itemfitPD <- function(GDINA.obj, lambda = 2/3, bootstrap = FALSE, R = 1000, Stone = FALSE,
                      init.parm = FALSE, p.adjust.method = "holm", person.sim = "post", cores = 2,
                      digits = 4, bound = 1e-10, seed = 123456){



  #------------------
  # Argument control
  #------------------
  s1 <- Sys.time()
  if(Stone & !bootstrap){
    warning("Stone's method can only be applied along with bootstrapping. Stone's method has not been applied.")
    Stone <- FALSE
  }
  if(!person.sim %in% c("post", "EAP", "MAP", "MLE")){stop("person.sim must be 'post', 'EAP', 'MAP', or 'MLE'")}


  #------------------------------
  # Call seq.itemfit.R1 function
  #------------------------------

  ifit <- itemfitPD1(GDINA.obj, Stone, lambda, p.adjust.method, digits = digits, bound=bound)

  #-----------------------
  # Without bootstrapping
  #-----------------------

  if(!bootstrap){
    res <- list(X2 = round(ifit$X2, digits), G2 = round(ifit$G2, digits), PD = round(ifit$PD, digits),
                N = ifit$N, O = round(ifit$O, digits), E = round(ifit$E, digits))

  }else{
    # Gather information
    sequential <- extract(GDINA.obj,"sequential")
    model <- extract(GDINA.obj,"models")
    catprob.parm <- extract(GDINA.obj,"catprob.parm")
    dat <- extract(GDINA.obj,"dat")
    N <- extract(GDINA.obj,"nobs")
    K <- extract(GDINA.obj,"natt")
    Q <- extract(GDINA.obj,"originalQ")
    J <- extract(GDINA.obj,"nitem")
    H <- rep(1, J)
    if(sequential) H <- table(Q[,1])

    post <- extract(GDINA.obj, "posterior.prob")
    patt <- extract(GDINA.obj,"attributepattern")
    if(init.parm){
      init.catprob <- GDINA.obj$catprob.parm
    } else {
      init.catprob <- NULL
    }

    # Bootstrapping
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("doParallel is required for parallel execution.")
    }
    if (!requireNamespace("doRNG", quietly = TRUE)) {
      stop("doRNG is required for reproducibility in parallel execution.")
    }
    num.cores <- parallel::detectCores()
    if(cores > num.cores){stop(paste("The maximum number of CPU cores supported by this machine is", num.cores))}

    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG(seed)

    boot.res <- foreach::foreach(r = 1:R, .packages = c("GDINA"), .combine = rbind, .inorder = FALSE) %dopar% {
      try({
        if(person.sim == "post"){
          gen.att <- patt[sample(1:length(post), N, T, post),]
        } else {
          gen.att <- personparm(GDINA.obj, what = person.sim)[sample(1:N, N, T), 1:K]
        }
        sim <- simGDINA(N, Q, catprob.parm = catprob.parm, model = model, attribute = gen.att, sequential = sequential)
        boot.fit <- GDINA(sim$dat, Q, model = model, sequential = sequential, catprob.parm = init.catprob, verbose = 0)

        boot.ifit <- itemfitPD1(boot.fit, Stone, lambda, p.adjust.method, digits, bound)
        return(c(boot.ifit$X2$stat, boot.ifit$G2$stat, boot.ifit$PD$stat))
      })
    }
    parallel::stopCluster(cl)
    boot.X2 <- boot.res[, 1:J]
    boot.G2 <- boot.res[, (J + 1):(2 * J)]
    boot.PD <- boot.res[, (2 * J + 1):(3 * J)]
    rownames(boot.X2) <- rownames(boot.G2) <- rownames(boot.PD) <- paste0("R", 1:R)
    colnames(boot.X2) <- colnames(boot.G2) <- colnames(boot.PD) <- paste0("J", 1:J)

    # Calculate statistical significance
    X2 <- G2 <- PD <- data.frame(item = 1:J, stat = rep(NA, J), p = rep(NA, J), adjp = rep(NA, J))
    X2$stat <- ifit$X2$stat
    X2$p <- sapply(1:J, function(j) mean(ifit$X2$stat[j] < boot.X2[,j]))
    X2$adjp <- p.adjust(X2$p, method = p.adjust.method)
    G2$stat <- ifit$G2$stat
    G2$p <- sapply(1:J, function(j) mean(ifit$G2$stat[j] < boot.G2[,j]))
    G2$adjp <- p.adjust(G2$p, method = p.adjust.method)
    PD$stat <- ifit$PD$stat
    PD$p <- sapply(1:J, function(j) mean(ifit$PD$stat[j] < boot.PD[,j]))
    PD$adjp <- p.adjust(PD$p, method = p.adjust.method)

    # Return results
    res <- list(X2 = round(X2, digits), boot.X2 = round(boot.X2, digits),
                G2 = round(G2, digits), boot.G2 = round(boot.G2, digits),
                PD = round(PD, digits), boot.PD = round(boot.PD, digits),
                N = round(ifit$N, digits), O = round(ifit$O, digits), E = round(ifit$E, digits))
  }

  res$options=list(lambda = lambda,bootstrap = bootstrap, R = R, Stone = Stone,
                   p.adjust.method = p.adjust.method, person.sim = person.sim, cores = cores)
  s2 <- Sys.time()
  res$time = s2 - s1
  class(res) <- "itemfitPD"
  return(res)
}

# Auxiliary functions
itemfitPD1 <- function(GDINA.obj, Stone = FALSE, lambda = 2/3, p.adjust.method = "holm", digits = 4, bound = 1e-10){

  # Argument control
  if(!p.adjust.method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")){stop("p.adjust.method must be 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', or 'none'")}

  # Gather information
  sequential <- extract(GDINA.obj,"sequential")
  model <- extract(GDINA.obj,"models")
  catprob.parm <- extract(GDINA.obj,"catprob.parm")
  dat <- extract(GDINA.obj,"dat")
  N <- extract(GDINA.obj,"nobs")
  K <- extract(GDINA.obj,"natt")
  Q <- extract(GDINA.obj,"originalQ")
  J <- extract(GDINA.obj,"nitem")
  H <- rep(1, J)
  if(sequential) H <- table(Q[,1])

  post <- extract(GDINA.obj, "posterior.prob")
  patt <- apply(extract(GDINA.obj,"attributepattern"), 1, paste, collapse = "")
  att <- apply(personparm(GDINA.obj, "EAP"), 1, paste, collapse = "")
  p <- GDINA.obj$LC.prob

  if(Stone){
    Nl <- post * N
  } else {
    Nl <- table(att)[patt]
  }
  empty.l <- which(is.na(Nl))
  L <- length(Nl) - length(empty.l)
  Lj <- rep(L, J)
  if(length(empty.l) > 0){
    Nl <- Nl[-empty.l]
    patt <- patt[-empty.l]
  }
 # cat("L=",L,"J=",J,"H=",H)
  # Calculate N, O, and E, and fit indices
  O.ljh <- E.ljh <- array(0, dim = c(L, J, max(H) + 1), dimnames = list(paste0("L", 1:L), paste0("J", 1:J), paste0("H", 0:max(H))))
  X2 <- G2 <- PD <- data.frame(item = 1:J, stat = rep(NA, J), df = rep(NA, J), p = rep(NA, J), adjp = rep(NA, J))
  npar <- sapply(GDINA.obj$delta.parm, length)
  npar[npar > L] <- L
  if(!sequential){
    m <- as.numeric(npar)
  } else {
    m <- sapply(1:J, function(j) sum(npar[grepl(paste0("Item ", j, " "), names(npar))]))
  }
  if(Stone){
    r.ljh <- array(0, dim = c(L, J, max(H)))
    for(h in 1:max(H)){
      dat.h <- matrix(0, nrow(dat), ncol(dat))
      dat.h[dat == h] <- 1
      r.ljh[,,h] <- t(t(dat.h) %*% exp(GDINA.obj$technicals$logposterior.i))
    }
  }
  for(j in 1:J){
    x2 <- g2 <- pd <- 0
    for(l in 1:L){
      dat.lj <- dat[which(att == patt[l]), j]
      if(!sequential){
        p.lj <- p[j, patt[l]]
      } else {
        p.lj <- p[grepl(paste0("Item ", j, " "), rownames(p)), patt[l]]
      }
      p.lj <- c(1 - sum(p.lj), p.lj)
      if(Stone){
        o.lj <- r.ljh[l, j, 1:H[j]]
        o.lj <- c(Nl[l] - sum(o.lj), o.lj)
      } else {
        o.lj <- table(dat.lj)
        oh0.lj <- setdiff(0:H[j], as.numeric(names(o.lj)))
        if(length(oh0.lj) > 0){
          o.lj <- c(o.lj, rep(0, length(oh0.lj)))
          names(o.lj) <- c(names(table(dat.lj)), oh0.lj)
          o.lj <- o.lj[order(names(o.lj))]
        }
      }
      e.lj <- p.lj * Nl[l]
      o.lj[o.lj == 0] <- bound; o.lj[is.na(o.lj)] <- bound
      e.lj[e.lj == 0] <- bound; e.lj[is.na(e.lj)] <- bound
      O.ljh[l, j, 1:(H[j] + 1)] <- o.lj
      E.ljh[l, j, 1:(H[j] + 1)] <- e.lj
    }
    O <- tmp.O <- O.ljh[, j, 1:(H[j] + 1)]
    E <- tmp.E <- E.ljh[, j, 1:(H[j] + 1)]
    Lj[j] <- nrow(O)
    X2$stat[j] <- sum((O - E)^2 / E)
    G2$stat[j] <- 2 * sum(O * log(O / E))
    PD$stat[j] <- 2 / (lambda * (lambda + 1)) * sum(O * ((O / E)^lambda - 1))
  }

  # Calculate statistical significance
  X2$df <- G2$df <- PD$df <- Lj * H - m
  X2$p <- 1 - pchisq(X2$stat, X2$df)
  X2$adjp <- p.adjust(X2$p, method = p.adjust.method)
  G2$p <- 1 - pchisq(G2$stat, G2$df)
  G2$adjp <- p.adjust(G2$p, method = p.adjust.method)
  PD$p <- 1 - pchisq(PD$stat, PD$df)
  PD$adjp <- p.adjust(PD$p, method = p.adjust.method)

  # Return results
  return(list(X2 = round(X2, digits), G2 = round(G2, digits), PD = round(PD, digits),
              N = round(Nl, digits), O = round(O.ljh, digits), E = round(E.ljh, digits)))
}



