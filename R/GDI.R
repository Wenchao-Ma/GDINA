#' Q-matrix validation
#'
#' Q-matrix validation for the G-DINA model based on de la Torre and Chiu (2016).
#'
#' @param GDINA.obj An estimated model object of class \code{GDINA}
#' @param eps cutoff value for PVAF. 0.95 is the default.
#' @param method which Q-matrix validation method is used?
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
#' @seealso \code{\link{GDINA}}, \code{\link{mesaplot}}
#' @export
#' @examples
#'\dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' Q[10,] <- c(0,1,0)
#' mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' out <- Qval(mod1,eps = 0.95)
#' out
#' #If many entries are modified, you may want to check
#' #the PVAF plot using the function plotPVAF or
#' #to change eps. eps = 0.9 or 0.8 seems another two
#' #reasonable choices.
#' extract(out,what = "PVAF")
#' #See also:
#' extract(out,what = "varsigma")
#' extract(out,what = "sug.Q")
#'
#' # Draw a mesa plot
#' mesaplot(out,item=10,type="all",no.qvector=5)
#'}




Qval <- function(GDINA.obj,
                 method = "PVAF",
                 eps = 0.95,
                 digits = 4) {
  if (class(GDINA.obj) != "GDINA")
    stop("GDINA.obj must be a GDINA object from GDINA function.", call. = FALSE)


  if (extract(GDINA.obj, "att.str"))
    stop("Q-matrix validation is not available if attributes are structured.",
         call. = FALSE)

  if (eps > 1 || eps <= 0)
    stop("eps must be greater than 0 and less than 1.", call. = FALSE)

  Y <- extract(GDINA.obj, "seq.dat")
  Q <- extract(GDINA.obj, "Q")

  if (max(Q) > 1)
    stop("Q-validation can only be used for dichotomous attribute G-DINA model.",
         call. = FALSE)


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
  est.p <- rc / rn
  patt <- attributepattern(K)[-1, ]
  loc <- eta(patt) #2^K-1 x 2^K
  vsg <- varsigma(as.matrix(t(loc)), as.matrix(est.p), c(w))
  PVAF <- vsg / vsg[, L - 1]
  if (method == "PVAF") {
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
  } else{
    stop("Only PVAF is available.",call. = FALSE)
    out <- PVAF[, which(Kj == 1)]
    maxPVAF <- maxPVAF.loc <- NULL
    for (k in 1:max(Kj)) {
      maxPVAF <-
        cbind(maxPVAF, apply(PVAF[, which(Kj == k), drop = FALSE], 1, max))
      maxPVAF.loc <-
        cbind(maxPVAF.loc, apply(PVAF[, which(Kj == k), drop = FALSE], 1, which.max) +
                sum(Kj < k))

    }
    maxPVAF.change <- t(apply(cbind(0, maxPVAF), 1, diff))
    elbow <- out / (1 - PVAF)

    val_q <- patt[apply(elbow[, -ncol(elbow)], 1, which.max), ]
  }

  out.vsg <- round(t(vsg), digits)
  out.PVAF <- round(t(PVAF), digits)
  rownames(out.vsg) <-
    rownames(out.PVAF) <-
    apply(patt, 1, paste, collapse = "")
  Q <- data.frame(Q, row.names = extract(GDINA.obj, "item.names"))


  qvalid <-
    list(
      sug.Q = val_q,
      Q = Q,
      varsigma = out.vsg,
      PVAF = out.PVAF,
      eps = eps,
      est.p = est.p
    )
  class(qvalid) <- "Qval"
  return(qvalid)
}
