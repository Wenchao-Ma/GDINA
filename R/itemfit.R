#' Item fit statistics
#'
#' Calculate item fit statistics (Chen, de la Torre, & Zhang, 2013) and draw heatmap plot for item pairs
#'
#' @param GDINA.obj An estimated model object of class \code{GDINA}
#' @param person.sim Simulate expected responses from the posterior or based on EAP, MAP and MLE estimates.
#' @param p.adjust.methods p-values for the proportion correct, transformed correlation, and log-odds ratio
#'  can be adjusted for multiple comparisons at test and item level. This is conducted using \code{p.adjust} function in \pkg{stats},
#'  and therefore all adjustment methods supported by \code{p.adjust} can be used, including \code{"holm"},
#'  \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"} and \code{"BY"}. See \code{p.adjust}
#'  for more details. \code{"holm"} is the default.
#' @param N.resampling the sample size of resampling. By default, it is the maximum of 1e+5 and ten times of current sample size.
#' @param randomseed random seed; This is used to make sure the results are replicable. The default random seed is 123456.
#' @param cor.use how to deal with missing values when calculating correlations? This argument will be passed to \code{use} when calling \code{stats::cor}.
#' @param digits How many decimal places in each number? The default is 4.
#' @return an object of class \code{itemfit} consisting of several elements that can be extracted using
#'  method \code{extract}. Components that can be extracted include:
#' \describe{
#' \item{p}{the proportion correct statistics, adjusted and unadjusted p values for each item}
#' \item{r}{the transformed correlations, adjusted and unadjusted p values for each item pair}
#' \item{logOR}{the log odds ratios, adjusted and unadjusted p values for each item pair}
#' \item{maxitemfit}{the maximum proportion correct, transformed correlation, and log-odds ratio for each item with associated item-level adjusted p-values}
#' }
#'
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @export
#' @references
#' Chen, J., de la Torre, J., & Zhang, Z. (2013). Relative and Absolute Fit Evaluation in Cognitive Diagnosis Modeling.
#' \emph{Journal of Educational Measurement, 50}, 123-140.
#'
#' @examples
#' \dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#'
#' mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' mod1
#' itmfit <- itemfit(mod1)
#'
#' # Print "test-level" item fit statistics
#' # p-values are adjusted for multiple comparisons
#' # for proportion correct, there are J comparisons
#' # for log odds ratio and transformed correlation,
#' # there are J*(J-1)/2 comparisons
#'
#' itmfit
#'
#' # The following gives maximum item fit statistics for
#' # each item with item level p-value adjustment
#' # For each item, there are J-1 comparisons for each of
#' # log odds ratio and transformed correlation
#' summary(itmfit)
#'
#' # use extract to extract various components
#' extract(itmfit,"r")
#'
#' mod2 <- GDINA(dat,Q,model="DINA")
#' itmfit2 <- itemfit(mod2)
#' #misfit heatmap
#' plot(itmfit2)
#' itmfit2
#'}


itemfit <- function(GDINA.obj,person.sim = "post",p.adjust.methods = "holm",
                    cor.use = "pairwise.complete.obs",
                    digits = 4,N.resampling = NULL,randomseed=123456){


  if(!class(GDINA.obj)=="GDINA") stop("GDINA.obj must be a GDINA estimate.",call. = FALSE)
  if(extract(GDINA.obj,"ngroup")>1) stop("Itemfit is not available for multiple group estimation.",call. = FALSE)
  if (extract(GDINA.obj, "sequential")) stop("Itemfit is not available for sequential models.", call. = FALSE)
  itemfitcall <- match.call()
  if (exists(".Random.seed", .GlobalEnv))
    oldseed <- .GlobalEnv$.Random.seed
  else
    oldseed <- NULL
  set.seed(randomseed)
  dat <- as.matrix(extract(GDINA.obj, "dat"))
  Q <- extract(GDINA.obj, "Q")
  Qc <- extract(GDINA.obj, "Qc")
  K <- extract(GDINA.obj, "natt")
  N <- extract(GDINA.obj, "nobs")
  J <- extract(GDINA.obj, "nitem")
  Pr <- t(extract(GDINA.obj, "LCprob.parm"))  #L x J


# -------------Item Fit----------------------#

  # Generating model-based responses
  if (is.null(N.resampling)) {
    Nfit <- max(1e+05, 10 * N)
  } else{
    Nfit <- N.resampling
  }
  Rep <- ceiling(Nfit / N)
  pattern <- t(extract(GDINA.obj, "attributepattern"))

  if (person.sim == "post") {
    post <- extract(GDINA.obj, "posterior.prob")
    att_group <- sample(1:length(post), Rep * N, replace = TRUE, prob = post)
  }else{
    if (max(Q) > 1)
      person.sim <- "EAP"

    att <- switch(
      tolower(person.sim),
      eap = personparm.GDINA(GDINA.obj, what = "EAP"),
      mle = personparm.GDINA(GDINA.obj, what = "MLE")[, 1:K],
      map = personparm.GDINA(GDINA.obj, what = "MAP")[, 1:K]
    )
    att_group <-
      apply(att, 1, function(x)
        which.max(colSums(x == pattern)))[rep(1:N, Rep)]
  }



  Yfit <- Pr[att_group, ] > matrix(runif(length(att_group) * J), ncol = J)

  if(any(is.na(dat))){
    fitstat <- fitstats(dat,Yfit,FALSE)
    fitstat$r <- cor(dat, use = cor.use)
    fitstat$rfit <- cor(Yfit, use = cor.use)
  }else{
     fitstat <- fitstats(dat,Yfit,TRUE)
  }

  itempair <- NULL
  for (i in 1:(J - 1)) {
    itempair <- rbind(itempair, cbind(i, seq(i + 1, J)))
  }
  r.pairs <-
    data.frame(expected.r = fitstat$rfit[lower.tri(fitstat$rfit, diag = FALSE)],
               observed.r = fitstat$r[lower.tri(fitstat$r, diag = FALSE)])
  r.pairs$expected.fisherZ <-
    0.5 * log((1 + r.pairs$expected.r) / (1 - r.pairs$expected.r))
  r.pairs$observed.fisherZ <-
    0.5 * log((1 + r.pairs$observed.r) / (1 - r.pairs$observed.r))
  r.pairs$rstat <-
    abs(r.pairs$observed.fisherZ - r.pairs$expected.fisherZ)
  r.pairs$rstat.SE <- sqrt(1 / (N - 3))
  r.pairs$zstat <- r.pairs$rstat / r.pairs$rstat.SE
  r.pairs$unadj.pvalue <- pnorm(r.pairs$zstat, lower.tail = FALSE) * 2
  r.pairs$test.adj.pvalue <-
    stats::p.adjust(r.pairs$unadj.pvalue, method = p.adjust.methods)
  # r.pairs$item.adj.pvalue <- stats::p.adjust(r.pairs$unadj.pvalue,method = p.adjust.methods)
  l.pairs <-
    data.frame(expected.logOR = log(fitstat$lfit[lower.tri(fitstat$lfit, diag = FALSE)]),
               observed.logOR = log(fitstat$l[lower.tri(fitstat$l, diag = FALSE)]))
  l.pairs$lstat <-
    abs(l.pairs$observed.logOR - l.pairs$expected.logOR)
  l.pairs$lstat.SE <-
    sqrt(fitstat$sefit[lower.tri(fitstat$sefit, diag = FALSE)])
  l.pairs$zstat <- l.pairs$lstat / l.pairs$lstat.SE

  l.pairs$unadj.pvalue <- pnorm(l.pairs$zstat, lower.tail = FALSE) * 2
  l.pairs$test.adj.pvalue <-
    stats::p.adjust(l.pairs$unadj.pvalue, method = p.adjust.methods)
  p <- data.frame(expected.p=c(fitstat$pfit),
                observed.p=colMeans(dat,na.rm = TRUE))
  p$pstat <- abs(p$expected.p - p$observed.p)
  p$pstat.SE <- sqrt(c(fitstat$pfit) * (1 - c(fitstat$pfit)) / N)
  p$zstat <- p$pstat / p$pstat.SE
  p$unadj.pvalue <- pnorm(p$zstat, lower.tail = FALSE) * 2
  p$test.adj.pvalue <-
    stats::p.adjust(p$unadj.pvalue, method = p.adjust.methods)

  max.itemlevel.fit <- matrix(NA, J, 6)
  for (j in 1:J) {
    loc <- which(apply(itempair == j, 1, any))
    max.itemlevel.fit[j, ] <-
      c(
        max(r.pairs$zstat[loc]),
        r.pairs$unadj.pvalue[loc][which.max(r.pairs$zstat[loc])],
        stats::p.adjust(r.pairs$unadj.pvalue[loc], method = p.adjust.methods)[which.max(r.pairs$zstat[loc])],
        max(l.pairs$zstat[loc]),
        l.pairs$unadj.pvalue[loc][which.max(l.pairs$zstat[loc])],
        stats::p.adjust(l.pairs$unadj.pvalue[loc], method = p.adjust.methods)[which.max(l.pairs$zstat[loc])]
      )
  }
  max.itemlevel.fit <- round(cbind(p$zstat, p$unadj.pvalue, max.itemlevel.fit), digits)
  colnames(max.itemlevel.fit) <-
    c(
      "z.prop",
      "pvalue[z.prop]",
      "max[z.r]",
      "pvalue.max[z.r]",
      "adj.pvalue.max[z.r]",
      "max[z.logOR]",
      "pvalue.max[z.logOR]",
      "adj.pvalue.max[z.logOR]"
    )
  rownames(max.itemlevel.fit) <- extract(GDINA.obj,"item.names")
  r.pairs <-
    data.frame(item.pair.1 = itempair[, 1],
               item.pair.2 = itempair[, 2],
               round(r.pairs, digits))
  l.pairs <-
    data.frame(item.pair.1 = itempair[, 1],
               item.pair.2 = itempair[, 2],
               round(l.pairs, digits))
  p <- data.frame(item = 1:J, round(p, digits),row.names = extract(GDINA.obj,"item.names"))

  if (!is.null(oldseed))
    .GlobalEnv$.Random.seed <- oldseed
  else
    rm(".Random.seed", envir = .GlobalEnv)

  output <-
    list(
      r = r.pairs,
      p = p,
      logOR = l.pairs,
      max.itemlevel.fit = max.itemlevel.fit,
      options = list(
        person.sim = person.sim,
        p.adjust.methods = p.adjust.methods,
        digits = digits,
        N.resampling = N.resampling,
        randomseed = randomseed,
        call = itemfitcall
      )
    )
  class(output) <- "itemfit"

  return(output)
}



