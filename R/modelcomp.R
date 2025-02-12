#' Item-level model comparison using Wald, LR or LM tests
#'
#' This function evaluates whether the saturated G-DINA model can be replaced by reduced
#' CDMs without significant loss in model data fit for each item using the Wald test, likelihood ratio (LR) test or Lagrange multiplier (LM) test.
#' For Wald test, see de la Torre (2011), de la Torre and Lee (2013), Ma, Iaconangelo and de la Torre (2016) and Ma & de la Torre (2018) for details.
#' For LR test and a two-step LR approximation procedure, see Sorrel, de la Torre, Abad, and Olea (2017), Ma (2017) and Ma & de la Torre (2019).
#' For LM test, which is only applicable for DINA, DINO and ACDM, see Sorrel, Abad, Olea, de la Torre, and Barrada (2017).
#' This function also calculates the dissimilarity
#' between the reduced models and the G-DINA model, which can be viewed as a measure of effect size (Ma, Iaconangelo & de la Torre, 2016).
#'
#' After the test statistics for each reduced CDM were calculated for each item, the
#' reduced models with p values less than the pre-specified alpha level were rejected.
#' If all reduced models were rejected for an item, the G-DINA model was used as the best model;
#' if at least one reduced model was retained, two diferent rules can be implemented for selecting
#' the best model specified in argument \code{decision.args}:
#'
#' (1) when \code{rule="simpler"},
#'
#'  If (a) the DINA or DINO model
#'  was one of the retained models, then the DINA or DINO model with the larger p
#'  value was selected as the best model; but if (b) both DINA and DINO were rejected, the reduced
#'  model with the largest p value was selected as the best model for this item. Note that
#'  when the p-values of several reduced models were greater than 0.05, the DINA and DINO models were
#'  preferred over the A-CDM, LLM, and R-RUM because of their simplicity.
#'
#'  (2) When \code{rule="largestp"} (default),
#'
#'  The reduced model with the largest p-values is selected as the most appropriate model.
#'
#' @param GDINA.obj An estimated model object of class \code{GDINA}
#' @param method method for item level model comparison; can be \code{wald}, \code{LR} or \code{LM}.
#' @param items a vector of items to specify the items for model comparsion
#' @param models a vector specifying which reduced CDMs are possible reduced CDMs for each
#'   item. The default is "DINA","DINO","ACDM","LLM",and "RRUM".
#' @param DS whether dissimilarity index should be calculated? \code{FALSE} is the default.
#' @param decision.args a list of options for determining the most appropriate models including (1) \code{rule} can be
#'   either \code{"simpler"} or \code{"largestp"}. See details;
#'   (2) \code{alpha.level} for the nominal level of decision; and (3) \code{adjusted} can be either \code{TRUE} or \code{FALSE}
#'   indicating whether the decision is based on p value (\code{adjusted = FALSE}) or adjusted p values.
#' @param Wald.args a list of options for Wald test including (1) \code{SE.type} giving the type of
#'   covariance matrix for the Wald test; by default, it uses outer product of gradient based on incomplete information matrix;
#'   (2) \code{varcov} for user specified variance-covariance matrix. If supplied, it must
#'   be a list, giving the variance covariance matrix of success probability for each item or category.
#'   The default is \code{NULL}, in which case, the estimated variance-covariance matrix from the GDINA function
#'   is used.
#' @param LR.args a list of options for LR test including for now only \code{LR.approx}, which is either \code{TRUE} or \code{FALSE},
#'   indicating whether a two-step LR approximation is implemented or not.
#' @param LM.args a list of options for LM test including \code{reducedMDINA}, \code{reducedMDINO}, and \code{reducedMACDM} for
#' DINA, DINO and ACDM estimates from the \code{GDINA} function; \code{SE.type} specifies the type of covariance matrix.
#' @param p.adjust.methods adjusted p-values for multiple hypothesis tests. This is conducted using \code{p.adjust} function in \pkg{stats},
#'  and therefore all adjustment methods supported by \code{p.adjust} can be used, including \code{"holm"},
#'  \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"} and \code{"BY"}. See \code{p.adjust}
#'  for more details. \code{"holm"} is the default, indicating the Holm method.
#'
#'
#' @return an object of class \code{modelcomp}. Elements that can be
#' extracted using \code{extract} method include
#' \describe{
#' \item{stats}{Wald or LR statistics}
#' \item{pvalues}{p-values associated with the test statistics}
#' \item{adj.pvalues}{adjusted p-values}
#' \item{df}{degrees of freedom}
#' \item{DS}{dissimilarity between G-DINA and other CDMs}
#' }
#' @author Wenchao Ma, The University of Minnesota, \email{wma@umn.edu}
#' Miguel A. Sorrel, Universidad Autonoma de Madrid
#' Jimmy de la Torre, The University of Hong Kong
#'
#' @seealso \code{\link{GDINA}}, \code{\link{autoGDINA}}
#'
#'@references
#' de la Torre, J., & Lee, Y. S. (2013). Evaluating the wald test for item-level comparison of
#' saturated and reduced models in cognitive diagnosis. \emph{Journal of Educational Measurement, 50}, 355-373.
#'
#' Ma, W., Iaconangelo, C., & de la Torre, J. (2016). Model similarity, model selection and attribute classification.
#' \emph{Applied Psychological Measurement, 40}, 200-217.
#'
#' Ma, W. (2017). \emph{A Sequential Cognitive Diagnosis Model for Graded Response: Model Development, Q-Matrix Validation,and Model Comparison. Unpublished doctoral dissertation.} New Brunswick, NJ: Rutgers University.
#'
#' Ma, W., & de la Torre, J. (2019). Category-Level Model Selection for the Sequential G-DINA Model. \emph{Journal of Educational and Behavioral Statistics}. 44, 61-82.
#'
#' Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive Diagnosis Modeling. \emph{Journal of Statistical Software, 93(14)}, 1-26.
#'
#' Sorrel, M. A., Abad, F. J., Olea, J., de la Torre, J., & Barrada, J. R. (2017). Inferential Item-Fit Evaluation in Cognitive Diagnosis Modeling. \emph{Applied Psychological Measurement, 41,} 614-631.
#'
#' Sorrel, M. A., de la Torre, J., Abad, F. J., & Olea, J. (2017). Two-Step Likelihood Ratio Test for Item-Level Model Comparison in Cognitive Diagnosis Models. \emph{Methodology, 13}, 39-47.
#'
#'
#'
#'@examples
#'\dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' # --- GDINA model ---#
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' fit
#'
#' ###################
#' #
#' # Wald test
#' #
#' ###################
#'
#' w <- modelcomp(fit)
#' w
#' # wald statistics
#' extract(w,"stats")
#' #p values
#' extract(w,"pvalues")
#' # selected models
#' extract(w,"selected.model")
#' ##########################
#' #
#' # LR and Two-step LR test
#' #
#' ##########################
#'
#' lr <- modelcomp(fit,method = "LR")
#' lr
#' TwostepLR <- modelcomp(fit,items =c(6:10),method = "LR",LR.args = list(LR.approx = TRUE))
#' TwostepLR
#'
#' ##########################
#' #
#' # LM test
#' #
#' ##########################
#'
#' dina <- GDINA(dat = dat, Q = Q, model = "DINA")
#' dino <- GDINA(dat = dat, Q = Q, model = "DINO")
#' acdm <- GDINA(dat = dat, Q = Q, model = "ACDM")
#' lm <- modelcomp(method = "LM",LM.args=list(reducedMDINA = dina,
#' reducedMDINO = dino, reducedMACDM  = acdm))
#' lm
#'
#'
#' }
## @import MASS
#' @export
modelcomp <- function(GDINA.obj=NULL,method = "Wald",items = "all", p.adjust.methods = "holm",
                      models=c("DINA","DINO","ACDM","LLM","RRUM"),
                      decision.args = list(rule = "largestp", alpha.level = 0.05, adjusted = FALSE), DS = FALSE,
                      Wald.args = list(SE.type = 2,varcov = NULL),
                      LR.args = list(LR.approx = FALSE),
                      LM.args = list(reducedMDINA = NULL, reducedMDINO = NULL, reducedMACDM  = NULL, SE.type = 2)){
  if(is.null(GDINA.obj)&method!="LM") stop("GDINA.obj must be provided unless LM test is requested.",call. = FALSE)

  if(method=="LM"){
    if(is.null(LM.args$reducedMDINA)&&is.null(LM.args$reducedMDINO)&&is.null(LM.args$reducedMACDM)){
    stop("LM.args must be specified when LM method is selected.",call. = FALSE)
    }else{
      reducedMDINA <- LM.args$reducedMDINA
      reducedMDINO <- LM.args$reducedMDINO
      reducedMACDM <- LM.args$reducedMACDM
      if(!is.null(reducedMDINA)){
        dat <- reducedMDINA$options$dat
        Q <- reducedMDINA$options$Q
        J <- nrow(Q)
        item.names <- extract(reducedMDINA,"item.names")
      }else if(!is.null(reducedMDINO)){
        dat <- reducedMDINO$options$dat
        Q <- reducedMDINO$options$Q
        J <- nrow(Q)
        item.names <- extract(reducedMACDM,"item.names")
      }else{
        dat <- reducedMACDM$options$dat
        Q <- reducedMACDM$options$Q
        J <- nrow(Q)
        item.names <- extract(reducedMACDM,"item.names")
      }
    if(is.null(LM.args$SE.type)) LM.args$SE.type <- 2
    }


  }else{
    stopifnot(isa(GDINA.obj,"GDINA"))
    if (!all(toupper(extract(GDINA.obj,"models"))%in%c("GDINA","LOGGDINA","LOGITGDINA")))
      stop ("Implementing the Wald and LR tests for item-level model comparison requires all items to be fitted by the G-DINA model.",call. = FALSE)
    if(!is.null(extract(GDINA.obj,"att.str")))
      stop("Item-level model comparison is not available for structured attributes.",call. = FALSE)
    Q <- 1*(extract(GDINA.obj,"Q")>0.5)
    item.names <- extract(GDINA.obj,"item.names")
  }

  my.decision.args = list(rule = "simpler", alpha.level = 0.05, adjusted = FALSE)
  my.Wald.args = list(SE.type = 2,varcov = NULL)
  my.LR.args = list(LR.approx = FALSE)
  my.LM.args = list(reducedMDINA = NULL, reducedMDINO = NULL, reducedMACDM  = NULL, SE.type = 2)

  decision.args <- utils::modifyList(my.decision.args,decision.args)
  Wald.args <- utils::modifyList(my.Wald.args,Wald.args)
  LR.args <- utils::modifyList(my.LR.args,LR.args)
  LM.args <- utils::modifyList(my.LM.args,LM.args)

  allModels <- c("DINA","DINO","ACDM","LLM","RRUM")
  model_num <- pmatch(toupper(models),allModels)
  J <- length(item.names)
  Kjs <- rowSums(Q)
  Ks <- cumsum(2^Kjs)

  if(length(items)==1){
    if(tolower(items)=="all") {
      items <- which(rowSums(Q)>1)
    }else {
      if (!is.numeric(items)) stop("item must be specified as item numbers.",call. = FALSE)
      if (sum(Q[items,])==1) {
        items <- NULL
        stop("Wald test can only be used for the items requiring more than 1 attributes.",call. = FALSE)
      }
    }
  }else {
    if (!all(sapply(items,is.numeric))) stop("Items must be specified as item numbers.",call. = FALSE)

    items <- intersect(which(Kjs>1),items)
  }
  stopifnot(!is.null(items)&&!is.null(models))
  neg2LL <- MCstat <- pvalues <- df <- matrix(NA,length(items),length(allModels))
  if (tolower(method)=="wald"){

    for (i in seq(length(items))){ # for each item
      j <- items[i]
      Kj <- rowSums(Q[items,,drop=FALSE]) # Kj for each item
      RDA <- RDINA(unique(Kj))
      RDO <- RDINO(unique(Kj))
      RAM <- RACDM(unique(Kj))

      # variance covariance matrix for item j
      if (is.null(Wald.args$varcov)) {
        cov <- extract(GDINA.obj,"catprob.cov",SE.type = Wald.args$SE.type)
        ind <- cov$index
        vcov <- cov$cov[ind[which(ind$cat==j),5],ind[which(ind$cat==j),5]]
        }
      else{vcov <- Wald.args$varcov[[j]]}

      if (any(model_num==1)){
        MCstat[i,1] <- t(RDA[[Kjs[j]]]%*%extract(GDINA.obj,"catprob.parm")[[j]])%*%
          MASS::ginv(RDA[[Kjs[j]]]%*%vcov%*%t(RDA[[Kjs[j]]]))%*%
          (RDA[[Kjs[j]]]%*%extract(GDINA.obj,"catprob.parm")[[j]])
        df[i,1] <- 2^Kjs[j]-2
        pvalues[i,1] <- pchisq(MCstat[i,1],df[i,1],lower.tail = F)
      }
      if (any(model_num==2)){
        MCstat[i,2] <- t(RDO[[Kjs[j]]]%*%extract(GDINA.obj,"catprob.parm")[[j]])%*%
          MASS::ginv(RDO[[Kjs[j]]]%*%vcov%*%t(RDO[[Kjs[j]]]))%*%
          (RDO[[Kjs[j]]]%*%extract(GDINA.obj,"catprob.parm")[[j]])
        df[i,2] <- 2^Kjs[j]-2
        pvalues[i,2] <- pchisq(MCstat[i,2],df[i,2],lower.tail = F)
      }
      if (any(model_num==3)){
        MCstat[i,3] <- t(RAM[[Kjs[j]]]%*%extract(GDINA.obj,"catprob.parm")[[j]])%*%
          MASS::ginv(RAM[[Kjs[j]]]%*%vcov%*%t(RAM[[Kjs[j]]]))%*%
          (RAM[[Kjs[j]]]%*%extract(GDINA.obj,"catprob.parm")[[j]])
        df[i,3] <- 2^Kjs[j]-Kjs[j]-1
        pvalues[i,3] <- pchisq(MCstat[i,3],df[i,3],lower.tail = F)

      }
      if (any(model_num==4)){
        fpj <- logit(extract(GDINA.obj,"catprob.parm")[[j]])
        grad_fpj <- solve(diag(extract(GDINA.obj,"catprob.parm")[[j]]*(1-extract(GDINA.obj,"catprob.parm")[[j]])))
        var_fpj <- grad_fpj%*%vcov%*%t(grad_fpj)
        MCstat[i,4] <- t(RAM[[Kjs[j]]]%*%fpj)%*% MASS::ginv(RAM[[Kjs[j]]]%*%var_fpj%*%t(RAM[[Kjs[j]]]))%*%(RAM[[Kjs[j]]]%*%fpj)
        df[i,4] <- 2^Kjs[j]-Kjs[j]-1
        pvalues[i,4] <- pchisq(MCstat[i,4],df[i,4],lower.tail = F)
      }
      if (any(model_num==5)){
        fpj <- log(extract(GDINA.obj,"catprob.parm")[[j]])
        grad_fpj <- solve(diag(extract(GDINA.obj,"catprob.parm")[[j]]))
        var_fpj <- grad_fpj%*%vcov%*%t(grad_fpj)
        MCstat[i,5] <- t(RAM[[Kjs[j]]]%*%fpj)%*%MASS::ginv(RAM[[Kjs[j]]]%*%var_fpj%*%t(RAM[[Kjs[j]]]))%*% (RAM[[Kjs[j]]]%*%fpj)
        df[i,5] <- 2^Kjs[j]-Kjs[j]-1
        pvalues[i,5] <- pchisq(MCstat[i,5],df[i,5],lower.tail = F)
      }

    }
  }else if(toupper(method)=="LR"){


    fun.args <- as.list(GDINA.obj$extra$call)[-c(1:3)]

    for (i in seq(length(items))){ # for each item
      j <- items[i]
      m <- rep(0,nrow(Q))
      if(LR.args$LR.approx){
        maxitr <- rep(0,nrow(Q))
        maxitr[j] <- 1
        att.dist <- "fixed"
      }else{
        maxitr <- rep(1000,nrow(Q))
        att.dist <- extract(GDINA.obj,"att.dist")
      }
      for(mm in model_num){
        m[j] <- mm
        updated.args <- utils::modifyList(fun.args,list(dat = quote(extract(GDINA.obj,"dat")),
                                                 Q = quote(extract(GDINA.obj,"originalQ")), att.dist = att.dist,
                                                 model=m,verbose=0, control = list(maxitr = maxitr),
                                                 catprob.parm=quote(extract(GDINA.obj,"catprob.parm")),
                                                 att.prior = quote(extract(GDINA.obj,"att.prior"))))
        tmp.est <- do.call(GDINA,updated.args)
        neg2LL[i,mm] <- deviance(tmp.est)
        MCstat[i,mm] <- deviance(tmp.est) - deviance(GDINA.obj)
        if(MCstat[i,mm]<=0){
          # MCstat[i,mm] <- .Machine$double.eps
          warning("LR statistic is negative.",call. = FALSE)
        }
        df[i,mm] <- npar(GDINA.obj)$`No. of total item parameters` - npar(tmp.est)$`No. of total item parameters`
        pvalues[i,mm] <- pchisq(MCstat[i,mm],df[i,mm],lower.tail = F)
      }
    }
    colnames(neg2LL) <- allModels
    rownames(neg2LL) <- extract(GDINA.obj,"item.names")[items]
  }else if(toupper(method)=="LM"){


    if ("DINA" %in% models) {
      if (is.null(reducedMDINA) | (all(reducedMDINA$model == "DINA") == FALSE)) {
        stop("Implementing the LM test for the DINA model requires all items are fitted by the DINA model.",
             call. = FALSE)
      }

      vDINA <- itemprob_se_M(reducedMDINA, type = LM.args$SE.type)
      varmat.probDINA <- list()
      for (jj in 1:nrow(Q)) {varmat.probDINA[[jj]]<- vDINA$cov[vDINA$ind$item ==jj , vDINA$ind$item ==jj]}

    }
    if ("DINO" %in% models) {
      if (is.null(reducedMDINO) | (all(reducedMDINO$model == "DINO") == FALSE)) {
        stop("Implementing the LM test for the DINO model requires all items are fitted by the DINO model.",
             call. = FALSE)
      }
      vDINO <- itemprob_se_M(reducedMDINO, type = LM.args$SE.type)
      varmat.probDINO <- list()
      for (jj in 1:nrow(Q)) {varmat.probDINO[[jj]]<- vDINO$cov[vDINO$ind$item ==jj , vDINO$ind$item ==jj]}

    }
    if ("ACDM" %in% models) {
      if (is.null(reducedMACDM) | (all(reducedMACDM$model == "ACDM") == FALSE)) {
        stop("Implementing the LM test for the ACDM model requires all items are fitted by the ACDM model.",
             call. = FALSE)
      }
      vACDM <- itemprob_se_M(reducedMACDM, type = LM.args$SE.type)
      varmat.probACDM <- list()
      for (jj in 1:nrow(Q)) {varmat.probACDM[[jj]]<- vACDM$cov[vACDM$ind$item ==jj , vACDM$ind$item ==jj]}

    }

    for (i in seq(length(items))){ # for each item
      j <- items[i]
Kj <- Kjs[j]
        if ("DINA" %in% models) {
          prob.j <- extract(reducedMDINA,"catprob.parm")[[j]]
          nj <- extract(reducedMDINA,"expectedTotal")[j,seq(prob.j)]
          rj <- extract(reducedMDINA,"expectedCorrect")[j,seq(prob.j)]
          gradient <- matrix((1/(prob.j*(1-prob.j)))*(rj-prob.j*nj),nrow = 1)
          MCstat[i,1] <- gradient%*%varmat.probDINA[[j]]%*%t(gradient)
          df[i,1] <- abs(2^Kj-2)
          pvalues[i,1] <- 1 - pchisq(MCstat[i,1], df = abs(2^Kj-2))
          # gradientDINA[[j]] <- gradient
        }

        if ("DINO" %in% models) {
          prob.j <- extract(reducedMDINO,"catprob.parm")[[j]]
          nj <- extract(reducedMDINO,"expectedTotal")[j,seq(prob.j)]
          rj <- extract(reducedMDINO,"expectedCorrect")[j,seq(prob.j)]
          gradient <- matrix((1/(prob.j*(1-prob.j)))*(rj-prob.j*nj),nrow = 1)
          MCstat[i,2] <- gradient%*%varmat.probDINO[[j]]%*%t(gradient)
          df[i,2] <- abs(2^Kj-2)
          pvalues[i,2] <- 1 - pchisq(MCstat[i,2], df = abs(2^Kj-2))
        }

        if ("ACDM" %in% models) {
          prob.j <- extract(reducedMACDM,"catprob.parm")[[j]]
          nj <- extract(reducedMACDM,"expectedTotal")[j,seq(prob.j)]
          rj <- extract(reducedMACDM,"expectedCorrect")[j,seq(prob.j)]
          gradient <- matrix(prob.j*(1-prob.j)^(-1)*(rj-prob.j*nj),nrow = 1)
          MCstat[i,3] <- gradient%*%varmat.probACDM[[j]]%*%t(gradient)
          df[i,3] <- length(prob.j)-(Kj+1)
          pvalues[i,3] <- 1 - pchisq(MCstat[i,3], df = df[i,3])

        }

      }
    }


  ##### dissimilarity index
  ds.f <- NULL
  if(!is.null(GDINA.obj)){

    if(DS){
      ds <- lapply(extract(GDINA.obj,"catprob.parm")[items],function(p){
        unlist(lapply(models,function(m) DS(p,m)$DS))
      })
      ds <- do.call(rbind,ds)
      ds.f <- matrix(NA,nrow(ds),5)
      if(length(models)<5) {
        ds.f[,model_num] <- ds
      }else{
        ds.f <- ds
      }
      colnames(ds.f) <- allModels
      # ds.f <- round(ds.f,digits)
    }

  }

  adj.pvalues <- pvalues
  adj.pvalues[!is.na(pvalues)] <- stats::p.adjust(pvalues[!is.na(pvalues)], method = p.adjust.methods)


  ## selected CDMs
  # print(pvalues)

  if (decision.args$adjusted == TRUE){
    ps <- adj.pvalues
  }else{
    ps <- pvalues
  }

  if (tolower(decision.args$rule) == "simpler") {

    if(any(toupper(models)%in%c("DINA","DINO"))){
      if(any(toupper(models)%in%c("LLM","RRUM","ACDM"))){
        #some DINO and DINA; some acdm
        m = apply(ps, 1, function(x) {
          if (max(x[1:2], na.rm = TRUE) > decision.args$alpha.level) {
            which.max.randomtie(x[1:2])
          } else if (max(x[3:5], na.rm = TRUE) > decision.args$alpha.level) {
            which.max.randomtie(x[3:5]) + 2
          } else{
            return(0)
          }
        })
      }else{
        #no acdm
        m <- apply(ps,1,which.max.randomtie)
        m[apply(ps,1,max,na.rm=TRUE)<decision.args$alpha.level] <- 0
      }
    }else{
      # no dina and dino
      m <- apply(ps,1,which.max.randomtie)
      m[apply(ps,1,max,na.rm=TRUE)<decision.args$alpha.level] <- 0
    }
  } else if(tolower(decision.args$rule) == "largestp"){
    m = apply(ps, 1, function(x) {
      if (max(x, na.rm = T) > decision.args$alpha.level) {
        which.max.randomtie(x)
      } else{
        return(0)
      }
    })
  }
# print(m)

  if (decision.args$adjusted == TRUE){
    p.na <- cbind(NA,pvalues)
    adj.na <- cbind(NA,ps)
  }else{
    p.na <- cbind(NA,ps)
    adj.na <- cbind(NA,adj.pvalues)
    }
m.adj <- adj.na[matrix(c(seq_len(nrow(ps)),m+1),ncol = 2)]
m.p <- p.na[matrix(c(seq_len(nrow(pvalues)),m+1),ncol = 2)]
  ret  <-  data.frame(models = rep("GDINA",J), pvalues = rep(NA,J), adj.pvalues = rep(NA,J),stringsAsFactors = FALSE)
  ret$models[items] <- c("GDINA","DINA","DINO","ACDM","LLM","RRUM")[m+1]
  ret$pvalues[items] <- round(m.p,4)
  ret$adj.pvalues[items] <- round(m.adj,4)
  row.names(ret) <- item.names


  colnames(MCstat) <- colnames(df) <- colnames(pvalues) <- colnames(adj.pvalues) <- allModels
  rownames(MCstat) <- rownames(df) <- rownames(pvalues) <- rownames(adj.pvalues) <- item.names[items]
  out <- list(stats=MCstat,pvalues=pvalues,adj.pvalues = adj.pvalues,df=df,DS=ds.f,models = models,method=method,
              DS = DS, Wald.args = Wald.args, LR.args = LR.args,neg2LL = neg2LL, decision.args = decision.args,
              selected.model=ret,p.adjust.methods = p.adjust.methods)
  class(out) <- "modelcomp"
  return(out)

}
