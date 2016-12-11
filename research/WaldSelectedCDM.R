#' Q-matrix validation, model selection and calibration in one run
#'
#' \code{autoGDINA} conducts a series of CDM analyses within the G-DINA framework. Particularly,
#' the GDINA model is fitted to the data first using the \code{GDINA} function;
#' then, the Q-matrix is validated using the function \code{Qval}.
#' Based on the suggested Q-matrix, the data is fitted by the G-DINA model again, followed
#' by an item level model selection via the Wald test (See \code{modelcomp}. Lastly,
#' the selected models are calibrated based on the suggested Q-matrix using the \code{GDINA} function.
#' The Q-matrix validation and item-level model selection can be disabled by the users.
#' Possible reduced CDMs for Wald test include the DINA model, the DINO model, A-CDM, LLM and RRUM.
#' Users can specify all or some of them for model selection. See \code{Details}
#' for the rules of item-level model selection.
#'
#'
#' @details
#'
#' This section describes how reduced models are selected at item level:
#'
#' After the Wald statistics for each reduced CDM were calculated for each item, the
#' reduced models with p values less than the pre-specified alpha level were rejected.
#' If all reduced models were rejected for an item, the G-DINA model was used as the best model;
#' if at least one reduced model was retained, two diferent rules can be implemented for selecting
#' the best model:
#'
#' when \code{modelselectionrule} is \code{simpler}:
#'
#'  If (a) the DINA or DINO model
#'  was one of the retained models, then the DINA or DINO model with the larger p
#'  value was selected as the best model; but if (b) both DINA and DINO were rejected, the reduced
#'  model with the largest p value was selected as the best model for this item. Note that
#'  when the p-values of several reduced models were greater than 0.05, the DINA and DINO models were
#'  preferred over the A-CDM, LLM, and R-RUM because of their simplicity. This procedure is originally
#'  proposed by Ma, Iaconangelo, and de la Torre (2016).
#'
#'  When \code{modelselectionrule} is \code{largestp}:
#'
#'  The reduced model with the largest p-values is selected as the most appropriate model.
#'
#'   When \code{modelselectionrule} is \code{DS}:
#'
#'  The reduced model with non-significant p-values but smallest dissimilarity index is selected as the most appropriate model.
#'  Dissimilarity index can be viewed as an effect size, which quatifies how dis-similar the reduced model is from the
#'  G-DINA model. See Ma, Iaconangelo, and de la Torre (2016).
#'
#' @param dat A required \eqn{N \times J} \code{matrix} or \code{data.frame} consisting of binary
#'    responses of \eqn{N} examinees to \eqn{J} items. Missing values need to be coded as \code{NA}.
#' @param Q A required \eqn{J \times K} item by attribute association matrix (Q-matrix; Tatsuoka, 1983),
#'    where J represents test
#'    length and K represents the number of attributes. For binary attributes,
#'    1 denotes attributes are measured by items and 0 means attributes are not
#'    necessary.
#' @param Qvalid logical; validate Q-matrix or not? \code{TRUE} is the default.
#' @param alpha.level nominal level for the Wald test. The default is 0.05.
#' @param modelselection logical; conducting model selection or not?
#' @param modelselectionrule how to conducted model selection? Possible options include
#'    \code{simpler} and \code{largestp}. See \code{Details}.
#' @param eps cut-off value for PVAF if \code{Qvalid=TRUE}. The default is 0.95.
#' @param reducedCDM a vector specifying which reduced CDMs are possible reduced CDMs for each
#'   item. The default is "DINA","DINO","ACDM","LLM",and "RRUM".
#' @param GDINA.option a list of options for \code{\link{GDINA}} model estimation
#' @param SelectedCDM.option a list of options to estimate selected CDMs
#' @seealso \code{\link{GDINA}}, \code{\link{modelcomp}}, \code{\link{Qval}}
#'
#' @return CDM.obj an object of class \code{GDINA}. It is based on selected CDMs using the Wald test and suggested Q-matrix if \code{Qvalid=TRUE}.
#' @return Qval.obj The results of Q-validation if \code{Qvalid=TRUE}. It is an object of class \code{Qvalidation} with components:
#' \describe{
#' \item{zeta}{zeta index of each item for each possible q-vector}
#' \item{PVAF}{PVAF of each item for each possible q-vector}
#' \item{sugg.Q}{suggested Q-matrix}
#' }
#' @return Wald.obj a list consisting of wald statistics and associated p-values for model selection using the Wald test
#' @return GDINA.1 initial GDINA calibration
#' @return GDINA.2 GDINA calibration based on validated Q-matrix
#' @return options options of this function
#'
#'@examples
#' # simulated responses
#' Q <- sim10GDINA$simQ
#' dat <- sim10GDINA$simdat
#'
#' #misspecified Q
#' misQ <- Q
#' misQ[10,] <- c(0,1,0)
#' out1 <- autoGDINA(dat,misQ,modelselectionrule="largestp")
#' out1
#' summary(out1)
#' coef(out1)
#'
#' #using the other selection rule
#' out11 <- autoGDINA(dat,misQ,modelselectionrule="simpler")
#' out11
#' summary(out11)
#' coef(out11)
#'
#' # disable model selection function
#' out12 <- autoGDINA(dat,misQ,modelselection=FALSE)
#' out12
#' summary(out12)
#' coef(out12)
#'
#' # -- Only consider some reduced CDMs
#' out2 <- autoGDINA(dat,misQ,reducedCDM = c("DINA","LLM"))
#' #fitted models
#' out2$CDM.obj$model
#' summary(out2)
#' coef(out2)
#'
#' # Disable Q-matrix validation
#' out3 <- autoGDINA(dat = dat, Q = misQ, Qvalid = FALSE,alpha.level=0.01)
#' out3
#' summary(out3)
#' coef(out3)
#'
#' @export

autoGDINA <-
  function(dat, Q, reducedCDM = c("DINA", "DINO", "ACDM", "LLM", "RRUM"),
           alpha.level = 0.05, modelselection = TRUE,
           modelselectionrule = "simpler",
           Qvalid = TRUE, eps = 0.95, GDINA1.option = list(),
           GDINA2.option = list(),
           SelectedCDM.option = list()) {
    ASEcall <- match.call()
    # options(warn = -1)
    GDINA.obj <- GDINA.obj2 <- NULL
    CDM.opt <- GDINA2.opt <- GDINA1.opt <-list(model = "GDINA", sequential = FALSE, item.names = NULL,
      higher.order = FALSE, higher.order.model ="Rasch", higher.order.method = "MMLE", higher.order.SE = FALSE,
      verbose = TRUE, itemprob.parm = NULL, higher.order.parm = NULL,
      mono.constraint = FALSE, SE = TRUE, SE.type = 2,
      empirical = TRUE, att.prior = NULL, att.str = FALSE,
      lower.p = 0.0001,upper.p = 0.9999, smallNcorrection = c(0.0005,0.001),
      nstarts = 1, conv.crit = 0.001, conv.type = "max.p.change",maxitr = 1000,
      digits = 4,diagnosis = 0,Mstep.warning = FALSE,optimizer = "all",
      randomseed = 123456)
    if(any(names(GDINA1.option)%in%names(GDINA1.opt))) stop("GDINA1.option must be a list consisting of arguments for GDINA function.",call. = FALSE)
    if(any(names(GDINA2.option)%in%names(GDINA2.opt))) stop("GDINA2.option must be a list consisting of arguments for GDINA function.",call. = FALSE)
    if(any(names(SelectedCDM.option)%in%names(CDM.opt))) stop("SelectedCDM.option must be a list consisting of arguments for GDINA function.",call. = FALSE)
    GDINA1.opt[names(GDINA1.option)] <- GDINA1.option
    GDINA2.opt[names(GDINA2.option)] <- GDINA2.option
    CDM.opt[names(SelectedCDM.option)] <- SelectedCDM.option
    if (GDINA1.opt$higher.order == TRUE) GDINA1.opt$empirical <- FALSE
    if (GDINA2.opt$higher.order == TRUE) GDINA2.opt$empirical <- FALSE
    if (CDM.opt$higher.order == TRUE) CDM.opt$empirical <- FALSE
    cat("Initial calibration...")
    #print(GDINA.opt)
    GDINA1.obj <-
      GDINA(dat, Q, model = GDINA1.opt$model, sequential = GDINA1.opt$sequential,
        higher.order = GDINA1.opt$higher.order,
        higher.order.model = GDINA1.opt$higher.order.model,
        higher.order.method = GDINA1.opt$higher.order.method,
        higher.order.SE = GDINA1.opt$higher.order.SE,
        mono.constraint = GDINA1.opt$mono.constraint,
        SE = GDINA1.opt$SE,lower.p = GDINA1.opt$lower.p,upper.p = GDINA1.opt$upper.p,
        smallNcorrection = GDINA1.opt$smallNcorrection,
        SE.type = GDINA1.opt$SE.type, conv.type=GDINA1.opt$conv.type,
        empirical = GDINA1.opt$empirical,Mstep.warning=GDINA1.opt$Mstep.warning,
        verbose = GDINA1.opt$verbose,optimizer=GDINA1.opt$optimizer,
        nstarts = GDINA1.opt$nstarts,randomseed=GDINA1.opt$randomseed,
        att.prior = GDINA1.opt$att.prior,
        att.str = GDINA1.opt$att.str,
        higher.order.parm = GDINA1.opt$higher.order.parm,
        digits = GDINA1.opt$digits,
        conv.crit = GDINA1.opt$conv.crit,
        diagnosis = GDINA1.opt$diagnosis,
        maxitr = GDINA1.opt$maxitr)
    Q <- show.GDINA(GDINA1.obj,what = "Q")
    cat("done.\n")
    Qv <- NULL
    if (Qvalid)  {
      cat("Q-matrix validation...")
      Qv <- Qval(GDINA1.obj, eps = eps)
      mQ <- show.Qval(Qv,what = "sug.Q")
      if (any(mQ != Q)) {
        GDINA2.obj <-
          GDINA(dat, Q, model = GDINA2.opt$model, sequential = GDINA2.opt$sequential,
                higher.order = GDINA2.opt$higher.order,
                higher.order.model = GDINA2.opt$higher.order.model,
                higher.order.method = GDINA2.opt$higher.order.method,
                higher.order.SE = GDINA2.opt$higher.order.SE,
                mono.constraint = GDINA2.opt$mono.constraint,
                SE = GDINA2.opt$SE,lower.p = GDINA2.opt$lower.p,upper.p = GDINA2.opt$upper.p,
                smallNcorrection = GDINA2.opt$smallNcorrection,
                SE.type = GDINA2.opt$SE.type, conv.type=GDINA2.opt$conv.type,
                empirical = GDINA2.opt$empirical,Mstep.warning=GDINA2.opt$Mstep.warning,
                verbose = GDINA2.opt$verbose,optimizer=GDINA2.opt$optimizer,
                nstarts = GDINA2.opt$nstarts,randomseed=GDINA2.opt$randomseed,
                att.prior = GDINA2.opt$att.prior,
                att.str = GDINA2.opt$att.str,
                higher.order.parm = GDINA2.opt$higher.order.parm,
                digits = GDINA2.opt$digits,
                conv.crit = GDINA2.opt$conv.crit,
                diagnosis = GDINA2.opt$diagnosis,
                maxitr = GDINA2.opt$maxitr)
      } else{
        GDINA2.obj <- GDINA1.obj
      }
      cat("done.\n")
    } else{
      GDINA2.obj <- GDINA1.obj
    }
    Wp <- NULL
    Q <- show.GDINA(GDINA2.obj,what = "Q")
    if (modelselection) {
      cat("Model selection... ")
if (tolower(modelselectionrule)=="ds"){
  Wp <- modelcomp.wald(GDINA2.obj, models = reducedCDM,DS=TRUE)
}else{
  Wp <- modelcomp.wald(GDINA2.obj, models = reducedCDM,DS=FALSE)
}

      ps <- show.modelcomp(Wp,what = "wald.p")

      if (tolower(modelselectionrule) == "simpler") {
        model = apply(ps, 1, function(x) {
          if (max(x[1:2], na.rm = T) > alpha.level) {
            which.max(x[1:2])
          } else if (max(x[3:5], na.rm = T) > alpha.level) {
            which.max(x[3:5]) + 2
          } else{
            return(0)
          }
        })
      } else if(tolower(modelselectionrule) == "largestp"){
        model = apply(ps, 1, function(x) {
          if (max(x, na.rm = T) > alpha.level) {
            which.max(x)
          } else{
            return(0)
          }
        })
      }else if(tolower(modelselectionrule) == "ds"){
        model = apply(cbind(ps,show.modelcomp(Wp,what = "DS")), 1, function(x) {
          if (max(x[1:5], na.rm = T) > alpha.level) {
            sample(which(x[5+which(x[1:5]>alpha.level)]==min(x[5+which(x[1:5]>alpha.level)])),1)
          } else{
            return(0)
          }
        })
        # print(model)
      }
      models = rep(0, nrow(Q))
      models[which(rowSums(Q)!= 1)] = model

    }else{models <- 0}
    cat("done. \n")

    cat("Final calibration... ")

    CDM.obj <-
      GDINA(dat, Q, model = CDM.opt$model, sequential = CDM.opt$sequential,
            higher.order = CDM.opt$higher.order,
            higher.order.model = CDM.opt$higher.order.model,
            higher.order.method = CDM.opt$higher.order.method,
            higher.order.SE = CDM.opt$higher.order.SE,
            mono.constraint = CDM.opt$mono.constraint,
            SE = CDM.opt$SE,lower.p = CDM.opt$lower.p,upper.p = CDM.opt$upper.p,
            smallNcorrection = CDM.opt$smallNcorrection,
            SE.type = CDM.opt$SE.type, conv.type=CDM.opt$conv.type,
            empirical = CDM.opt$empirical,Mstep.warning=CDM.opt$Mstep.warning,
            verbose = CDM.opt$verbose,optimizer=CDM.opt$optimizer,
            nstarts = CDM.opt$nstarts,randomseed=CDM.opt$randomseed,
            att.prior = CDM.opt$att.prior,
            att.str = CDM.opt$att.str,
            higher.order.parm = CDM.opt$higher.order.parm,
            digits = CDM.opt$digits, conv.crit = CDM.opt$conv.crit,
            diagnosis = CDM.opt$diagnosis,
            maxitr = CDM.opt$maxitr)
    cat("done.\n")
    out <-list( CDM.obj = CDM.obj,
        Qval.obj = Qv,
        Wald.obj = Wp,
        GDINA1.obj = GDINA1.obj,
        GDINA2.obj = GDINA2.obj,
        options = list(
          reducedCDM = reducedCDM,
          alpha.level = alpha.level,
          Qvalid = Qvalid, eps = eps,
          modelselection = modelselection,
          modelselectionrule =modelselectionrule,
          GDINA1.option = GDINA1.option,
          GDINA2.option = GDINA2.option,
          SelectedCDM.option = SelectedCDM.option,
          ASEcall = ASEcall
        )
      )
    class(out) = "autoGDINA"
    return(out)
  }
