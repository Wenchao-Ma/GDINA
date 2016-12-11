#' Q-matrix validation, model selection and calibration in one run
#'
#' \code{autoGDINA} conducts a series of CDM analyses within the G-DINA framework. Particularly,
#' the GDINA model is fitted to the data first using the \code{\link{GDINA}} function;
#' then, the Q-matrix is validated using the function \code{\link{Qval}}.
#' Based on the suggested Q-matrix, the data is fitted by the G-DINA model again, followed
#' by an item level model selection via the Wald test using \code{\link{modelcomp}}. Lastly,
#' the selected models are calibrated based on the suggested Q-matrix using the \code{\link{GDINA}} function.
#' The Q-matrix validation and item-level model selection can be disabled by the users.
#' Possible reduced CDMs for Wald test include the DINA model, the DINO model, A-CDM, LLM and RRUM.
#' See \code{Details} for the rules of item-level model selection.
#'
#'
#' @details
#'
#' After the Wald statistics for each reduced CDM were calculated for each item, the
#' reduced models with p values less than the pre-specified alpha level were rejected.
#' If all reduced models were rejected for an item, the G-DINA model was used as the best model;
#' if at least one reduced model was retained, three diferent rules can be implemented for selecting
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
#'  The reduced model with non-significant p-values but the smallest dissimilarity index is selected as the most appropriate model.
#'  Dissimilarity index can be viewed as an effect size measure, which quatifies how dis-similar the reduced model is from the
#'  G-DINA model. See Ma, Iaconangelo, and de la Torre (2016).
#'
#' @note Returned \code{GDINA1.obj}, \code{GDINA2.obj} and \code{CDM.obj} are objects of class \code{GDINA},
#'  and all S3 methods suitable for \code{GDINA} objects can be applied. See \code{\link{GDINA}} and \code{\link{extract}}.
#'  Similarly, returned \code{Qval.obj} and \code{Wald.obj} are objects of class \code{\link{Qval}} and \code{\link{modelcomp}}.
#'
#' @param Qvalid logical; validate Q-matrix or not? \code{TRUE} is the default.
#' @param alpha.level nominal level for the Wald test. The default is 0.05.
#' @param modelselection logical; conducting model selection or not?
#' @param modelselectionrule how to conducted model selection? Possible options include
#'    \code{simpler}, \code{largestp} and \code{DS}. See \code{Details}.
#' @param eps cut-off value for PVAF if \code{Qvalid=TRUE}. The default is 0.95.
#' @param reducedCDM a vector specifying which reduced CDMs are possible reduced CDMs for each
#'   item. The default is "DINA","DINO","ACDM","LLM",and "RRUM".
#' @param GDINA1.option options for initial G-DINA calibration
#' @param GDINA2.option options for second G-DINA calibration
#' @param CDM.option options for final calibration
#' @inheritParams GDINA
#' @seealso \code{\link{GDINA}}, \code{\link{modelcomp}}, \code{\link{Qval}}
#'
#' @return a list consisting of the following elements:
#' \describe{
#' \item{GDINA1.obj}{initial GDINA calibration of class \code{GDINA}}
#' \item{GDINA2.obj}{second GDINA calibration of class \code{GDINA}}
#' \item{Qval.obj}{Q validation object of class \code{Qval}}
#' \item{Wald.obj}{model comparison object of class \code{modelcomp}}
#' \item{CDM.obj}{Final CDM calibration of class \code{GDINA}}
#' }
#'
#' @author {Wenchao Ma, Rutgers University, \email{wenchao.ma@@rutgers.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#'
#' @references
#' Ma, W., Iaconangelo, C., & de la Torre, J. (2016). Model similarity, model selection and attribute classification.
#' \emph{Applied Psychological Measurement, 40}, 200-217.
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
#' AIC(out1$CDM.obj)
#'
#'
#'
#' #using the other selection rule
#' out11 <- autoGDINA(dat,misQ,modelselectionrule="simpler",reducedCDM = c("DINO","DINA"))
#' out11
#' summary(out11)
#'
#' # disable model selection function
#' out12 <- autoGDINA(dat,misQ,modelselection=FALSE)
#' out12
#' summary(out12)
#'
#' # -- Only consider some reduced CDMs
#' out2 <- autoGDINA(dat,misQ,reducedCDM = c("RRUM","LLM"))
#'
#' # Disable Q-matrix validation
#' out3 <- autoGDINA(dat = dat, Q = misQ, Qvalid = FALSE,alpha.level=0.01)
#' out3
#' summary(out3)
#'
#' @export

autoGDINA <-
  function(dat, Q, modelselection = TRUE, Qvalid = TRUE,
           reducedCDM = c("DINA", "DINO", "ACDM", "LLM", "RRUM"),
           alpha.level = 0.05, modelselectionrule = "simpler",
           eps = 0.95, GDINA1.option = list(),
           GDINA2.option = list(), CDM.option = list()) {
    ASEcall <- match.call()
    # options(warn = -1)
    CDM.opt <- GDINA2.opt <- GDINA1.opt <-
      list(model = "GDINA", sequential = FALSE, item.names = NULL,
           higher.order = FALSE, higher.order.model ="Rasch",
           higher.order.method = "MMLE",
           verbose = 0, catprob.parm = NULL, higher.order.struc.parm = NULL,
           mono.constraint = FALSE,
           empirical = TRUE, att.prior = NULL, att.str = FALSE,
           lower.p = 0.0001,upper.p = 0.9999, smallNcorrection = c(0.0005,0.001),
           nstarts = 1, conv.crit = 0.001, conv.type = "max.p.change",maxitr = 1000,
           diagnosis = 0,Mstep.warning = FALSE,optimizer = "all",
           randomseed = 123456)
    if(GDINA1.opt$sequential) stop("autoGDINA is not available for sequential models.",call. = FALSE)
    GDINA1.opt[names(GDINA1.option)] <- GDINA1.option
    # print(GDINA1.opt)
    GDINA2.opt[names(GDINA2.option)] <- GDINA2.option
    CDM.opt[names(CDM.option)] <- CDM.option
    if (GDINA1.opt$higher.order == TRUE) GDINA1.opt$empirical <- FALSE
    if (GDINA2.opt$higher.order == TRUE) GDINA2.opt$empirical <- FALSE
    if (CDM.opt$higher.order == TRUE) CDM.opt$empirical <- FALSE

    cat("Initial calibration...")
    if (length(GDINA1.option)>0){
      GDINA1.opt$dat <- dat
      GDINA1.opt$Q <- Q
      GDINA1.obj <- do.call("GDINA",GDINA1.opt)
    }else{
      GDINA1.obj <- GDINA(dat = dat, verbose = 0, Q = Q)
    }

    Q <- extract.GDINA(GDINA1.obj,what = "Q")
    cat("done.\n")
    Qv <- NULL
    if (Qvalid)  {
      cat("Q-matrix validation...")
      Qv <- Qval(GDINA1.obj, eps = eps)
      mQ <- extract.Qval(Qv,what = "sug.Q")
      if (any(mQ != Q)) {
        if (length(GDINA2.option)>0){
          GDINA2.opt$dat <- dat
          GDINA2.opt$Q <- mQ
          GDINA2.obj <- do.call("GDINA",GDINA2.opt)

        }else{
          GDINA2.obj <- GDINA(dat = dat, verbose = 0, Q = mQ)
        }
      } else{
        GDINA2.obj <- GDINA1.obj
      }
      cat("done.\n")
    } else{
      GDINA2.obj <- GDINA1.obj
    }
    Wp <- NULL
    Q <- extract.GDINA(GDINA2.obj,what = "Q")
    if (modelselection) {
      cat("Model selection... ")
if (tolower(modelselectionrule)=="ds"){
  Wp <- modelcomp(GDINA2.obj, models = reducedCDM,DS=TRUE)
}else{
  Wp <- modelcomp(GDINA2.obj, models = reducedCDM,DS=FALSE)
}

      ps <- extract.modelcomp(Wp,what = "wald.p")

      if (tolower(modelselectionrule) == "simpler") {

        if(any(toupper(reducedCDM)%in%c("DINA","DINO"))){
          if(any(toupper(reducedCDM)%in%c("LLM","RRUM","ACDM"))){
            #some DINO and DINA; some acdm
        model = apply(ps, 1, function(x) {
          if (max(x[1:2], na.rm = T) > alpha.level) {
            which.max(x[1:2])
          } else if (max(x[3:5], na.rm = T) > alpha.level) {
            which.max(x[3:5]) + 2
          } else{
            return(0)
          }
        })
          }else{
            #no acdm
            model <- apply(ps,1,which.max)
            model[apply(ps,1,max,na.rm=TRUE)<alpha.level] <- 0
          }
        }else{
          # no dina and dino
          model <- apply(ps,1,which.max)
          model[apply(ps,1,max,na.rm=TRUE)<alpha.level] <- 0
        }
      } else if(tolower(modelselectionrule) == "largestp"){
        model = apply(ps, 1, function(x) {
          if (max(x, na.rm = T) > alpha.level) {
            which.max(x)
          } else{
            return(0)
          }
        })
      }else if(tolower(modelselectionrule) == "ds"){
        model = apply(cbind(ps,extract.modelcomp(Wp,what = "DS")), 1, function(x) {
          if (max(x[1:5], na.rm = T) > alpha.level) {
            which.min.randomtie(x[5+which(x[1:5]>alpha.level)])
          } else{
            return(0)
          }
        })
        # print(model)
      }
      models  <-  rep(0, nrow(Q))
      models[which(rowSums(Q)!= 1)] <- model
      cat("done. \n")
    }else{models <- 0}


    cat("Final calibration... ")
    if (length(CDM.option)>0){
      CDM.opt$dat <- dat
      CDM.opt$Q <- Q
      CDM.opt$model <- models
      CDM.obj <- do.call("GDINA",CDM.opt)
    }else{
      CDM.obj <- GDINA(dat = dat, Q = Q, verbose = 0, model = models)
    }
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
          CDM.option = CDM.option,
          ASEcall = ASEcall
        )
      )
    class(out) = "autoGDINA"
    return(out)
  }
