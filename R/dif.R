#' @include GDINA.R
#' @title  Differential item functioning for cognitive diagnosis models
#'
#' @description   This function is used to detect differential item functioning using the Wald test (Hou, de la Torre, & Nandakumar, 2014; Ma, Terzi, & de la Torre, 2021) and the likelihood ratio
#' test (Ma, Terzi, & de la Torre, 2021). The forward anchor item search procedure developed in Ma, Terzi, and de la Torre (2021) was implemented. Note that it can only detect DIF for two groups currently.
#'
#' @param dat item responses from two groups; missing data need to be coded as \code{NA}
#' @param Q Q-matrix specifying the association between items and attributes
#' @param model model for each item.
#' @param group a factor or a vector indicating the group each individual belongs to. Its length must be equal to the number of individuals.
#' @param method DIF detection method; It can be \code{"wald"} for Hou, de la Torre, and Nandakumar's (2014)
#' Wald test method, and \code{"LR"} for likelihood ratio test (Ma, Terzi, Lee,& de la Torre, 2017).
#' @param p.adjust.methods adjusted p-values for multiple hypothesis tests. This is conducted using \code{p.adjust} function in \pkg{stats},
#'  and therefore all adjustment methods supported by \code{p.adjust} can be used, including \code{"holm"},
#'  \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"} and \code{"BY"}. See \code{p.adjust}
#'  for more details. \code{"holm"} is the default.
#' @param anchor.items which items will be used as anchors? Default is \code{NULL}, which means none of the items are used as anchors.
#'  For LR method, it can also be an integer vector giving the item numbers for anchors or \code{"all"}, which means all items are treated as anchor items.
#' @param dif.items which items are subject to DIF detection? Default is \code{"all"}. It can also be an integer vector giving the item numbers.
#' @param approx Whether an approximated LR test is implemented? If TRUE, parameters of items except the studied one will not be re-estimated.
#' @param SE.type Type of standard error estimation methods for the Wald test.
#' @param FS.args arguments for the forward anchor item search procedure developed in Ma, Terzi, and de la Torre (2021). A list with the following elements:
#'  \itemize{
#'    \item \code{on} - logical; \code{TRUE} if activate the forward anchor item search procedure. Default = \code{FALSE}.
#'    \item \code{alpha.level} - nominal level for Wald or LR test. Default = .05.
#'    \item \code{maxit} - maximum number of iterations allowed. Default = 10.
#'    \item \code{verbose} - logical; print information for each iteration or not? Default = \code{FALSE}.
#'    }
#' @param ... arguments passed to GDINA function for model calibration
#' @return A data frame giving the Wald statistics and associated p-values.
#'
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @seealso \code{\link{GDINA}}
#' @export
#' @examples
#' \dontrun{
#' set.seed(123456)
#' N <- 3000
#' Q <- sim30GDINA$simQ
#' gs <- matrix(.2,ncol = 2, nrow = nrow(Q))
#' # By default, individuals are simulated from uniform distribution
#' # and deltas are simulated randomly
#' sim1 <- simGDINA(N,Q,gs.parm = gs,model="DINA")
#' sim2 <- simGDINA(N,Q,gs.parm = gs,model=c(rep("DINA",nrow(Q)-1),"DINO"))
#' dat <- rbind(extract(sim1,"dat"),extract(sim2,"dat"))
#' gr <- rep(c("G1","G2"),each=N)
#'
#' # DIF using Wald test
#' dif.wald <- dif(dat, Q, group=gr, method = "Wald")
#' dif.wald
#' # DIF using LR test
#' dif.LR <- dif(dat, Q, group=gr, method="LR")
#' dif.LR
#' # DIF using Wald test + forward search algorithm
#' dif.wald.FS <- dif(dat, Q, group=gr, method = "Wald", FS.args = list(on = TRUE, verbose = TRUE))
#' dif.wald.FS
#' # DIF using LR test + forward search algorithm
#' dif.LR.FS <- dif(dat, Q, group=gr, method = "LR", FS.args = list(on = TRUE, verbose = TRUE))
#' dif.LR.FS
#'}
#' @references
#' Hou, L., de la Torre, J., & Nandakumar, R. (2014). Differential item functioning assessment in cognitive diagnostic modeling: Application of the Wald test to
#' investigate DIF in the DINA model. \emph{Journal of Educational Measurement, 51}, 98-125.
#'
#' Ma, W., Terzi, R., & de la Torre, J. (2021). Detecting differential item functioning using multiple-group cognitive diagnosis models. \emph{Applied Psychological Measurement}.
#'


dif <- function(dat, Q, group, model = "GDINA", method = "wald", anchor.items = NULL, dif.items = "all", p.adjust.methods = "holm", approx = FALSE,
                SE.type = 2, FS.args = list(on = FALSE, alpha.level = .05, maxit = 10, verbose = FALSE),...){

  if (!is.matrix(dat))
    dat <- as.matrix(dat)

  if (!is.matrix(Q))
    Q <- as.matrix(Q)

  if (nrow(dat) != length(group))
    stop("The length of group variable must be equal to the number of individuals.", call. = FALSE)

  dat <- dat[order(group),]
  group <- sort(group)

  if(is.factor(group)){
    gr.label <- levels(group)
  }else if(is.vector(group)){
    gr.label <- unique(group)
  }

  if (length(gr.label) != 2)
    stop("Only two group DIF can be examined.", call. = FALSE)

  J <- nrow(Q)

  ### Anchor items
  if (length(anchor.items) == 1 && tolower(anchor.items) == "all")
    anchor.items <- seq_len(J)

  myFS <- list( on = FALSE, alpha.level = .05, maxit = 10, verbose = FALSE  )

  FS.args <- utils::modifyList(myFS,FS.args)


  log.purification <- NULL

  if(tolower(method)=="wald"){

    if(FS.args$on){

      anchor.items <- NULL
      dif.items <- seq_len(J)

      x <- purif.WaldDIF(dat = dat,Q = Q, group = group, model = model, SE.type = SE.type,
                         alpha.level = FS.args$alpha.level, maxit = FS.args$maxit, progress = FS.args$verbose, ...)
      output <- x$output
      log.purification <- x$log
      output <- as.data.frame(output)
      colnames(output) <- c("Wald stat.","df","p.value")
      rownames(output) <- paste("Item",seq_len(J))

      }else{


      ### DIF items => a numeric vector
      if (length(dif.items) == 1 && tolower(dif.items) == "all") {
        dif.items <- seq_len(J)
      } else if(any(!is.numeric(dif.items))){
        stop("dif.items needs to be correctly specified.", call. = FALSE)
      }else if (min(dif.items) <= 0 || max(dif.items) > J){
        stop("dif.items needs to be correctly specified.", call. = FALSE)
      }


      if(all(1:J %in% anchor.items))
          stop("At least one item needs to be non-anchor items when Wald test is used.",call. = FALSE)

      if(any(dif.items %in% anchor.items))
          stop("dif.items must be different from anchor.items.",call. = FALSE)

      nonstudied.items <- NULL

      if(!identical(sort(union(anchor.items,dif.items)),seq_len(J)))
        nonstudied.items <- setdiff(seq_len(J),union(anchor.items,dif.items))

      output <- WaldDIF(dat = dat,Q = Q, group = group, anchor.items = anchor.items, dif.items = dif.items, nonstudied.items = nonstudied.items,
                        model = model, SE.type = SE.type, ...)

      output <- as.data.frame(output)
      colnames(output) <- c("Wald stat.","df","p.value")
      rownames(output) <- paste("Item",dif.items)
      }



  }else if(method=="LR"){

    if(FS.args$on) {

      anchor.items <- NULL
      dif.items <- seq_len(J)

      output <- NULL
      it <- 0

      while(it<FS.args$maxit){
        output <- LRDIF(dat = dat,Q = Q, group = group, model = model, anchor.items = anchor.items, dif.items = dif.items, LR.approx = approx,...)
        it <- it + 1
        log.purification[[it]] <- output
        if(FS.args$verbose){
          cat("Iter = ", it,"Anchor items = Items ",anchor.items,"\n")
          # rownames(output) <- paste("Item",seq_len(J))
          print(output)
        }
        new.anchoritems.loc <- which(output$p.value > FS.args$alpha.level)
        if(length(new.anchoritems.loc)==0)
          new.anchoritems.loc <- NULL
        if(identical(anchor.items,new.anchoritems.loc))
          break

        anchor.items <- new.anchoritems.loc
      }
      # rownames(output) <- paste("Item",seq_len(J))
    }else{

      output <- LRDIF(dat = dat,Q = Q, group = group, model = model, anchor.items = anchor.items, dif.items = dif.items, LR.approx = approx,...)

    }

    # rownames(output) <- extract(est,"item.names")[dif.items]
    # output <- lr.out
  }
  output$'adj.pvalue' <- stats::p.adjust(output$'p.value', method = p.adjust.methods)
output <- list(test=output,group=group,p.adjust.methods=p.adjust.methods,log.purification = log.purification)
class(output) <- "dif"
invisible(output)

}

  WaldDIF <-
    function(dat, Q, group, model, anchor.items, dif.items, nonstudied.items = NULL, SE.type = 2, ...) {


      if(length(model)==1)
        model <- rep(model, ncol(dat))

      m <- model[c(dif.items, dif.items, anchor.items, nonstudied.items, nonstudied.items)]
      JD <- 2 * length(dif.items)
      JA <- length(anchor.items)
      JN <- 2 * length(nonstudied.items)
      J <- JD + JA + JN

      Data <- matrix(0, nrow(dat), J)
      QQ <- matrix(0, J, ncol(Q))

      gr.label <- unique(group)


      Data[, seq_len(JD)] <- GDINA::bdiagMatrix(list(dat[which(group == gr.label[1]), dif.items],
                                                     dat[which(group == gr.label[2]), dif.items]), NA)
      QQ[seq_len(JD), ] <- Q[rep(dif.items, 2), ]


      if (!is.null(anchor.items)) {
        Data[, (JD + 1):(JD + JA)] <- dat[, anchor.items]
        QQ[(JD + 1):(JD + JA), ] <- Q[anchor.items, ]
      }

      if (!is.null(nonstudied.items)) {

        Data[, (JD + JA + 1):J] <- GDINA::bdiagMatrix(list(dat[which(group == gr.label[1]), nonstudied.items],
                                                           dat[which(group == gr.label[2]), nonstudied.items]), NA)
        QQ[(JD + JA + 1):J, ] <- Q[rep(nonstudied.items, 2), ]
      }

      est <- GDINA::GDINA(dat = Data, Q = QQ, group = group, verbose = 0,model = m, ...)

      output <- matrix(0, JD / 2, 3)

      item.parm <- extract(est, "catprob.parm")
      dcov <- extract(est, "delta.cov", SE.type = SE.type)
      for (j in seq_len(JD / 2)) {
        x <- c(extract(est, "delta.parm")[[j]],
               extract(est, "delta.parm")[[j + JD / 2]])
        R <- cbind(diag(length(x) / 2), -1 * diag(length(x) / 2))
        vcov <-
          bdiagMatrix(list(dcov$cov[dcov$index$loc[dcov$index$item == j],
                                    dcov$index$loc[dcov$index$item == j]],
                           dcov$cov[dcov$index$loc[dcov$index$item == j + JD / 2],
                                    dcov$index$loc[dcov$index$item == j + JD / 2]]))
        output[j, 1] <-
          t(R %*% x) %*% MASS::ginv(R %*% vcov %*% t(R)) %*% (R %*% x)
        output[j, 2] <- nrow(R)
        output[j, 3] <- pchisq(output[j, 1], nrow(R), lower.tail = FALSE)
      }

output

    }


  purif.WaldDIF <-
    function(dat, Q, group, model, SE.type = 2, alpha.level = 0.05, maxit = 10, progress = FALSE, ...) {

      anchor.items <- NULL
      it <- 0
      J <- ncol(dat)
      dif.items <- seq_len(J)
      log.purification <- list()
      while(it < maxit){

        if(is.null(anchor.items)){ # it = 0
          output <- WaldDIF(dat = dat,Q = Q, group = group, model=model, SE.type = SE.type, anchor.items = anchor.items,dif.items = dif.items, ...)
        }else{ # it = 1, 2, ...
          output <- matrix(0, J, 3)
          nonanchor <- setdiff(seq_len(J),anchor.items)
          if(length(nonanchor)==0){
            nonanchor <- NULL
          }else{
            output[nonanchor,] <- WaldDIF(dat = dat, Q = Q, group = group, model=model, SE.type = SE.type, anchor.items = anchor.items,dif.items = nonanchor, ...)
          }

          l.anchor <- length(anchor.items)

          if(l.anchor==1){ # single anchor item
            output[anchor.items,] <- log.purification[[1]][anchor.items,]
          }else{
            for(j in seq_len(l.anchor)){
              output[anchor.items[j],] <- WaldDIF(dat = dat, Q = Q, group = group, model=model, SE.type = SE.type,
                                                     anchor.items = anchor.items[-j],dif.items = anchor.items[j],nonstudied.items = nonanchor,...)
              # print(output)
            }
          }


        }

        it <- it + 1
        log.purification[[it]] <- output
        if(progress){
          if(is.null(anchor.items)){
            cat("Iter = ", it,"No anchor items\n")
            print(output)
          }else{
            cat("Iter = ", it,"Anchor items = Items ",anchor.items,"\n")
            print(output)
          }
       }
        new.anchoritems.loc <- which(output[,3] > alpha.level)
        if(length(new.anchoritems.loc)==0)
          new.anchoritems.loc <- NULL
        if(identical(anchor.items,new.anchoritems.loc))
          break

        anchor.items <- new.anchoritems.loc
      }

  list(output = output, log = log.purification)

}

LRDIF <- function(dat, Q, group, model, anchor.items, dif.items, LR.approx = FALSE,...){

  J <- ncol(dat)
  JJ <- seq_len(J)

  est <- NULL

  if(length(model)==1) model <- rep(model, J)
  J <- nrow(Q)

  if (length(dif.items) == 1 && tolower(dif.items) == "all")
      dif.items <- JJ

  JD <- length(dif.items)
  gr.label <- unique(group)
  G1 <- which(group == gr.label[1])
  G2 <- which(group == gr.label[2])

  maxit <- 2000
  output <- data.frame(neg2LL=rep(NA,JD),
                       LRstat=rep(NA,JD),
                       df=rep(NA,JD),
                       'p.value'=rep(NA,JD))
  rownames(output) <- paste("Item",dif.items)

  if(identical(JJ,sort(anchor.items))) {
    # dif item has dif pars; all other items have same pars for two gr
    if(LR.approx)
      maxit <- c(rep(0,J-1),2000,2000)

    #est: all items have the same pars across groups
    est <- GDINA::GDINA(dat, Q, group = group, model = model, verbose = 0,...)
    item.parm <- extract(est,"catprob.parm")
    for (j in seq_len(JD)){

      locj <- setdiff(JJ,dif.items[j]) #item no. except the studied one

      estj <- GDINA::GDINA(dat = cbind(dat[,locj],
                                       bdiagMatrix(list(dat[G1,dif.items[j]],
                                                        dat[G2,dif.items[j]]),NA)),
                           model = model[c(locj,dif.items[j],dif.items[j])],
                           Q = Q[c(locj,dif.items[j],dif.items[j]),],
                           group = group,control=list(maxitr = maxit),verbose = 0,
                           catprob.parm = item.parm[c(locj,dif.items[j],dif.items[j])],
                           att.prior =  t(extract(est,"posterior.prob")),...)
      output$LRstat[j] <-  deviance(est) - deviance(estj)
      output$df[j] <- npar(estj)$`No. of total item parameters` - npar(est)$`No. of total item parameters` # distribution parameters identical
      output$neg2LL[j] <- deviance(estj)
      output$'p.value'[j] <- pchisq(output$LRstat[j],output$df[j],lower.tail = FALSE)
    }
  }else if(is.null(anchor.items)){#dif item has the same par; all other items have dif pars
    # est: unique item parameters for each group on each item
    est <- GDINA::GDINA(dat = bdiagMatrix(list(dat[G1,],
                                               dat[G2,]),NA),
                        model = rep(model, 2), Q = rbind(Q,Q), group = group,verbose = 0, ...)
    item.parm <- extract(est,"catprob.parm")
    if(LR.approx)
      maxit <- c(rep(0,2*J-2),2000)


    for (j in seq_len(length(dif.items))){
      locj <- setdiff(JJ,dif.items[j]) #item no. except the studied one

      estj <- GDINA::GDINA(dat = cbind(bdiagMatrix(list(dat[G1,locj],
                                                        dat[G2,locj]),NA),
                                       dat[,dif.items[j]]),
                           Q = Q[c(locj,locj,dif.items[j]),],
                           model = model[c(locj,locj,dif.items[j])],
                           group = group,control=list(maxitr = maxit),verbose = 0,
                           catprob.parm = item.parm[c(locj,(locj+J),dif.items[j])],
                           att.prior = t(extract(est,"posterior.prob")),...)
      output$LRstat[j] <- deviance(estj) - deviance(est)
      output$df[j] <- npar(est)$`No. of total item parameters` - npar(estj)$`No. of total item parameters`
      output$neg2LL[j] <- deviance(estj)
      output$'p.value'[j] <- pchisq(output$LRstat[j],output$df[j],lower.tail = FALSE)
    }


  }else{
    variant.items <- setdiff(JJ, anchor.items)
    # est: anchor items + non-anchor for g1 + non-anchor for g2
    est <- GDINA::GDINA(dat = cbind(dat[,anchor.items],
                                    bdiagMatrix(list(dat[G1,variant.items],
                                                     dat[G2,variant.items]),
                                                NA)),
                        model = model[c(anchor.items,variant.items,variant.items)],
                        Q = Q[c(anchor.items,variant.items,variant.items),],
                        group = group, verbose = 0,...)
    item.parm <- extract(est,"catprob.parm")
    item.set <- c(anchor.items,variant.items,variant.items)


    for (j in seq_len(length(dif.items))){
      # if the studied item is an anchor item, it will be treated as a non-anchor item
      if(dif.items[j]%in%anchor.items){
        if(LR.approx)
          maxit <- c(rep(0,length(item.set)-1),rep(2000,2))

        locj <- which(dif.items[j]==item.set)
        if(length(anchor.items)==1){
          estj <- GDINA::GDINA(dat = cbind(bdiagMatrix(list(dat[G1,variant.items],
                                                            dat[G2,variant.items]),
                                                       NA),
                                           bdiagMatrix(list(dat[G1,dif.items[j]],
                                                            dat[G2,dif.items[j]]),
                                                       NA)),
                               model = model[c(variant.items,variant.items,
                                               dif.items[j],dif.items[j])],
                               Q = Q[c(variant.items,variant.items,
                                       dif.items[j],dif.items[j]),],
                              group = group, verbose = 0,
                              control=list(maxitr = maxit),
                              catprob.parm = item.parm[c(seq_len(length(item.parm))[-locj],locj,locj)],
                              att.prior = t(extract(est,"posterior.prob")),...)
        }else{
          updated.anchor.items <- anchor.items[-locj]

          estj <- GDINA::GDINA(dat = cbind(dat[,updated.anchor.items],
                                           bdiagMatrix(list(dat[G1,variant.items],
                                                            dat[G2,variant.items]),
                                                       NA),
                                           bdiagMatrix(list(dat[G1,dif.items[j]],
                                                            dat[G2,dif.items[j]]),
                                                       NA)),
                               model = model[c(updated.anchor.items,variant.items,variant.items,
                                               dif.items[j],dif.items[j])],
                               Q = Q[c(updated.anchor.items,variant.items,variant.items,
                                       dif.items[j],dif.items[j]),],
                               group = group,verbose = 0, control=list(maxitr = maxit),
                               catprob.parm = item.parm[c(seq_len(length(item.parm))[-locj],locj,locj)],
                               att.prior = t(extract(est,"posterior.prob")),...)
        }

        output$LRstat[j] <- deviance(est) - deviance(estj)
        output$df[j] <- npar(estj)$`No. of total item parameters` - npar(est)$`No. of total item parameters`
        output$neg2LL[j] <- deviance(estj)
        output$'p.value'[j] <- pchisq(output$LRstat[j],output$df[j],lower.tail = FALSE)
      }else{

        # if the studied item is not an anchor item, it will be treated as an anchor item
        if(LR.approx)
          maxit <- c(rep(0,length(item.set)-2),2000)

        if(length(variant.items)>1){

          locj <- which(variant.items!=dif.items[j])

          estj <- GDINA::GDINA(dat = cbind(dat[,anchor.items],
                                           bdiagMatrix(list(dat[G1,variant.items[locj]],
                                                            dat[G2,variant.items[locj]]),
                                                       NA),
                                           dat[,dif.items[j]]),
                               model = model[c(anchor.items,variant.items[locj],variant.items[locj],dif.items[j])],
                               Q = Q[c(anchor.items,variant.items[locj],variant.items[locj],dif.items[j]),],
                               group = group,verbose = 0, control=list(maxitr = maxit),
                               catprob.parm = item.parm[c(which(item.set!=dif.items[j]),which(item.set==dif.items[j])[1])],
                               att.prior = t(extract(est,"posterior.prob")),...)
          output$LRstat[j] <- deviance(estj) - deviance(est)
          output$df[j] <- npar(est)$`No. of total item parameters` - npar(estj)$`No. of total item parameters`
          output$neg2LL[j] <- deviance(estj)
          output$'p.value'[j] <- pchisq(output$LRstat[j],output$df[j],lower.tail = FALSE)
        }else{ # only one variant item

          estj <- GDINA::GDINA(dat = dat[,c(anchor.items,variant.items)],
                               model = model[c(anchor.items,variant.items)],
                               Q = Q[c(anchor.items,variant.items),],
                               group = group,verbose = 0, control=list(maxitr = maxit),
                               catprob.parm = item.parm[c(which(item.set!=dif.items[j]),which(item.set==dif.items[j])[1])],
                               att.prior = t(extract(est,"posterior.prob")),...)
          output$LRstat[j] <- deviance(estj) - deviance(est)
          output$df[j] <- npar(est)$`No. of total item parameters` - npar(estj)$`No. of total item parameters`
          output$neg2LL[j] <- deviance(estj)
          output$'p.value'[j] <- pchisq(output$LRstat[j],output$df[j],lower.tail = FALSE)
        }

      }

  }


  }
  output
}

