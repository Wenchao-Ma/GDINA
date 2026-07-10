#' @include GDINA.R
#' @title  Differential item functioning for cognitive diagnosis models
#'
#' @description   This function is used to detect differential item functioning using the Wald test (Hou, de la Torre, & Nandakumar, 2014; Ma, Terzi, & de la Torre, 2021) and the likelihood ratio
#' test (Ma, Terzi, & de la Torre, 2021). The forward anchor item search procedure developed in Ma, Terzi, and de la Torre (2021) was implemented.
#'
#' @param dat item responses from two or more groups; missing data need to be coded as \code{NA}
#' @param Q Q-matrix specifying the association between items and attributes
#' @param model model for each item.
#' @param sequential Logical; whether a sequential model is fit to the data. Default is \code{FALSE}.
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
#' @author Wenchao Ma, The University of Minnesota, \email{wma@umn.edu}
#'  Jimmy de la Torre, The University of Hong Kong
#' @seealso \code{\link{GDINA}}
#' @export
#' @examples
#' \dontrun{
#'
#'####################################
#'#          Example 1.              #
#'#        two-group DIF             #
#'####################################
#'
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
#'
#'####################################
#'#          Example 2.              #
#'#        two-group DIF             #
#'#    with attribute structure.     #
#'####################################
#' # --- User-specified attribute structure ----#
#' Q <- sim30GDINA$simQ
#' K <- ncol(Q)
#' # linear structure A1->A2->A3->A4->A5
#' linear <- list(c(1,2),
#'                c(2,3),
#'                c(3,4),
#'                c(4,5))
#' struc <- att.structure(linear,K)
#' set.seed(123)
#' # data simulation
#' N <- 1000
#' true.lc <- sample(c(1:2^K),N,replace=TRUE,prob=struc$att.prob)
#' table(true.lc) #check the sample
#' true.att <- attributepattern(K)[true.lc,]
#'  gs <- matrix(rep(0.1,2*nrow(Q)),ncol=2)
#'  # data simulation
#'  simD <- simGDINA(N,Q,gs.parm = gs, model = "ACDM",attribute = true.att)
#'  dat <- extract(simD,"dat")
#' gr <- rep(1:2,each=N/2)
#' dif.wald <- dif(dat, Q, group=gr, method = "Wald",
#'                 att.str = diverg, att.dist = "saturated")
#' dif.wald
#'
#'
#'####################################
#'#          Example 3.              #
#'#        Three-group DIF           #
#'####################################
#'
#' set.seed(123456)
#' N <- 500
#' Q <- sim30GDINA$simQ
#' gs <- matrix(.2,ncol = 2, nrow = nrow(Q))
#' # By default, individuals are simulated from uniform distribution
#' # and deltas are simulated randomly
#' sim1 <- simGDINA(N,Q,gs.parm = gs,model="DINA")
#' sim2 <- simGDINA(N,Q,gs.parm = gs,model=c(rep("DINA",nrow(Q)-1),"DINO"))
#' sim3 <- simGDINA(N,Q,gs.parm = gs,model="DINA")
#' dat <- rbind(extract(sim1,"dat"),extract(sim2,"dat"),extract(sim3,"dat"))
#' gr <- rep(c("G1","G2","G3"),each=N)
#'
#' # DIF using Wald test - omnibus test
#' dif.wald <- dif(dat, Q, group=gr, method = "Wald", FS.args = list(on = TRUE, verbose = TRUE))
#' dif.wald
#' # pairwise comparison between all pair of groups
#' dif.pair = pairwiseDIF(dif.wald)
#' dif.pair
#' # draw plot to show the DIF
#' plot(dif.pair, withSE = TRUE)
#'}
#' @references
#' Hou, L., de la Torre, J., & Nandakumar, R. (2014). Differential item functioning assessment in cognitive diagnostic modeling: Application of the Wald test to
#' investigate DIF in the DINA model. \emph{Journal of Educational Measurement, 51}, 98-125.
#'
#' Ma, W., Terzi, R., & de la Torre, J. (2021). Detecting differential item functioning using multiple-group cognitive diagnosis models. \emph{Applied Psychological Measurement}.
#'


dif <- function(dat, Q, group, model = "GDINA", sequential = FALSE, method = "wald", anchor.items = NULL, dif.items = "all", p.adjust.methods = "holm", approx = FALSE,
                SE.type = 2, FS.args = list(on = FALSE, alpha.level = .05, maxit = 10, verbose = FALSE),...){

  if(sequential){
    dat <- seq_coding(dat, Q)
    originalQ <- Q
    Q <- Q[,-c(1, 2)]
    all.item.names <- paste("Item", originalQ[, 1], "Cat", originalQ[, 2])
  }else{
    all.item.names <- paste("Item", seq(nrow(Q)))
  }



  if (!is.matrix(dat))
    dat <- as.matrix(dat)

  if (!is.matrix(Q))
    Q <- as.matrix(Q)

  group.info <- normalize_dif_group(group = group, nobs = nrow(dat))
  dat <- dat[group.info$order,,drop = FALSE]
  group <- group.info$group
  gr.label <- group.info$label
  ngroup <- group.info$ngroup

  if (ngroup < 2)
    stop("At least two observed groups are required.", call. = FALSE)

  J <- nrow(Q)

  ### Anchor items
  if (length(anchor.items) == 1 && tolower(anchor.items) == "all")
    anchor.items <- seq_len(J)

  myFS <- list( on = FALSE, alpha.level = .05, maxit = 10, verbose = FALSE  )

  FS.args <- utils::modifyList(myFS,FS.args)


  log.purification <- NULL

  method <- tolower(method)

  if(method=="wald"){

    if(FS.args$on){

      anchor.items <- NULL
      dif.items <- seq_len(J)

      x <- purif.WaldDIF(dat = dat,Q = Q, group = group, model = model, SE.type = SE.type,
                         alpha.level = FS.args$alpha.level, maxit = FS.args$maxit, progress = FS.args$verbose, ...)
      output <- x$output
      log.purification <- x$log
      output <- as.data.frame(output)
      colnames(output) <- c("Wald stat.","df","p.value")
      rownames(output) <- all.item.names

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

      output <- WaldDIF(dat = dat,Q = Q, group = group, gr.label = gr.label,
            anchor.items = anchor.items,
                        dif.items = dif.items, nonstudied.items = nonstudied.items,
                        model = model, SE.type = SE.type, ...)

      output <- as.data.frame(output)
      colnames(output) <- c("Wald stat.","df","p.value")
      rownames(output) <- all.item.names[dif.items]
      }



  }else if(method=="lr"){

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
          cat("Iter = ", it,"Anchor items =  ",all.item.names[anchor.items],"\n")
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
  }else{
    stop("method must be either 'wald' or 'LR'.", call. = FALSE)
  }
  tested.items <- dif.items
  if (length(tested.items) == 1 && is.character(tested.items) && tolower(tested.items) == "all")
    tested.items <- seq_len(J)
  output$'adj.pvalue' <- stats::p.adjust(output$'p.value', method = p.adjust.methods)
output <- list(test=output,group=group,p.adjust.methods=p.adjust.methods,
               log.purification = log.purification,
               input = list(dat = dat,
                            Q = Q,
                            group = group,
                            gr.label = gr.label,
                            ngroup = ngroup,
                            model = model,
                            all.item.names = all.item.names,
                            tested.items = tested.items,
                            sequential = sequential,
                            method = method,
                            fit.args = list(...),
                            SE.type = SE.type,
                            call = match.call()))
class(output) <- "dif"
invisible(output)

}

normalize_dif_group <- function(group, nobs){
  if (nobs != length(group))
    stop("The length of group variable must be equal to the number of individuals.", call. = FALSE)

  ord <- order(group)
  group <- group[ord]

  if(is.factor(group)){
    group <- droplevels(group)
    gr.label <- levels(group)
  }else{
    gr.label <- unique(group)
  }

  list(group = group, order = ord, label = gr.label, ngroup = length(gr.label))
}

build_group_specific_block <- function(dat, items, group, gr.label){
  if(is.null(items) || length(items) == 0)
    return(NULL)

  GDINA::bdiagMatrix(lapply(gr.label, function(g) dat[group == g, items, drop = FALSE]), NA)
}

build_group_specific_index <- function(items, ngroup){
  if(is.null(items) || length(items) == 0)
    return(integer(0))

  rep(items, ngroup)
}

build_lr_item_config <- function(shared.items, variant.items, ngroup,
                                 gr.label = NULL, item.names = NULL){
  shared.items <- if(is.null(shared.items)) integer(0) else shared.items
  variant.items <- if(is.null(variant.items)) integer(0) else variant.items

  config <- data.frame(item = c(shared.items, rep(variant.items, ngroup)),
                       group = c(rep(0L, length(shared.items)),
                                 rep(seq_len(ngroup), each = length(variant.items))))

  config$refit.item <- seq_len(nrow(config))
  config$original.item <- config$item
  config$item.type <- ifelse(config$group == 0L, "shared", "group-specific DIF")

  group.no <- config$group
  group.no[group.no == 0L] <- NA_integer_
  config$group.no <- group.no

  if(is.null(gr.label)){
    config$group.label <- "shared"
    config$group.label[config$group > 0L] <- as.character(config$group[config$group > 0L])
  }else{
    config$group.label <- rep("shared", nrow(config))
    config$group.label[config$group > 0L] <- as.character(gr.label[config$group[config$group > 0L]])
  }

  if(is.null(item.names)){
    config$item.label <- paste("Item", config$item)
  }else{
    config$item.label <- as.character(item.names[config$item])
  }

  config
}

build_lr_model_frame <- function(dat, Q, group, gr.label, model, shared.items, variant.items,
                                 item.names = NULL){
  ngroup <- length(gr.label)
  if(length(model) == 1)
    model <- rep(model, nrow(Q))
  item.index <- c(shared.items, build_group_specific_index(variant.items, ngroup))
  data.parts <- list()

  if(!is.null(shared.items) && length(shared.items) > 0)
    data.parts[[length(data.parts) + 1L]] <- dat[, shared.items, drop = FALSE]

  if(!is.null(variant.items) && length(variant.items) > 0)
    data.parts[[length(data.parts) + 1L]] <- build_group_specific_block(dat = dat, items = variant.items,
                                                                        group = group, gr.label = gr.label)

  Data <- if(length(data.parts) == 1L) {
    data.parts[[1L]]
  }else{
    do.call(cbind, data.parts)
  }

  list(dat = Data,
       Q = Q[item.index, , drop = FALSE],
       model = model[item.index],
       config = build_lr_item_config(shared.items = shared.items,
                                     variant.items = variant.items,
                  ngroup = ngroup,
                  gr.label = gr.label,
                  item.names = item.names))
}

map_lr_init_params <- function(source.item.parm, source.config, target.config){
  if(nrow(target.config) == 0)
    return(list())

  idx <- integer(nrow(target.config))
  for(i in seq_len(nrow(target.config))){
    exact <- which(source.config$item == target.config$item[i] &
                     source.config$group == target.config$group[i])
    if(length(exact) == 0)
      exact <- which(source.config$item == target.config$item[i] & source.config$group == 0L)
    if(length(exact) == 0)
      exact <- which(source.config$item == target.config$item[i])
    idx[i] <- exact[1]
  }

  source.item.parm[idx]
}

fit_lr_model <- function(dat, Q, group, gr.label, model, shared.items, variant.items,
          item.names = NULL,
                         init.parm = NULL, att.prior = NULL, control.maxitr = NULL, ...){
  frame <- build_lr_model_frame(dat = dat, Q = Q, group = group, gr.label = gr.label,
                                model = model, shared.items = shared.items,
            variant.items = variant.items,
            item.names = item.names)

  fit.args <- list(dat = frame$dat,
                   Q = frame$Q,
                   group = group,
                   model = frame$model,
                   verbose = 0)

  if(!is.null(control.maxitr))
    fit.args$control <- list(maxitr = control.maxitr)
  if(!is.null(init.parm))
    fit.args$catprob.parm <- init.parm
  if(!is.null(att.prior))
    fit.args$att.prior <- att.prior

  fit.args <- c(fit.args, list(...))

  list(est = do.call(GDINA::GDINA, fit.args), config = frame$config)
}

build_pairwise_wald_table <- function(est, config, dif.items, gr.label, item.names,
                                      SE.type = 2, p.adjust.methods = "holm"){
  if(length(dif.items) == 0)
    return(data.frame(item = character(0), group1 = character(0), group2 = character(0),
                      check.names = FALSE))

  delta.parm <- extract(est, "delta.parm")
  dcov <- extract(est, "delta.cov", SE.type = SE.type)
  group.pairs <- utils::combn(seq_along(gr.label), 2, simplify = FALSE)
  pairwise.output <- vector("list", length(dif.items) * length(group.pairs))
  idx <- 0L

  for(dif.item in dif.items){
    item.loc <- which(config$item == dif.item & config$group > 0L)
    item.loc <- item.loc[order(config$group[item.loc])]
    x <- unlist(delta.parm[item.loc], use.names = FALSE)
    nparm <- length(delta.parm[[item.loc[1]]])
    vcov <- bdiagMatrix(lapply(item.loc, function(loc) {
      dcov$cov[dcov$index$loc[dcov$index$item == loc],
               dcov$index$loc[dcov$index$item == loc], drop = FALSE]
    }))

    for(pair in group.pairs){
      idx <- idx + 1L
      R <- matrix(0, nrow = nparm, ncol = length(x))
      loc1 <- ((pair[1] - 1L) * nparm + 1L):(pair[1] * nparm)
      loc2 <- ((pair[2] - 1L) * nparm + 1L):(pair[2] * nparm)
      R[, loc1] <- diag(nparm)
      R[, loc2] <- -diag(nparm)
      stat <- as.numeric(t(R %*% x) %*% MASS::ginv(R %*% vcov %*% t(R)) %*% (R %*% x))
      pairwise.output[[idx]] <- data.frame(item = item.names[dif.item],
                                           group1 = as.character(gr.label[pair[1]]),
                                           group2 = as.character(gr.label[pair[2]]),
                                           'Wald stat.' = stat,
                                           df = nrow(R),
                                           p.value = pchisq(stat, nrow(R), lower.tail = FALSE),
                                           check.names = FALSE)
    }
  }

  pairwise.output <- do.call(rbind, pairwise.output)
  pairwise.output$adj.pvalue <- unsplit(lapply(split(pairwise.output$p.value, pairwise.output$item),
                                              stats::p.adjust,
                                              method = p.adjust.methods),
                                        pairwise.output$item)
  pairwise.output
}



  WaldDIF <-
    function(dat, Q, group, gr.label = unique(group), model, anchor.items, dif.items,
             nonstudied.items = NULL, SE.type = 2, ...) {


      if(length(model)==1)
        model <- rep(model, ncol(dat))

      ngroup <- length(gr.label)
      ndif <- length(dif.items)

      m <- model[c(build_group_specific_index(dif.items, ngroup),
                   anchor.items,
                   build_group_specific_index(nonstudied.items, ngroup))]
      JD <- ngroup * ndif
      JA <- length(anchor.items)
      JN <- ngroup * length(nonstudied.items)
      J <- JD + JA + JN

      Data <- matrix(0, nrow(dat), J)
      QQ <- matrix(0, J, ncol(Q))


      Data[, seq_len(JD)] <- build_group_specific_block(dat = dat, items = dif.items,
                                                        group = group, gr.label = gr.label)
      QQ[seq_len(JD), ] <- Q[build_group_specific_index(dif.items, ngroup), , drop = FALSE]


      if (!is.null(anchor.items)) {
        Data[, (JD + 1):(JD + JA)] <- dat[, anchor.items, drop = FALSE]
        QQ[(JD + 1):(JD + JA), ] <- Q[anchor.items, , drop = FALSE]
      }

      if (!is.null(nonstudied.items)) {

        Data[, (JD + JA + 1):J] <- build_group_specific_block(dat = dat, items = nonstudied.items,
                                                              group = group, gr.label = gr.label)
        QQ[(JD + JA + 1):J, ] <- Q[build_group_specific_index(nonstudied.items, ngroup), , drop = FALSE]
      }

      est <- GDINA::GDINA(dat = Data, Q = QQ, group = group, verbose = 0,model = m, ...)

      output <- matrix(0, ndif, 3)

      delta.parm <- extract(est, "delta.parm")
      dcov <- extract(est, "delta.cov", SE.type = SE.type)
      for (j in seq_len(ndif)) {
        item.loc <- j + ndif * (seq_len(ngroup) - 1)
        x <- unlist(delta.parm[item.loc], use.names = FALSE)
        nparm <- length(delta.parm[[item.loc[1]]])
        R <- cbind(kronecker(matrix(-1, nrow = ngroup - 1, ncol = 1), diag(nparm)),
                   kronecker(diag(ngroup - 1), diag(nparm)))
        vcov <- bdiagMatrix(lapply(item.loc, function(loc) {
          dcov$cov[dcov$index$loc[dcov$index$item == loc],
                   dcov$index$loc[dcov$index$item == loc], drop = FALSE]
        }))
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

  if(length(model)==1) model <- rep(model, J)
  J <- nrow(Q)

  if (length(dif.items) == 1 && tolower(dif.items) == "all")
      dif.items <- JJ

  JD <- length(dif.items)
  gr.label <- unique(group)

  output <- data.frame(neg2LL=rep(NA,JD),
                       LRstat=rep(NA,JD),
                       df=rep(NA,JD),
                       'p.value'=rep(NA,JD))
  rownames(output) <- paste("Item",dif.items)

  if(identical(JJ,sort(anchor.items))) {
    shared.items <- JJ
    variant.items <- integer(0)
  }else if(is.null(anchor.items)){
    shared.items <- integer(0)
    variant.items <- JJ
  }else{
    shared.items <- sort(anchor.items)
    variant.items <- setdiff(JJ, anchor.items)
  }

  base.fit <- fit_lr_model(dat = dat, Q = Q, group = group, gr.label = gr.label,
                           model = model, shared.items = shared.items,
                           variant.items = variant.items, ...)
  base.est <- base.fit$est
  base.config <- base.fit$config
  base.item.parm <- extract(base.est,"catprob.parm")
  base.att.prior <- t(extract(base.est,"posterior.prob"))
  base.npar <- npar(base.est)$`No. of total item parameters`
  base.deviance <- deviance(base.est)

  for (j in seq_len(JD)){
    dif.item <- dif.items[j]
    if(dif.item %in% shared.items){
      target.shared <- setdiff(shared.items, dif.item)
      target.variant <- sort(unique(c(variant.items, dif.item)))
      unrestricted <- TRUE
    }else{
      target.shared <- sort(unique(c(shared.items, dif.item)))
      target.variant <- setdiff(variant.items, dif.item)
      unrestricted <- FALSE
    }

    target.frame <- build_lr_model_frame(dat = dat, Q = Q, group = group, gr.label = gr.label,
                                         model = model, shared.items = target.shared,
                                         variant.items = target.variant)
    target.init <- map_lr_init_params(source.item.parm = base.item.parm,
                                      source.config = base.config,
                                      target.config = target.frame$config)
    target.maxit <- NULL
    if(LR.approx)
      target.maxit <- ifelse(target.frame$config$item == dif.item, 2000, 0)

    target.fit <- fit_lr_model(dat = dat, Q = Q, group = group, gr.label = gr.label,
                               model = model, shared.items = target.shared,
                               variant.items = target.variant, init.parm = target.init,
                               att.prior = base.att.prior, control.maxitr = target.maxit, ...)

    target.est <- target.fit$est
    target.npar <- npar(target.est)$`No. of total item parameters`
    target.deviance <- deviance(target.est)

    if(unrestricted){
      output$LRstat[j] <- base.deviance - target.deviance
      output$df[j] <- target.npar - base.npar
    }else{
      output$LRstat[j] <- target.deviance - base.deviance
      output$df[j] <- base.npar - target.npar
    }
    output$neg2LL[j] <- target.deviance
    output$'p.value'[j] <- pchisq(output$LRstat[j],output$df[j],lower.tail = FALSE)
  }
  output
}

