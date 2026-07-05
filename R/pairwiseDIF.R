#' @title Pairwise post hoc DIF analysis
#'
#' @description Conduct pairwise Wald follow-up tests after omnibus DIF analysis
#' with three or more groups.
#' The input can be a \code{dif} object returned by \code{dif()}, in which case
#' DIF items are selected using adjusted omnibus p-values, or raw data plus
#' user-specified DIF items and optional anchor items. \code{plot()} function
#' can be used to draw grouped bar chart for DIF items.
#'
#' @param object an object returned by \code{dif()}.
#' @param dat item responses from two or more groups; missing data need to be coded as \code{NA}.
#' @param Q Q-matrix specifying the association between items and attributes.
#' @param group a factor or a vector indicating the group each individual belongs to.
#' @param model model for each item.
#' @param sequential Logical; whether a sequential model is fit to the data. Default is \code{FALSE}.
#' @param dif.items which items are subject to pairwise DIF follow-up.
#' @param anchor.items optional anchor items. If omitted, all non-DIF items are treated as shared across groups.
#' @param alpha.level adjusted omnibus p-value cutoff used to flag DIF items when \code{object} is supplied.
#' @param p.adjust.methods adjusted p-values for pairwise Wald tests within each item.
#' @param SE.type Type of standard error estimation methods for the Wald test.
#' @param ... arguments passed to \code{GDINA()} for the final post hoc refit.
#' @return A \code{pairwiseDIF} object containing the pairwise Wald test table and final post hoc fit.
#' @seealso \code{\link{dif}}
#'
#' @author Wenchao Ma, The University of Minnesota, \email{wma@umn.edu}
#'
#' @export
#' @examples
#' \dontrun{
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
#' dif.wald <- dif(dat, Q, group=gr, method = "Wald")
#' dif.wald
#' # pairwise comparison between all pair of groups
#' dif.pair = pairwiseDIF(dif.wald)
#' dif.pair
#' # draw plot to show the DIF
#' plot(dif.pair, withSE = TRUE)
#'}
pairwiseDIF <- function(object = NULL, dat = NULL, Q = NULL, group = NULL, model = "GDINA",
                        sequential = FALSE, dif.items = NULL, anchor.items = NULL,
                        alpha.level = 0.05, p.adjust.methods = "holm", SE.type = NULL, ...){

  object.mode <- inherits(object, "dif")

  if(object.mode){
    if(is.null(object$input))
      stop("The supplied dif object does not contain the stored inputs needed for pairwiseDIF. Refit dif() with the current version first.", call. = FALSE)

    dat <- object$input$dat
    Q <- object$input$Q
    group <- object$input$group
    gr.label <- object$input$gr.label
    ngroup <- object$input$ngroup
    model <- object$input$model
    all.item.names <- object$input$all.item.names
    sequential <- object$input$sequential
    J <- nrow(Q)
    if(is.null(SE.type))
      SE.type <- object$input$SE.type
    fit.args <- utils::modifyList(object$input$fit.args, list(...))

    if(is.null(dif.items)){
      tested.items <- object$input$tested.items
      dif.items <- tested.items[object$test$adj.pvalue <= alpha.level]
    }
  }else{
    if(is.null(dat) || is.null(Q) || is.null(group))
      stop("dat, Q, and group must be supplied when object is not a dif object.", call. = FALSE)

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
    J <- nrow(Q)

    if(is.null(SE.type))
      SE.type <- 2
    fit.args <- list(...)
  }

  if(ngroup < 2)
    stop("At least two observed groups are required.", call. = FALSE)

  if(is.null(dif.items))
    stop("dif.items must be supplied when object is not a dif object.", call. = FALSE)

  if (length(dif.items) == 1 && is.character(dif.items) && tolower(dif.items) == "all")
    dif.items <- seq_len(J)

  if(any(!is.numeric(dif.items)) || min(dif.items) <= 0 || max(dif.items) > J)
    stop("dif.items needs to be correctly specified.", call. = FALSE)

  dif.items <- sort(unique(as.integer(dif.items)))

  if(length(dif.items) == 0){
    warning("No items were flagged for DIF for pairwise follow-up.", call. = FALSE)
    output <- list(test = data.frame(), dif.items = integer(0), anchor.items = seq_len(J),
                   group.labels = gr.label, p.adjust.methods = p.adjust.methods,
                   posthoc.fit = NULL, posthoc.config = NULL,
                   item.names = all.item.names, sequential = sequential,
                   source = if(object.mode) "dif" else "manual")
    class(output) <- "pairwiseDIF"
    return(invisible(output))
  }

  if(length(anchor.items) == 1 && is.character(anchor.items) && tolower(anchor.items) == "all")
    anchor.items <- setdiff(seq_len(J), dif.items)

  if(!is.null(anchor.items)){
    if(any(!is.numeric(anchor.items)) || min(anchor.items) <= 0 || max(anchor.items) > J)
      stop("anchor.items needs to be correctly specified.", call. = FALSE)
    anchor.items <- sort(unique(as.integer(anchor.items)))
    if(any(anchor.items %in% dif.items))
      stop("anchor.items must be different from dif.items.", call. = FALSE)
  }

  shared.items <- sort(unique(c(setdiff(seq_len(J), dif.items), anchor.items)))

  fit.call <- c(list(dat = dat, Q = Q, group = group, gr.label = gr.label,
                     model = model, shared.items = shared.items,
                     variant.items = dif.items),
                fit.args)
  posthoc.fit <- do.call(fit_lr_model, fit.call)

  output <- list(test = build_pairwise_wald_table(est = posthoc.fit$est,
                                                  config = posthoc.fit$config,
                                                  dif.items = dif.items,
                                                  gr.label = gr.label,
                                                  item.names = all.item.names,
                                                  SE.type = SE.type,
                                                  p.adjust.methods = p.adjust.methods),
                 dif.items = dif.items,
                 anchor.items = shared.items,
                 group.labels = gr.label,
                 p.adjust.methods = p.adjust.methods,
                 posthoc.fit = posthoc.fit$est,
                 posthoc.config = posthoc.fit$config,
                 item.names = all.item.names,
                 sequential = sequential,
                 source = if(object.mode) "dif" else "manual")
  class(output) <- "pairwiseDIF"
  invisible(output)
}
