#' @include GDINA.R
#' @title  Create plots for GDINA estimates
#'
#' @description   Create various plots for GDINA estimates
#'
#' @param x model object of class \code{\link{GDINA}}
#' @param what type of plot. Can be \code{"IRF"} for item/category response function plot,
#'      or \code{"mp"} for mastery probabilities for individuals.
#' @param item A scalar or vector specifying the item(s) for IRF plots.
#' @param withSE logical; Add error bar (estimate - SE, estimate + SE) to the IRF plots?
#' @param SE.type How is SE estimated. By default, it's based on OPG using incomplete information.
#' @param person A scalar or vector specifying the number of individuals for mastery plots.
#' @param att.names Optional; a vector for attribute names.
#' @param ... additional arguments
#' @seealso \code{\link{GDINA}}, \code{\link{autoGDINA}}
#' @export
#' @examples
#' \dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' #plot item response functions for item 10
#' plot(mod1, item = 10)
#' plot(mod1, what = "IRF", item = 10,withSE = TRUE)
#'
#' # plot mastery probabilities for individuals 4 and 10
#' plot(mod1, what = "mp", person = c(4,10))
#' plot(mod1, what = "mp", person = c(4,10,15),
#' att.names = c("addition","subtraction","multiplication"))
#'}

#' @export
plot.GDINA <-
  function(x, what = "IRF", item = "all", withSE = FALSE, SE.type = 2,
           person = 1, att.names = NULL,...)
  {
    stopifnot(isa(x,"GDINA"))

    if(toupper(what)=="IRF"){
      if(extract(x,"sequential")){
        tit <- "Processing functions"
      }else{
        tit <- "Item success probabilities"
      }

      lc <- p <- upper <- lower <- NULL
      if (withSE) se <- extract(x,what = "catprob.se",SE.type=SE.type)
      ip <- extract(x,what = "catprob.parm")
      if(length(item) == 1){
        if(item == "all"){
          item <- seq_len(length(ip))
        }else{
          item <- item
        }
      } else{
          item <- item
        }
      for (j in item){
        tmp.obj <- ip[[j]]

        tmp.name <- gsub("P\\(","",names(tmp.obj))
        tmp.name <- gsub("\\)","",tmp.name)
        # names(tmp.obj) <- tmp.name
        lc <- factor(tmp.name,levels = tmp.name)

        if(withSE){
          lower=tmp.obj-se[[j]]
          lower[lower<0] <- 0
          upper=tmp.obj+se[[j]]
          upper[upper>1] <- 1
          dat <- data.frame(lc = lc,p = tmp.obj,lower=lower,upper=upper)
          print(ggplot2::ggplot(data = dat, aes(x = lc, y = p)) +
                  geom_bar(stat = "identity", position = "dodge") +
                  geom_errorbar(aes(ymax=upper,ymin=lower), position = "dodge", width = 0.15) +
                  ylim(0,1)+
                  labs(x = "Latent group", y = "Probability of success",
                       title = paste(tit,"for", extract(x,"item.names")[j])))
        }else{
          dat <- data.frame(lc = lc,p = tmp.obj)
          print(ggplot2::ggplot(data = dat, aes(x = lc, y = p)) +
                  geom_bar(stat = "identity", position = "dodge") +
                  ylim(0,1)+
                  labs(x = "Latent group", y = "Probability of success",
                       title = paste(tit,"for", extract(x,"item.names")[j])))
        }



      }
    }else if(tolower(what)=="mp"){
      df <- personparm(x,"mp")
      if(is.null(att.names)){
        att.names <- colnames(df)
      }else{
        att.names <- att.names
      }
      np <- length(person)
      if(np>1){
        dff <- c(t(df[person,]))
        dat <- data.frame(att = rep(att.names,np),mp = dff,person = factor(rep(person,each = ncol(df))))
        print(ggplot2::ggplot(data = dat, ggplot2::aes_string(x = "att", y = "mp")) +
                geom_bar(stat = "identity", position = "dodge",ggplot2::aes_string(fill = "person")) +
                ylim(0,1)+
                labs(x = "Attribute", y = "Mastery probability",
                     title = paste("Mastery probability")))
      }else{
        dff <- c(df[person,])
        dat <- data.frame(att = att.names,mp = dff,person = factor(rep(person,ncol(df))))
        print(ggplot2::ggplot(data = dat, ggplot2::aes_string(x = "att", y = "mp")) +
                geom_bar(stat = "identity", position = "dodge") +
                ylim(0,1)+
                labs(x = "Attribute", y = "Mastery probability",
                     title = paste("Mastery probability for individual",person)))
      }

    }


  }

#' Grouped bar plots for pairwise DIF posthoc analysis
#'
#' Create grouped bar charts of group-specific item/category response probabilities
#' for DIF items identified by \code{pairwiseDIF()}.
#'
#' @param x model object of class \code{\link{pairwiseDIF}}
#' @param item A scalar or vector specifying the DIF item(s) to plot.
#' @param withSE logical; Add error bar (estimate - SE, estimate + SE) to the plots?
#' @param SE.type How is SE estimated. By default, it's based on OPG using incomplete information.
#' @param ... additional arguments
#' @seealso \code{\link{pairwiseDIF}}, \code{\link{dif}}
#' @export
plot.pairwiseDIF <- function(x, item = "all", withSE = FALSE, SE.type = 2, ...){
  lc <- p <- upper <- lower <- group <- NULL

  stopifnot(inherits(x, "pairwiseDIF"))

  if(length(x$dif.items) == 0L || is.null(x$posthoc.fit))
    stop("No DIF items are available for plotting.", call. = FALSE)

  if(is.null(x$posthoc.config) || is.null(x$item.names))
    stop("The supplied pairwiseDIF object does not contain the stored plotting metadata. Refit pairwiseDIF() with the current version first.", call. = FALSE)

  if(length(item) == 1L && is.character(item) && tolower(item) == "all")
    item <- x$dif.items

  if(any(!is.numeric(item)) || any(!item %in% x$dif.items))
    stop("item must be 'all' or a numeric vector of DIF items in the pairwiseDIF object.", call. = FALSE)

  item <- unique(as.integer(item))

  if(isTRUE(x$sequential)){
    tit <- "Group-specific processing functions"
    ylab <- "Probability"
  }else{
    tit <- "Group-specific item success probabilities"
    ylab <- "Probability of success"
  }

  ip <- extract(x$posthoc.fit, what = "catprob.parm")
  if(withSE)
    se <- extract(x$posthoc.fit, what = "catprob.se", SE.type = SE.type)

  for (j in item){
    item.loc <- which(x$posthoc.config$item == j & x$posthoc.config$group > 0L)
    item.loc <- item.loc[order(x$posthoc.config$group[item.loc])]
    plot.dat <- vector("list", length(item.loc))

    for(k in seq_along(item.loc)){
      loc <- item.loc[k]
      tmp.obj <- ip[[loc]]
      tmp.name <- gsub("P\\(", "", names(tmp.obj))
      tmp.name <- gsub("\\)", "", tmp.name)
      tmp.group <- as.character(x$group.labels[x$posthoc.config$group[loc]])

      if(withSE){
        lower <- tmp.obj - se[[loc]]
        lower[lower < 0] <- 0
        upper <- tmp.obj + se[[loc]]
        upper[upper > 1] <- 1
        plot.dat[[k]] <- data.frame(lc = tmp.name, p = tmp.obj, lower = lower,
                                    upper = upper, group = tmp.group)
      }else{
        plot.dat[[k]] <- data.frame(lc = tmp.name, p = tmp.obj, group = tmp.group)
      }
    }

    plot.dat <- do.call(rbind, plot.dat)
    plot.dat$lc <- factor(plot.dat$lc, levels = unique(plot.dat$lc))
    plot.dat$group <- factor(plot.dat$group, levels = as.character(x$group.labels))

    g <- ggplot2::ggplot(data = plot.dat, ggplot2::aes(x = lc, y = p, fill = group)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::ylim(0, 1) +
      ggplot2::labs(x = "Latent group", y = ylab,
                    fill = "Group",
                    title = paste(tit, "for", x$item.names[j]))

    if(withSE){
      g <- g + ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
                                      position = ggplot2::position_dodge(width = 0.9),
                                      width = 0.15)
    }

    print(g)
  }

  invisible(x)
}


#' Item fit plots
#'
#' Create plots of bivariate heatmap for item fit
#'
#' @param x model object of class \code{itemfit}
#' @param type type of heatmap plot
#' @param adjusted logical; plot adjusted or unadjusted p-values?
#' @param ... additional arguments
#' @seealso \code{\link{GDINA}}, \code{\link{itemfit}}
# #' @describeIn itemfit create bivariate heatmap plots
#' @examples
#' \dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#'
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' ift <- itemfit(fit)
#' # plot the adjusted p values for log odds or transformed correlation
#' plot(ift)
#' # plot unadjusted p values for log odds
#' plot(ift,adjusted = FALSE, type = "logOR")
#'}
#' @export
plot.itemfit <- function(x,type="all",adjusted=TRUE,...){
  item.pair.1 <- item.pair.2 <- unadj.pvalue <- test.adj.pvalue <- NULL
  if(type=="all"||toupper(type)=="LOGOR"){
    if(adjusted==FALSE){
      print(ggplot2::ggplot(extract.itemfit(x,"logOR"),
                            aes(x=factor(item.pair.2),
                                y=factor(item.pair.1),
                                fill=unadj.pvalue))+
              geom_tile()+ scale_fill_gradient(low="red",
                                               high="gray",
                                               limits=c(0,0.05))+
              theme_bw() +
              labs(x = "Item", y = "Item",
                   title = "Heatmap plot for unadjusted p-values of log odds ratio"))
    }else{
      print(ggplot2::ggplot(extract.itemfit(x,"logOR"),
                            aes(x=factor(item.pair.2),
                                y=factor(item.pair.1),
                                fill=test.adj.pvalue))+
              geom_tile()+ scale_fill_gradient(low="red",
                                               high="gray",
                                               limits=c(0,0.05))+
              theme_bw() +
              labs(x = "Item", y = "Item",
                   title = "Heatmap plot for adjusted p-values of log odds ratio"))

    }

  }
  if(type=="all"||toupper(type)=="R"){
      if(adjusted){
        print(ggplot2::ggplot(extract.itemfit(x,"r"),
                              aes(x=factor(item.pair.2),
                                  y=factor(item.pair.1),
                                  fill=test.adj.pvalue))+
                geom_tile()+ scale_fill_gradient(low="red",
                                                 high="gray",
                                                 limits=c(0,0.05))+
                theme_bw() +
                labs(x = "Items", y = "Items",
                     title = "Heatmap plot for adjusted p-values of transformed correlation"))

      }else{

        print(ggplot2::ggplot(extract.itemfit(x,"r"),
                              aes(x=factor(item.pair.2),
                                  y=factor(item.pair.1),
                                  fill=unadj.pvalue))+
                geom_tile()+ scale_fill_gradient(low="red",
                                                 high="gray",
                                                 limits=c(0,0.05))+
                theme_bw() +
                labs(x = "Item", y = "Item",
                     title = "Heatmap plot for unadjusted p-values of transformed correlation"))


      }
    }



}


#' Mesa plot for Q-matrix validation
#'
#' The mesa plot was first proposed by de la Torre and Ma (2016) for graphically illustrating the best q-vector(s) for each item.
#' The q-vector on the edge of the mesa is likely to be the best q-vector.
#'
#' @param x model object of class \code{Qvalidation}
#' @param item a vector specifying which item(s) the plots are drawn for
#' @param type types of the plot. It can be \code{"best"} or \code{"all"}. If \code{"best"},
#'     for all q-vectors requiring the same number of attributes, only the one with the largest PVAF
#'     is plotted, which means \eqn{K_j} q-vectors are plotted; If \code{"all"}, all q-vectors
#'     will be plotted.
#' @param eps the cutoff for PVAF. If not \code{NULL}, it must be a value between 0 and 1. A horizontal line will be drawn accordingly.
#' @param no.qvector the number of q vectors that need to be plotted when \code{type="all"}. The default is 10,
#'        which means the 10 q vectors with the largest PVAFs are plotted.
#' @param data.label logical; To show data label or not?
#' @param original.q.label logical; print the label showing the original q-vector or not?
#' @param auto.ylim logical; create y range automatically or not?
#' @param ... additional arguments passed to \code{plot} function
#' @seealso \code{\link{Qval}}, \code{\link{autoGDINA}}
#' @examples
#'\dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#' Q[1,] <- c(0,1,0)
#' mod1 <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' out <- Qval(mod1,eps = 0.9)
#' item <- c(1,2,10)
#' plot(out,item=item,data.label=FALSE,type="all")
#' plot(out,item=10,type="best",eps=0.95)
#' plot(out,item=10,type="all",no.qvector=6)
#'}
#'
#' @references
#'
#' de la Torre, J., & Ma, W. (2016, August). Cognitive diagnosis modeling: A general framework approach and its implementation in R. A Short Course at the Fourth Conference on Statistical Methods in Psychometrics, Columbia University, New York.
#'

#' @export
plot.Qval <-
  function(x, item, type = "best", no.qvector = 10,
           data.label = TRUE,eps = "auto",
           original.q.label = FALSE,auto.ylim = TRUE,...)
  {
    type <- match.arg(tolower(type), c("best", "all"))
    if (identical(eps, "auto")) eps <- round(x$eps, 2)
    if (x$sequential) {
      Q <- extract.Qval(x, "Q")[, -c(1:2), drop = FALSE]
    } else {
      Q <- extract.Qval(x, "Q")
    }

    K <- ncol(Q)
    L <- 2^K - 1
    patt <- attributepattern(K)
    fullPVAF <- extract.Qval(x, "PVAF")
    q.labels <- rownames(fullPVAF)
    if (is.null(q.labels)) {
      q.labels <- apply(patt[-1, , drop = FALSE], 1, paste0, collapse = "")
    }
    best.q.labels <- c(
      paste0(patt[1, ], collapse = ""),
      q.labels
    )

    for (j in item) {
      original.loc <- which(apply(patt[-1, , drop = FALSE], 1,
                                  function(pattern) all(pattern == Q[j, ])))

      if (type == "all") {
        n.plot <- min(no.qvector, L)
        locations <- order(fullPVAF[, j], decreasing = FALSE)
        locations <- locations[(L - n.plot + 1):L]
        plot.dat <- data.frame(
          rank = seq_len(n.plot),
          q.vector = factor(q.labels[locations], levels = q.labels[locations]),
          PVAF = fullPVAF[locations, j],
          original = locations == original.loc
        )
      } else {
        pvaf.with.zero <- c(0, fullPVAF[, j])
        n.attributes <- rowSums(patt)
        locations.by.size <- split(seq_along(pvaf.with.zero), n.attributes)
        best.locations <- vapply(
          locations.by.size,
          function(locations) locations[which.max(pvaf.with.zero[locations])],
          integer(1)
        )
        best.values <- pvaf.with.zero[best.locations]
        plot.dat <- data.frame(
          rank = seq_along(best.values),
          q.vector = factor(best.q.labels[best.locations],
                            levels = best.q.labels[best.locations]),
          PVAF = best.values,
          original = best.locations == original.loc + 1L
        )
      }

      if (auto.ylim) {
        lower <- max(0, round(min(plot.dat$PVAF) - 0.1, 1))
      } else {
        lower <- 0
      }

      g <- ggplot2::ggplot(plot.dat, ggplot2::aes(x = ggplot2::.data$rank, y = ggplot2::.data$PVAF)) +
        ggplot2::geom_line(color = "#3B536D", linewidth = 0.7) +
        ggplot2::geom_point(color = "#3B536D", size = 2.4) +
        ggplot2::scale_x_continuous(breaks = plot.dat$rank,
                                    labels = plot.dat$q.vector) +
        ggplot2::scale_y_continuous(limits = c(lower, 1),
                                    expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
        ggplot2::labs(
          x = "Candidate q-vector",
          y = "PVAF",
          title = paste("Mesa plot for item", j),
          subtitle = if (type == "best") "Best candidate by number of attributes" else
            paste("Top", n.plot, "candidate q-vectors")
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          plot.title = ggplot2::element_text(face = "bold"),
          plot.subtitle = ggplot2::element_text(color = "#687078")
        )

      if (any(plot.dat$original)) {
        g <- g + ggplot2::geom_point(
          data = plot.dat[plot.dat$original, , drop = FALSE],
          color = "#C44E52", size = 3.2
        )
      }
      if (data.label) {
        g <- g + ggplot2::geom_text(ggplot2::aes(label = round(ggplot2::.data$PVAF, 3)),
                                    vjust = -0.8, size = 3, color = "#3B536D")
      }
      if (!is.null(eps) && eps > 0 && eps < 1) {
        g <- g + ggplot2::geom_hline(yintercept = eps, linetype = "dashed",
                                     color = "#C44E52") +
          ggplot2::annotate("text", x = 1, y = eps, label = paste("eps =", eps),
                            hjust = 0, vjust = -0.5, color = "#C44E52", size = 3)
      }
      if (original.q.label) {
        g <- g + ggplot2::labs(caption = paste(
          "Original q-vector:", q.labels[original.loc],
          "| PVAF:", round(fullPVAF[original.loc, j], 3)
        ))
      }
      print(g)
    }
    invisible(x)
  }


