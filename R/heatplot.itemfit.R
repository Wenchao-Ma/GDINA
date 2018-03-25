#' Item fit plots
#'
#' Create plots of bivariate heatmap for item fit
#'
#' @param object model object of class \code{itemfit}
#' @param ... additional arguments
#' @seealso \code{\link{GDINA}}, \code{\link{itemfit}}
#'
#' @examples
#' \dontrun{
#' dat <- sim10GDINA$simdat
#' Q <- sim10GDINA$simQ
#'
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' ift <- itemfit(fit)
#' heatplot(ift)
#'}
#' @export
#'
heatplot <- function(object,...){
  UseMethod("heatplot")
}
#' @describeIn itemfit plot bivariate heatmap for misfit detection
#' @export
heatplot.itemfit <- function(object,...){
item.pair.1 <- item.pair.2 <- unadj.pvalue <- test.adj.pvalue <- NULL
  print(ggplot2::ggplot(extract.itemfit(object,"logOR"),
                        aes(x=factor(item.pair.2),
                            y=factor(item.pair.1),
                            fill=unadj.pvalue))+
          geom_tile()+ scale_fill_gradient(low="red",
                                           high="gray",
                                           limits=c(0,0.05))+
          theme_bw() +
          labs(x = "Items", y = "Items",
               title = "Heatmap plot for unadjusted p-values of log odds ratio"))

  print(ggplot2::ggplot(extract.itemfit(object,"logOR"),
                  aes(x=factor(item.pair.2),
                      y=factor(item.pair.1),
                      fill=test.adj.pvalue))+
    geom_tile()+ scale_fill_gradient(low="red",
                                     high="gray",
                                     limits=c(0,0.05))+
    theme_bw() +
    labs(x = "Items", y = "Items",
         title = "Heatmap plot for adjusted p-values of log odds ratio"))

  print(ggplot2::ggplot(extract.itemfit(object,"r"),
                  aes(x=factor(item.pair.2),
                      y=factor(item.pair.1),
                      fill=unadj.pvalue))+
    geom_tile()+ scale_fill_gradient(low="red",
                                     high="gray",
                                     limits=c(0,0.05))+
    theme_bw() +
    labs(x = "Items", y = "Items",
         title = "Heatmap plot for unadjusted p-values of transformed correlation"))

  print(ggplot2::ggplot(extract.itemfit(object,"r"),
                  aes(x=factor(item.pair.2),
                      y=factor(item.pair.1),
                      fill=test.adj.pvalue))+
    geom_tile()+ scale_fill_gradient(low="red",
                                     high="gray",
                                     limits=c(0,0.05))+
    theme_bw() +
    labs(x = "Items", y = "Items",
         title = "Heatmap plot for adjusted p-values of transformed correlation"))

}
