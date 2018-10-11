#' Graphical user interface of the GDINA function
#'
#' An experimental interactive Shiny application for running GDINA function
#'
#' @author {Wenchao Ma, The University of Alabama, \email{wenchao.ma@@ua.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @examples
#' \dontrun{
#' library(shiny)
#' library(shinydashboard)
#' startGDINA()
#' }
#'
#' @export
startGDINA <- function() {
  if (!requireNamespace(c("shiny","shinydashboard"), quietly = TRUE)) {
    stop("shiny and shinydashboard needed for startGDINA. Please install them.",
         call. = FALSE)
  }
  cat("Please wait while loading...\n")
  shiny::runApp(appDir = system.file("shiny", package="GDINA"),
                launch.browser = TRUE)

}
