#' Graphical user interface of the GDINA function
#'
#' An experimental interactive Shiny application for running GDINA function
#'
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
  shiny::runApp(appDir = system.file("shiny", package="GDINA"))

}
