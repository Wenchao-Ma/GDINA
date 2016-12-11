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
#' @import shinydashboard
#' @import shiny
#' @export
startGDINA <- function() {
  shiny::runApp(appDir = system.file("shiny", package="GDINA"))

}
