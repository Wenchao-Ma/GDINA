#' Graphical user interface of the GDINA function
#'
#' An interactive Shiny application for running GDINA function. See Ma and de la Torre (2019) and de la Torre and Akbay (2019) for tutorials.
#'
#' @author Wenchao Ma, The University of Minnesota, \email{wma@umn.edu}
#' @references
#'
#' de la Torre, J & Akbay, L. (2019). Implementation of Cognitive Diagnosis Modeling using the GDINA R Package. \emph{Eurasian Journal of Educational Research, 80}, 171-192.
#'
#' Ma, W., & de la Torre, J. (2019). Digital Module 05: Diagnostic measurement-The G-DINA framework.\emph{ Educational Measurement: Issues and Practice, 39}, 114-115.
#'
#' Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive Diagnosis Modeling. \emph{Journal of Statistical Software, 93(14)}, 1-26.
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
  cat("Please wait while loading...\n")
  shiny::runApp(appDir = system.file("shiny", package="GDINA"),
                launch.browser = TRUE)

}
