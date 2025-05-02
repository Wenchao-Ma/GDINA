#' Mental health symptom profiles data (polytomous data)
#'
#' Mental health data used in Tan, et al. (2023).
#'
#' The data consists of polytomous responses of 719 college students (34.6\% men, 83.9\% White, and 16.3\% first-year students)
#' to 40 items that measure four attributes. The item categories were collapsed to make sure each category has at least 50 responses if possible.
#'
#'
#' @format A list of dichotmous/polytomous responses and Qc-matrix with components:
#' \describe{
#' \item{\code{dat}}{Responses of 719 students to 40 items.}
#' \item{\code{Qc}}{The Qc-matrix.}
#' }
#'
#' @seealso \code{tan2023} for binary version of the data and \code{tan2023p25} for polytomous data with each cateogry having 25 or more responses, when possible.
#'
#' @references
#'
#' Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive Diagnosis Modeling. \emph{Journal of Statistical Software, 93(14)}, 1-26.
#'
#' @examples
#' \dontrun{
#' mod1 <- GDINA(tan2023p50$dat,tan2023p50$Qc,sequential = TRUE)
#' mod1
#' summary(mod1)

#'}
#'
#' @references
#'
#' Tan, Z., de la Torre, J., Ma, W., Huh, D., Larimer, M. E., & Mun, E. Y. (2023). A Tutorial on Cognitive Diagnosis Modeling for Characterizing Mental Health Symptom Profiles Using Existing Item Responses. \emph{Prevention science : the official journal of the Society for Prevention Research, 24}, 480–492. https://doi.org/10.1007/s11121-022-01346-8
"tan2023p50"
