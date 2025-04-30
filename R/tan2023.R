#' Mental health symptom profiles data
#'
#' Mental health data used in Tan, et al. (2023).
#'
#' The data consists of dichotomized responses of 719 college students (34.6% men, 83.9% White, and 16.3% first-year students)
#' to 40 items that measure four attributes.
#'
#'
#' @format A list of binary responses and Q-matrix with components:
#' \describe{
#' \item{\code{dat}}{Responses of 719 students to 40 items.}
#' \item{\code{Q}}{The \eqn{40 \times 4} Q-matrix.}
#' }
#'
#' @references
#'
#' Ma, W., & de la Torre, J. (2020). GDINA: An R Package for Cognitive Diagnosis Modeling. \emph{Journal of Statistical Software, 93(14)}, 1-26.
#'
#' @examples
#' \dontrun{
#' mod1 <- GDINA(tan2023$dat,tan2023$Q)
#' mod1
#' summary(mod1)

#'}
#'
#' @references
#'
#' Tan, Z., de la Torre, J., Ma, W., Huh, D., Larimer, M. E., & Mun, E. Y. (2023). A Tutorial on Cognitive Diagnosis Modeling for Characterizing Mental Health Symptom Profiles Using Existing Item Responses. \emph{Prevention science : the official journal of the Society for Prevention Research, 24}, 480–492. https://doi.org/10.1007/s11121-022-01346-8
"tan2023"
