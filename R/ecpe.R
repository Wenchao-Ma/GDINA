#' Examination for the Certificate of Proficiency in English (ECPE) data
#'
#' Examination for the Certificate of Proficiency in English (ECPE) data (the grammar section) has been used
#' in Henson and Templin (2007), Templin and Hoffman (2013), Feng, Habing, and Huebner (2013), and
#' Templin and Bradshaw (2014), among others.
#'
#' The data consists of responses of 2922 examinees to 28 items involving 3 attributes.
#' Attribute 1 is morphosyntactic rules, Attribute 2 is cohesive rules and
#' Attribute 3 is lexical rules.
#'
#'
#' @format A list of responses and Q-matrix with components:
#' \describe{
#' \item{\code{dat}}{Responses of 2922 examinees to 28 items.}
#' \item{\code{Q}}{The \eqn{28 \times 3} Q-matrix.}
#' }
#'
#' @examples
#' \dontrun{
#' mod1 <- GDINA(ecpe$dat,ecpe$Q)
#' mod1
#' summary(mod1)
#'
#' mod2 <- GDINA(ecpe$dat,ecpe$Q,model="RRUM")
#' mod2
#' anova(mod1,mod2)
#' # You may compare the following results with Feng, Habing, and Huebner (2013)
#' itemparm(mod2,"rrum")
#'}
#'
#' @references
#'
#' Feng, Y., Habing, B. T., & Huebner, A. (2013). Parameter estimation of the reduced RUM using the EM algorithm. \emph{Applied Psychological Measurement}, 0146621613502704.
#'
#' Henson, R. A., & Templin, J. (2007, April). Large-scale language assessment using cognitive diagnosis models. Paper presented at the annual meeting of the National Council for Measurement in Education in Chicago, Illinois.
#'
#' Templin, J., & Bradshaw, L. (2014). Hierarchical diagnostic classification models: A family of models for estimating and testing attribute hierarchies. \emph{Psychometrika, 79}, 317-339.
#'
#' Templin, J., & Hoffman, L. (2013). Obtaining diagnostic classification model estimates using Mplus. \emph{Educational Measurement: Issues and Practice, 32}, 37-50.
"ecpe"
