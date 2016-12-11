#' Design matrix for parameter transformation
#'
#' This function calculates the design matrix \eqn{M_j} in de la Torre (2011),
#' which can be useful for parameter transformation between probability of success and
#' delta.
#'
#' @param Kj the number of attributes for item j
#' @param model the model fitted to item j; it can be "GDINA","DINA","DINO","ACDM","LLM",or "RRUM".
#'        The default is "GDINA".
#'
#' @return a design matrix (Mj) which plays a critical role of transforming probability of success with
#'         delta parameters. See de la Torre (2011) for details.
#' @references
#'
#' de la Torre, J. (2011). The generalized DINA model framework. \emph{Psychometrika, 76}, 179-199.
#'
#' @export
#'
#' @examples
#'
#' # transform probability of success to delta
#' # based on saturated GDINA model
#' # assuming an item with 2 attributes has the following
#' # probability of success for 00, 10, 01 and 11
#' pj <- c(0.2,0.4,0.6,0.8)
#' Mj <- designmatrix(2)
#' # delta parameters can be calculated in this way:
#' deltaj <- solve(Mj)%*%pj
#' # for reduced CDMs, OLS or WLS may be used
#'
#'
designmatrix <- function(Kj,model="GDINA"){
  if (!is.positiveInteger(Kj)) stop('Kj must be positive integer.',call. = FALSE)
  if (toupper(model)=="GDINA"|model==0){
    return(designM_GDINA(Kj))
  }else if (toupper(model)%in%c("DINA","DINO","ACDM","LLM","RRUM")|model>=1){
    if(is.character(model)) {
      m <- which(c("DINA","DINO","ACDM","LLM","RRUM")==toupper(model))
    }else if(is.numeric(model)){
      if (!is.positiveInteger(model)|model>5) {
        stop('model must be "DINA","DINO","ACDM","LLM",or "RRUM".',call. = FALSE)
      }else{
        m <- model
      }
    }
    return(designM(alpha(Kj),m))
  }
}
