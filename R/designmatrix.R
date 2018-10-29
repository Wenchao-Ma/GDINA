#' Generate design matrix
#'
#' This function generates the design matrix for an item
#'
#' @param Kj Required except for the MS-DINA model; The number of attributes required for item j
#' @param model the model associated with the design matrix; It can be "GDINA","DINA","DINO", "ACDM" or "MSDINA".
#'        The default is "GDINA". Note that models "LLM" and "RRUM" have the same design matrix as the ACDM.
#' @param Qj the Q-matrix for item j; This is required for "MSDINA" model; The number of rows is equal to the number
#'        of strategies and the number of columns is equal to the number of attributes.
#'
#' @return a design matrix (Mj). See de la Torre (2011) for details.
#' @references
#'
#' de la Torre, J. (2011). The generalized DINA model framework. \emph{Psychometrika, 76}, 179-199.
#'
#'
#' @examples
#' \dontrun{
#' designmatrix(Kj = 2, model = "GDINA")
#' designmatrix(Kj = 3, model = "DINA")
#' msQj <- matrix(c(1,0,0,1,
#'                  1,1,0,0),nrow=2,byrow=TRUE)
#' designmatrix(model = "MSDINA",Qj = msQj)
#' }
#'
#' @export
#'
#'
#'
designmatrix <- function(Kj = NULL, model = "GDINA", Qj = NULL) {

  m <- model.transform(model,J = 1)
  if(m < 0) # log GDINA, logit GDINA => GDINA
    m <- 0
  if (m == 6) {
    # MSDINA
    # Kj is not necessary
    if (is.null(Qj) || max(Qj) > 1){
      stop("Qj is not correctly specified for the MS-DINA model.", call. = F)
    }else if(nrow(Qj)==1){
      m <- 1
      Kj <- sum(Qj)
    }else{
      Qj <- as.matrix(Qj)
      if (any(colSums(Qj) == 0))
        Qj <- Qj[, -which(colSums(Qj) == 0)]
      Ks <- rowSums(Qj)
      Kj <- sum(apply(Qj, 2, max))
      patt <- attributepattern(Kj)

      D <- matrix(c(rep(1, nrow(patt)), colSums(Qj %*% t(patt) == Ks) > 0), ncol = 2)
    }

  }

  if(m!=6){
    D <- designM(Kj, m)
  }
  return(D)

}
