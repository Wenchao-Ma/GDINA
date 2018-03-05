#' Generate design matrix
#'
#' This function generates the design matrix \eqn{M_j} in de la Torre (2011)
#'
#' @param Kj the number of attributes for item j
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
#' @export
#'
#'
#'
designmatrix <- function(Kj, model = "GDINA", Qj = NULL) {
  stopifnot(is.nonNegativeInteger(Kj))
  if (is.character(model)) {
    stopifnot(toupper(model) %in% c("GDINA", "DINA", "DINO", "ACDM", "LLM", "RRUM", "MSDINA"))
    m <-
      which(c("GDINA", "DINA", "DINO", "ACDM", "LLM", "RRUM", "MSDINA") == toupper(model))
  } else if (is.numeric(model)) {
    if (!is.nonNegativeInteger(model) | model > 7) {
      stop('model must be "GDINA", "DINA","DINO","ACDM","LLM", "RRUM" or "MSDINA".',
           call. = FALSE)
    } else{
      m <- model + 1
    }
  } else{
    stop("model is not correctly specified.", call. = FALSE)
  }
  if (m == 7) {
    # MSDINA
    # Kj is not necessary
    if (is.null(Qj) ||
        nrow(Qj) < 2 ||
        max(Qj) > 1)
      stop("Qj is not correctly specified for the MS-DINA model.", call. = F)
    Qj <- as.matrix(Qj)
    if (any(colSums(Qj) == 0))
      Qj <- Qj[, -which(colSums(Qj) == 0)]
    Ks <- rowSums(Qj)
    Kj <- sum(apply(Qj, 2, max))
    patt <- attributepattern(Kj)

    D <-
      matrix(c(rep(1, nrow(patt)), colSums(Qj %*% t(patt) == Ks) > 0), ncol = 2)

  } else{
    D <- designM(Kj, m - 1)
  }
  return(D)

}
