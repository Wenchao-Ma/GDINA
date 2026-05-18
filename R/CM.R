#' @title Calculate Misclassification Matrices
#'
#' @description
#' This function computes profile-level and attribute-level misclassification
#' matrices from a fitted \code{\link{GDINA}} model. This function is only applicable
#' to models with binary attributes.
#'
#' @param object An estimated GDINA object returned from \code{\link{GDINA}}.
#' @param classification A character string specifying the classification rule.
#' Supported values are \code{"MAP"}, \code{"MLE"}, and \code{"EAP"}.
#' Alternatively, a matrix of user-supplied classifications can be provided,
#' with one row per respondent and one column per attribute.
#' @param matrixtype A character string specifying which matrix to return.
#' Supported values are \code{"profile"}, \code{"attribute"}, and
#' \code{"both"}.
#'
#' @return Depending on \code{matrixtype}, the function returns one of the
#' following:
#' \describe{
#' \item{\code{"profile"}}{A class-by-class matrix of estimated
#'   \eqn{P(\hat{\alpha} = c' \mid \alpha = c)}.}
#' \item{\code{"attribute"}}{A list of \eqn{2 \times 2} matrices, one for
#'   each attribute, with estimated
#'   \eqn{P(\hat{\alpha}_k = a' \mid \alpha_k = a)}.}
#' \item{\code{"both"}}{A list with elements \code{profile_classification}
#'   and \code{att_classification}.}
#' }
#' For both pattern-level and attribute-level matrices, rows correspond to 
#' true classes and columns correspond to classified classes.
#' @details
#' The profile-level matrix is computed from posterior latent class
#' probabilities. The attribute-level matrices are computed from marginal
#' attribute mastery probabilities and the chosen classifications.
#' 
#'
#' @author Wenchao Ma, The University of Minnesota, \email{wma@umn.edu}
#'
#' @examples
#' \dontrun{
#' dat <- realdata_ECPE$dat
#' Q <- realdata_ECPE$Q
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' CM(fit)
#' CM(fit, matrixtype = "both")
#' }
#' @export



CM <- function(object,classification="MAP",matrixtype="profile"){

  matrixtype <- match.arg(matrixtype, c("profile", "attribute", "both"))

  post <- exp(GDINA::indlogPost(object))
  pattern <- extract(object,"attributepattern")

  if (is.character(classification)) {
    classification <- match.arg(classification, c("MAP", "MLE", "EAP"))
  }

  if(length(classification)==1 & (classification=="MAP")){
    gr <- max.col(extract(object,what = "logposterior.i"))
    pp <- pattern[gr,]
  } else if(length(classification)==1 & (classification=="MLE")){
    gr <- max.col(extract(object,what = "loglikelihood.i"))
    pp <- pattern[gr,]
  } else if (length(classification)==1 & classification=="EAP"){
    pp <- personparm(object, what="EAP")
     gr <- whichrow_AinB(pp, pattern)
  } else if (length(classification)!=1){
    stopifnot(is.matrix(classification))
    pp <- classification
     gr <- whichrow_AinB(pp, pattern)
  }
  misclassification <- att_classification <- NULL
  if(matrixtype %in% c("profile", "both")) {
  n_classes <- extract(object,what="nLC")
  assignment <- diag(n_classes)[gr, , drop = FALSE] # N x C
  misclassification <- crossprod(post, assignment) # C x C P(alpha_c hat|alpha_c)
  misclassification <- misclassification / pmax(rowSums(misclassification), .Machine$double.eps) # normalization per row
  colnames(misclassification) <- rownames(misclassification) <- apply(pattern, 1, paste, collapse = "") # class labels
  }
  if(matrixtype %in% c("attribute","both")){
    mp <- personparm(object,what="mp")
    att_classification <- list()
    for(k in seq_len(ncol(mp))){
      mpk <- cbind(1-mp[,k], mp[, k]) #N x 2 matrix of P(alpha_k=0) and P(alpha_k=1)
      att_classification_k <- cbind(1-pp[,k], pp[, k]) #N x 2 matrix of classified alpha_k=0 and alpha_k=1
      att_classification[[k]] <- crossprod(mpk, att_classification_k) # 2 x 2 matrix of P(alpha_k hat|alpha_k)
      att_classification[[k]] <- att_classification[[k]] / pmax(rowSums(att_classification[[k]]), .Machine$double.eps)
      colnames(att_classification[[k]]) <- rownames(att_classification[[k]]) <- c(paste0("alpha", k, "=0"), paste0("alpha", k, "=1"))
    }
    names(att_classification) <- paste0("Attribute ", seq_len(ncol(mp)))
  }


  list(profile_classification = misclassification, att_classification = att_classification)
}
