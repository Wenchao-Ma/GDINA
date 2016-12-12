#' Generate hierarchical attribute structures
#'
#' This function can be used to generate hierarchical attributes structures, and
#' to provide prior joint attribute distribution with hierarchical structures.
#'
#' @param hierarchy.list a list specifying the hierarchical structure between attributes. Each
#'     element in this list specifies a DIRECT prerequisite relation between two or more attributes.
#'     See \code{example} for more information.
#' @param K the number of attributes involved in the assessment
#' @param att.prob How are the probabilities for latent classes simulated? It can be \code{"random"} or \code{"uniform"}.
#'
#' @return att.str reduced latent classes under the specified hierarchical structure
#' @return impossible.latentclass impossible latent classes under the specified hierarchical structure
#' @return att.prob probabilities for all latent classes; 0 for impossible latent classes
#'
#' @seealso \code{\link{GDINA}}, \code{\link{autoGDINA}}
#' @author {Wenchao Ma, Rutgers University, \email{wenchao.ma@@rutgers.edu} \cr Jimmy de la Torre, The University of Hong Kong}
#' @export
#'
#' @examples
#' \dontrun{
#' #################
#' #
#' # Leighton et al. (2004, p.210)
#' #
#' ##################
#' # linear structure A1->A2->A3->A4->A5->A6
#' K <- 6
#' linear=list(c(1,2),c(2,3),c(3,4),c(4,5),c(5,6))
#' att.structure(linear,K)
#'
#' # convergent structure A1->A2->A3->A5->A6;A1->A2->A4->A5->A6
#' K <- 6
#' converg <- list(c(1,2),c(2,3),c(2,4),
#'                c(3,4,5), #this is how to show that either A3 or A4 is a prerequisite to A5
#'                c(5,6))
#'att.structure(converg,K)
#'
#' # convergent structure [the difference between this one and the previous one is that
#' #                       A3 and A4 are both needed in order to master A5]
#' K <- 6
#' converg2 <- list(c(1,2),c(2,3),c(2,4),
#'                c(3,5), #this is how to specify that both A3 and A4 are needed for A5
#'                c(4,5), #this is how to specify that both A3 and A4 are needed for A5
#'                c(5,6))
#'att.structure(converg2,K)
#'
#' # divergent structure A1->A2->A3;A1->A4->A5;A1->A4->A6
#' diverg <- list(c(1,2),
#'                c(2,3),
#'                c(1,4),
#'                c(4,5),
#'                c(4,6))
#'att.structure(diverg,K)
#'
#' # unstructured A1->A2;A1->A3;A1->A4;A1->A5;A1->A6
#' unstru <- list(c(1,2),c(1,3),c(1,4),c(1,5),c(1,6))
#' att.structure(unstru,K)
#'
#' ## See Example 4 and 5 in GDINA function
#'}
att.structure <- function(hierarchy.list=NULL,K,att.prob="uniform"){
  patt <- alpha(K)
  if (!is.null(hierarchy.list)){

        impos.id <- lapply(hierarchy.list,function(x){
          k <- length(x)
          if(max(x)>K) stop("Maximum element of hierarchy list cannot be greater than K.",call. = FALSE)
          which(rowSums(patt[,x[1:(k-1)],drop=FALSE])<patt[,x[k]])

    })


    if (is.list(impos.id)) {
      impos.id <- unique(c(unlist(impos.id)))
    }else{
      impos.id <- unique(c(impos.id))
    }
    pos.id <- setdiff(seq(1,nrow(patt)),impos.id)
    red.patt <- patt[-impos.id,]
  }else{
    red.patt <- patt
    pos.id <- seq(1,nrow(patt))
  }
  prior <- rep(0,nrow(patt))
  if(tolower(att.prob)=="uniform"){
    prior[pos.id] <- 1/length(pos.id)
    }else if(tolower(att.prob)=="random"){
    prior[pos.id] <- runif(length(pos.id))
    prior <- prior/sum(prior)
    }
    names(prior) <- apply(patt,1,function(x){paste(x,collapse = "")})

  return(list(att.str=red.patt,impossible.latentclass=sort(impos.id),att.prob=prior))
}
