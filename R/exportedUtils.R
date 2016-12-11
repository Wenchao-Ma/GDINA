#' Create a block diagonal matrix
#'
#' @param mlist a list of matrices
#' @param fill value to fill the non-diagnoal elements
#'
#' @return a block diagonal matrix
#' @export
#' @seealso \code{bdiag} in \pkg{Matrix}
#' @examples
#'
#' m1 <- bdiagMatrix(list(matrix(1:4, 2), diag(3)))
#' m2 <- bdiagMatrix(list(matrix(1:4, 2), diag(3)),fill = NA)
#'
bdiagMatrix <- function(mlist,fill=0){

  loc <- sapply(mlist,dim)
  out <- matrix(fill,rowSums(loc)[1],rowSums(loc)[2])
  cr <- cc <- 1
  for(i in 1:length(mlist)){
    out[cr:(cr+nrow(mlist[[i]])-1),cc:(cc+ncol(mlist[[i]])-1)] <- mlist[[i]]
    cr <- cr + nrow(mlist[[i]])
    cc <- cc + ncol(mlist[[i]])
  }
  out
}


#' Combine R Objects by Columns
#'
#' Combine a sequence of vector, matrix or data-frame arguments by columns. Vector is treated as a column matrix.
#'
#' @param ... vectors or matrices
#' @param fill a scalar used when these objects have different number of rows.
#'
#' @return a data frame
#' @export
#' @seealso \code{\link[base]{cbind}}
#' @examples
#' cjoint(2,c(1,2,3,4),matrix(1:6,2,3))
#' cjoint(v1 = 2, v2 = c(3,2), v3 = matrix(1:6,3,2),
#'        v4 = data.frame(c(3,4,5,6,7),rep("x",5)),fill = 99)
#'
cjoint <- function(...,fill = NA){
  arg <- list(...)
  if(length(arg) == 0) {
    return(NULL)
  }else{
    maxrow <- max(sapply(arg,function(x){nrow(as.matrix(x))}))
    r.arg <- lapply(arg,function(x){
      df <- data.frame(x)
      nrow.x <- nrow(df)
      if(nrow.x < maxrow) {
        rbind(df,matrix(fill,maxrow-nrow.x,ncol(df),dimnames = list(NULL,names(df))))
        }else{
          data.frame(x)
        }})
    out <- do.call(cbind,r.arg)
    names(out) <- NULL
    return(out)
  }

}

