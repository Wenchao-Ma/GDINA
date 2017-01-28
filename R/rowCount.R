#' Count and label unique rows in data frame
#'
#' @param df a data frame or matrix
#'
#' @return freq the number of rows
#' @return group the data frame with a column named row.no giving unique labels for all unique rows
#' @export
#' @import data.table
#' @examples
#'
#' df <- data.frame(V1=c(1L,2L),V2=LETTERS[1:3],V3=rep(1,12))
#' rowCount(df)
#'
#'
#'
rowCount <- function(df){
  DT <- data.table(df)
  varb <- colnames(DT)
  freq <- DT[,.N,by=c(varb)]
  DT$gr <- as.numeric(factor(apply(DT,1,paste0,collapse="")))
  return(list(freq=data.frame(freq),group=data.frame(DT)))
}

#' Count the frequency of a row vector in a data frame
#'
#' @param df a data frame or matrix
#' @param vec the vector for matching
#' @return count the number of vector vec in the data frame
#' @return row.no row numbers of the vector vec in the data frame
#' @export
#'
#' @examples
#'
#' df <- data.frame(V1=c(1L,2L),V2=LETTERS[1:3],V3=rep(1,12))
#' rowMatch(df,c(2,"B",1))
#'
#'
#'
rowMatch <- function(df,vec=NULL){
  logicalvec <- apply(df,1,paste0,collapse="")==paste0(vec,collapse = "")
  return(list(count=sum(logicalvec),row.no=which(logicalvec)))
}

#' Unique values in a vector
#'
#' @param vec a vector
#' @return sorted unique values
#' @export
#'
#' @examples
#'
#' vec <- c(4,2,3,5,4,4,4)
#' unique_only(vec)
#' # see the difference from unique
#' unique(vec)
#'
#' vec <- letters[1:5]
#' unique_only(vec)
#'
#'
#'
#' @seealso \link[base]{unique}
#'
unique_only <- function(vec){
  vec <- sort(vec)
  unique(vec)[which(table(vec)==1)]
}
