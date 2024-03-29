% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExportedFuncs.R
\name{designmatrix}
\alias{designmatrix}
\title{Generate design matrix}
\usage{
designmatrix(Kj = NULL, model = "GDINA", Qj = NULL, no.bugs = 0)
}
\arguments{
\item{Kj}{Required except for the MS-DINA model; The number of attributes required for item j}

\item{model}{the model associated with the design matrix; It can be "GDINA","DINA","DINO", "ACDM","LLM", "RRUM", "MSDINA", "BUGDINO", and "SISM".
The default is "GDINA". Note that models "LLM" and "RRUM" have the same design matrix as the "ACDM".}

\item{Qj}{the Q-matrix for item j; This is required for "MSDINA", and "SISM" models; The number of rows is equal to the number
of strategies for "MSDINA", and the number of columns is equal to the number of attributes.}

\item{no.bugs}{the number of bugs (or misconceptions). Note that bugs must be given in the last no.bugs columns.}
}
\value{
a design matrix (Mj). See de la Torre (2011) for details.
}
\description{
This function generates the design matrix for an item
}
\examples{
\dontrun{
designmatrix(Kj = 2, model = "GDINA")
designmatrix(Kj = 3, model = "DINA")
msQj <- matrix(c(1,0,0,1,
                 1,1,0,0),nrow=2,byrow=TRUE)
designmatrix(model = "MSDINA",Qj = msQj)
}

}
\references{
de la Torre, J. (2011). The generalized DINA model framework. \emph{Psychometrika, 76}, 179-199.
}
