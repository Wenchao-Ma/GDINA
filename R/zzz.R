#  zzz.R
#  Copy from R package MICE: Multivariate Imputation by Chained Equations
#    Copyright (c) 1999-2010 TNO Quality of Life, Leiden
#
#	 This file is part of the R package MICE.
#
#------------------------------.onAttach-------------------------------
.onAttach <- function(...){
  d <- utils::packageDescription("GDINA")
  base::packageStartupMessage("==============================\n",
                        paste(d$Package,"package - version",d$Version),"\n",
                        paste("          ",d$Date,"\n"),
                        "==============================\n")
  return()
}

