#' Shared utility functions for estimation
#' 
#' This file contains helper functions shared between SingleGroup_Estimation.R
#' and MultipleGroup_Estimation.R to reduce code duplication.


#' Initialize control parameters for EM estimation
#' 
#' Sets up default control values and validates user-provided overrides.
#' Handles both single-group and multiple-group estimation contexts.
#' 
#' @param control User-provided control list
#' @param ncat Number of categories (nonzero categories in sequential models)
#' @param model Numeric model vector
#' @param is_mg Logical: is this multiple-group estimation?
#' @return Modified control list with defaults applied
#' @keywords internal
init_control <- function(control, ncat, model, is_mg = FALSE) {
  myControl <- list(
    maxitr = 2000,
    conv.crit = 1e-4,
    conv.type = c("ip","mp"),
    nstarts = 3L,
    lower.p = 1e-4,
    upper.p = 1 - 1e-4,
    lower.prior = .Machine$double.eps,
    randomseed = 123456,
    smallNcorrection = c(.0005, .001),
    MstepMessage = FALSE,
    Cpp = FALSE,
    countitemparm = 0
  )
  
  control <- utils::modifyList(myControl, control)
  
  # Expand lower.p and upper.p if needed
  if (length(control$lower.p) == 1) {
    control$lower.p <- rep(control$lower.p, ncat)
  } else if (length(control$lower.p) != ncat) {
    stop("lower.p must have length of 1 or number of nonzero categories", call. = FALSE)
  }
  
  if (length(control$upper.p) == 1) {
    control$upper.p <- rep(control$upper.p, ncat)
  } else if (length(control$upper.p) != ncat) {
    stop("upper.p must have length of 1 or number of nonzero categories", call. = FALSE)
  }
  
  # Handle maxitr expansion
  if (length(control$maxitr) == 1L) {
    control$vmaxitr <- rep(control$maxitr, ncat)
  } else if (length(control$maxitr) != ncat) {
    stop("maxitr must have length of 1 or number of nonzero categories", call. = FALSE)
  } else {
    control$vmaxitr <- control$maxitr
    control$maxitr <- max(control$maxitr)
  }
  
  # Multiple-group specific: no Cpp fast implementation
  if (is_mg && control$Cpp) {
    warning("No fast implementation is available for multiple-group analysis.", call. = FALSE)
    control$Cpp <- FALSE
  }
  
  # Models 7+ (SISM, BUGDINO) use only 1 starting value
  control$nstarts <- ifelse(any(model > 6), 1L, 3L)
  
  set.seed(control$randomseed)
  
  control
}


#' Initialize solver arguments
#' 
#' Sets up default solver parameters for nloptr-based optimization.
#' 
#' @param auglag_args User-provided auglag arguments
#' @param solnp_args User-provided solnp arguments  
#' @param nloptr_args User-provided nloptr arguments
#' @return List containing initialized solver arguments
#' @keywords internal
init_solver_args <- function(auglag_args, solnp_args, nloptr_args) {
  myAuglag_args <- list(
    control.outer = list(
      trace = FALSE,
      method = "nlminb",
      kkt2.check = FALSE,
      eps = 1e-6
    )
  )
  auglag_args <- modifyList(myAuglag_args, auglag_args)
  
  mySolnp_args <- list(trace = 0)
  solnp_args <- modifyList(mySolnp_args, solnp_args)
  
  Mynloptr_args <- list(xtol_rel = 1e-4)
  nloptr_args <- modifyList(Mynloptr_args, nloptr_args)
  
  list(auglag_args = auglag_args, solnp_args = solnp_args, nloptr_args = nloptr_args)
}


#' Initialize solver vector
#' 
#' Expands solver specification to per-item vector.
#' 
#' @param solver Character string or vector specifying solver
#' @param ncat Number of categories
#' @return Character vector of length ncat
#' @keywords internal
init_solver <- function(solver, ncat) {
  if (is.null(solver)) {
    solver <- rep("auto", ncat)
  } else if (length(solver) == 1) {
    solver <- rep(solver, ncat)
  } else if (length(solver) != ncat) {
    stop("solver must have length of 1 or number of items.", call. = FALSE)
  }
  solver
}


#' Check EM convergence
#' 
#' Calculates maximum change across convergence criteria and checks 
#' whether convergence has been achieved.
#' 
#' @param dif_parm List with convergence differences: ip, prior, neg2LL, delt
#' @param conv_type Character vector of convergence types to check
#' @param conv_crit Numeric convergence criterion
#' @param neg2LL_current Current deviance value
#' @return List with maxchg and converged flag
#' @keywords internal
check_em_convergence <- function(dif_parm, conv_type, conv_crit, neg2LL_current) {
  maxchg <- 0
  
  if (any(tolower(conv_type) == "ip")) {
    maxchg <- max(maxchg, dif_parm$ip)
  }
  if (any(tolower(conv_type) == "delta")) {
    maxchg <- max(maxchg, dif_parm$delt)
  }
  if (any(tolower(conv_type) == "mp")) {
    maxchg <- max(maxchg, dif_parm$prior)
  }
  if (any(tolower(conv_type) == "neg2ll")) {
    maxchg <- max(maxchg, abs(dif_parm$neg2LL))
  }
  if (any(tolower(conv_type) == "relneg2ll")) {
    maxchg <- max(maxchg, abs(dif_parm$neg2LL) / neg2LL_current)
  }
  
  list(maxchg = maxchg, converged = maxchg < conv_crit)
}


#' Print EM iteration progress
#' 
#' Prints iteration information with appropriate formatting.
#' 
#' @param itr Iteration number
#' @param maxchg Maximum change
#' @param neg2LL Current deviance
#' @param verbose Verbosity level (1 = inline, 2 = newline)
#' @keywords internal
print_em_progress <- function(itr, maxchg, neg2LL, verbose) {
  if (verbose == 1L) {
    cat('\rIter =', itr, ' Max. abs. change =', formatC(maxchg, digits = 5, format = "f"),
        ' Deviance  =', formatC(neg2LL, digits = 2, format = "f"),
        '                                                                                 ')
  } else if (verbose == 2L) {
    cat('Iter =', itr, ' Max. abs. change =', formatC(maxchg, digits = 5, format = "f"),
        ' Deviance  =', formatC(neg2LL, digits = 2, format = "f"),
        '                                                                                \n')
  }
}


#' Initialize constraint matrices
#' 
#' Creates constraint matrices for monotonicity constraints across items.
#' 
#' @param ncat Number of categories
#' @param mono.constraint Vector of logicals indicating which items have monotonicity constraints
#' @param Kj Vector of attribute counts per item
#' @param reduced.LG List of reduced latent group matrices
#' @param ConstrPairs Optional pre-specified constraint pairs
#' @return List with ConstrType vector and ConstrMatrix list
#' @keywords internal
init_constraint_matrices <- function(ncat, mono.constraint, Kj, reduced.LG, ConstrPairs = NULL) {
  ConstrType <- rep(1, ncat)
  ConstrMatrix <- vector("list", ncat)
  
  if (is.null(ConstrPairs)) {
    ConstrPairs <- vector("list", ncat)
    for (j in seq_len(ncat)) {
      if (mono.constraint[[j]]) {
        ConstrType[j] <- 3
        ConstrPairs[[j]] <- partial_order2(Kj[j], reduced.LG[[j]])
        nctj <- nrow(ConstrPairs[[j]])
        tmp <- matrix(0, nctj, nrow(reduced.LG[[j]]))
        tmp[matrix(c(seq_len(nctj), ConstrPairs[[j]][, 1]), ncol = 2)] <- 1
        tmp[matrix(c(seq_len(nctj), ConstrPairs[[j]][, 2]), ncol = 2)] <- -1
        ConstrMatrix[[j]] <- tmp
      }
    }
  }
  
  list(ConstrType = ConstrType, ConstrMatrix = ConstrMatrix, ConstrPairs = ConstrPairs)
}


#' Initialize design matrices
#' 
#' Creates design matrices for all items based on model type and Q-matrix.
#' 
#' @param ncat Number of categories
#' @param model Vector of model codes
#' @param rule Vector of condensation rules
#' @param Kj Vector of attribute counts per item
#' @param reduced.LG List of reduced latent group matrices
#' @param Q Q-matrix (attributes only, no item/category columns)
#' @param originalQ Original Q-matrix with item/category columns for MS-DINA
#' @param no.bugs Number of bugs for BUGDINO model
#' @return List of design matrices
#' @keywords internal
init_design_matrices <- function(ncat, model, rule, Kj, reduced.LG, Q, originalQ, no.bugs = NULL) {
  if (any(model == -1)) {
    stop("design.matrix must be provided for user-defined models.", call. = FALSE)
  }
  
  DesignMatrices <- vector("list", ncat)
  for (j in seq_len(ncat)) {
    if (model[j] == 6) {  # MS-DINA
      DesignMatrices[[j]] <- designmatrix(
        model = model[j],
        Qj = originalQ[which(originalQ[, 1] == j), -c(1:2), drop = FALSE]
      )
    } else if (model[j] %in% c(7, 8)) {  # BUGDINO, SISM
      DesignMatrices[[j]] <- designmatrix(model = model[j], Qj = Q[j, ], no.bugs = no.bugs)
    } else if (rule[j] >= 0 && rule[j] <= 3) {
      DesignMatrices[[j]] <- designM(Kj[j], rule[j], reduced.LG[[j]])
    }
  }
  
  DesignMatrices
}
