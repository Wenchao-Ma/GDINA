#' @title Three-step ML correction for covariate regression
#'
#' @description
#' Fit covariate regression models for estimated attribute classifications or
#' latent profile classifications from a fitted
#' \code{\link{GDINA}} object using the maximum-likelihood three-step
#' correction (Vermunt, 2010). The correction uses the estimated misclassification matrix
#' \eqn{P(W=s\mid X=t)} obtained from \code{\link{CM}} function.
#'
#' @param object An estimated \code{GDINA} object returned from
#'   \code{\link{GDINA}}.
#' @param formula A one-sided or two-sided formula specifying the covariates.
#'   Only the right-hand side is used.
#' @param data A data frame containing the covariates in \code{formula}. It
#'   must contain one row per respondent in \code{object}.
#' @param level A character string specifying whether the corrected regression
#'   is carried out at the attribute level (\code{"attribute"}) or the latent
#'   profile level (\code{"profile"}).
#' @param attribute Optional integer vector giving the attributes to analyze when
#'   \code{level = "attribute"}. If \code{NULL}, all attributes are analyzed.
#' @param classification A character string specifying the classification rule.
#'   Supported values are \code{"MAP"}, \code{"MLE"}, and \code{"EAP"}.
#'   Alternatively, a matrix of user-supplied classifications can be provided,
#'   with one row per respondent and one column per attribute.
#' @param reference A character string specifying the reference class for
#'   profile-level regression. Supported values are \code{"last"}
#'   (default) and \code{"first"}. Ignored when
#'   \code{level = "attribute"}.
#' @param conf.level Confidence level used for Wald intervals.
#'
#' @details
#' With notations from Vermunt (2010), let \eqn{W} denote the estimated class,
#' \eqn{X} the unobserved true class,
#' and \eqn{Z} the covariates. The corrected likelihood implemented here is
#' based on
#' \deqn{\sum_i \log\left\{\sum_t P(W_i=s_i\mid X_i=t)P(X_i=t\mid Z_i)\right\},}
#' which corresponds to the three-step ML correction. The posterior
#' class weights used after fitting are based on
#' \deqn{P(X_i=t\mid W_i=s_i,Z_i) \propto P(W_i=s_i\mid X_i=t)P(X_i=t\mid Z_i),}
#' At the attribute level, a separate binary logistic
#' regression is fitted for each selected attribute. At the profile level, a
#' multinomial logistic regression is fitted using the latent profiles as the
#' outcome categories.
#'
#' This function is currently limited to models with binary attributes.
#'
#' @return
#' A list containing the corrected and naive regression fits. For
#' \code{level = "attribute"}, \code{results} is a named list with one entry
#' per fitted attribute. For \code{level = "profile"}, \code{results} contains
#' the multinomial regression fit. Each fit includes coefficient estimates,
#' Wald summaries, fitted probabilities, posterior class weights, and the
#' misclassification matrix used for correction.
#'
#' @references
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved three-step approaches. \emph{Political analysis, 18}(4), 450-469.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' N <- 1000
#' Z <- data.frame(
#'   z_cont = rnorm(N),
#'   z_bin = rbinom(N, 1, 0.45),
#'   z_cat = factor(sample(c("A", "B", "C"), N, replace = TRUE))
#' )
#'
#' p1 <- plogis(-0.3 + 0.9 * Z$z_cont + 0.8 * Z$z_bin - 0.5 * (Z$z_cat == "B") +
#'   0.4 * (Z$z_cat == "C"))
#' p2 <- plogis(0.2 - 0.8 * Z$z_cont + 0.5 * Z$z_bin + 0.6 * (Z$z_cat == "B") -
#'   0.3 * (Z$z_cat == "C"))
#' alpha <- cbind(rbinom(N, 1, p1), rbinom(N, 1, p2))
#'
#' Q <- matrix(c(
#'   1, 0,
#'   1, 0,
#'   1, 0,
#'   0, 1,
#'   0, 1,
#'   0, 1,
#'   1, 1,
#'   1, 1
#' ), byrow = TRUE, ncol = 2)
#' gs <- data.frame(guess = rep(0.15, nrow(Q)), slip = rep(0.15, nrow(Q)))
#' sim <- simGDINA(N, Q, gs.parm = gs, model = "DINA", attribute = alpha)
#' fit <- GDINA(sim$dat, sim$Q, model = "DINA", verbose = 0)
#'
#' # Attribute-level logistic regression
#' ts_att <- ThreeStepCov(fit, ~ z_cont + z_bin + z_cat, data = Z, level = "attribute")
#' print(ts_att)
#'
#' # Profile-level multinomial logistic regression
#' ts_profile <- ThreeStepCov(
#'   fit,
#'   ~ z_cont + z_bin + z_cat,
#'   data = Z,
#'   level = "profile",
#'   reference = "first"
#' )
#' print(ts_profile)
#' }
#' @export
ThreeStepCov <- function(object,
                         formula,
                         data,
                         level = c("attribute", "profile"),
                         attribute = NULL,
                         classification = "MAP",
                         reference = c("last", "first"),
                         conf.level = 0.95) {
  if (!inherits(object, "GDINA")) {
    stop("object must be a GDINA object.", call. = FALSE)
  }

  level <- match.arg(level)
  reference <- match.arg(reference)

  design <- .three_step_design_matrix(
    formula = formula,
    data = data,
    nobs = extract(object, "nobs")
  )
  class_info <- .three_step_classification(object, classification)

  if (any(class_info$hard_pattern > 1)) {
    stop("ThreeStepCov is currently only available for binary attributes.", call. = FALSE)
  }

  if (level == "attribute") {
    misclassification <- CM(
      object,
      classification = class_info$classification,
      matrixtype = "attribute"
    )$att_classification
    natt <- ncol(class_info$hard_pattern)
    if (is.null(attribute)) {
      attribute <- seq_len(natt)
    }
    attribute <- as.integer(attribute)
    if (any(is.na(attribute)) || any(attribute < 1L) || any(attribute > natt)) {
      stop("attribute must contain valid attribute indices.", call. = FALSE)
    }

    results <- lapply(attribute, function(k) {
      observed <- class_info$hard_pattern[, k]
      naive_fit <- .three_step_binary_fit_core(
        design = design,
        observed = observed,
        misclassification = diag(2),
        conf.level = conf.level
      )
      corrected_fit <- .three_step_binary_fit_core(
        design = design,
        observed = observed,
        misclassification = misclassification[[k]],
        start = naive_fit$coefficients,
        conf.level = conf.level
      )

      list(
        attribute = k,
        observed = observed,
        misclassification = misclassification[[k]],
        naive = naive_fit,
        corrected = corrected_fit
      )
    })
    names(results) <- paste0("Attribute ", attribute)
  } else {
    misclassification <- CM(
      object,
      classification = class_info$classification,
      matrixtype = "profile"
    )$profile_classification
    naive_fit <- .three_step_multinomial_fit_core(
      design = design,
      observed = class_info$hard_class,
      misclassification = diag(ncol(misclassification)),
      class_labels = colnames(misclassification),
      reference = reference,
      conf.level = conf.level
    )
    corrected_fit <- .three_step_multinomial_fit_core(
      design = design,
      observed = class_info$hard_class,
      misclassification = misclassification,
      class_labels = colnames(misclassification),
      start = c(naive_fit$coefficients),
      reference = reference,
      conf.level = conf.level
    )
    results <- list(
      observed = class_info$hard_class,
      profile = class_info$hard_pattern,
      misclassification = misclassification,
      reference = reference,
      naive = naive_fit,
      corrected = corrected_fit
    )
  }

  ret <- list(
    call = match.call(),
    formula = formula,
    level = level,
    classification = classification,
    reference = if (level == "profile") reference else NULL,
    design = design,
    results = results
  )
  class(ret) <- c("ThreeStepCov")
  return(ret)
}


.three_step_cov_print_table <- function(x) {
  if (x$level == "attribute") {
    out <- do.call(rbind, lapply(x$results, function(result) {
      data.frame(
        attribute = paste0("Attribute ", result$attribute),
        result$corrected$table,
        row.names = NULL,
        check.names = FALSE
      )
    }))
    rownames(out) <- NULL
    return(out)
  }

  x$results$corrected$table
}

.three_step_format_print_table <- function(x, digits = 4) {
  numeric_columns <- vapply(x, is.numeric, logical(1L))
  x[numeric_columns] <- lapply(x[numeric_columns], round, digits = digits)
  x
}

#' @export
print.ThreeStepCov <- function(x, ...) {
  cat("\nThree-step covariate regression\n")
  cat("Level:", x$level, "\n")

  if (x$level == "attribute") {
    cat("Attributes:", paste(names(x$results), collapse = ", "), "\n")
  } else {
    cat("Reference profile:", x$results$corrected$reference, "\n")
  }

  cat("\nCorrected coefficients\n")
  print(.three_step_format_print_table(.three_step_cov_print_table(x)), row.names = FALSE)
  invisible(x)
}


.three_step_design_matrix <- function(formula, data, nobs) {
  if (missing(formula) || !inherits(formula, "formula")) {
    stop("formula must be a one-sided or two-sided formula.", call. = FALSE)
  }
  if (missing(data) || !is.data.frame(data)) {
    stop("data must be a data.frame.", call. = FALSE)
  }

  terms_obj <- stats::terms(formula, data = data)
  if (attr(terms_obj, "response") > 0L) {
    terms_obj <- stats::delete.response(terms_obj)
  }
  mf <- stats::model.frame(terms_obj, data = data, na.action = stats::na.fail)
  design <- stats::model.matrix(terms_obj, data = mf)

  if (nrow(design) != nobs) {
    stop("data must have one row per respondent in object.", call. = FALSE)
  }

  design
}

.three_step_classification <- function(object, classification) {
  pattern <- extract(object, "attributepattern")

  if (is.character(classification)) {
    classification <- match.arg(classification, c("MAP", "MLE", "EAP"))
  }

  classification_input <- classification

  if (is.character(classification) && length(classification) == 1L) {
    if (classification == "MAP") {
      hard_class <- max.col(extract(object, "logposterior.i"))
      hard_pattern <- pattern[hard_class, , drop = FALSE]
    } else if (classification == "MLE") {
      hard_class <- max.col(extract(object, "loglikelihood.i"))
      hard_pattern <- pattern[hard_class, , drop = FALSE]
    } else if (classification == "EAP"){
      hard_pattern <- as.matrix(personparm(object, what = "EAP"))
      hard_class <- matchMatrix(pattern, hard_pattern)
    }
  } else {
    if (!is.matrix(classification) && !is.data.frame(classification)) {
      stop("classification must be one of MAP, MLE, EAP, or a classification matrix.", call. = FALSE)
    }

    classification_input <- hard_pattern <- as.matrix(classification)
  }

  if (nrow(hard_pattern) != extract(object, "nobs")) {
    stop("classification must have one row per respondent.", call. = FALSE)
  }
  if (ncol(hard_pattern) != ncol(pattern)) {
    stop("classification must have one column per attribute.", call. = FALSE)
  }

  if (anyNA(hard_class)) {
    stop("classification must contain valid attribute patterns.", call. = FALSE)
  }

  list(
    pattern = pattern,
    hard_pattern = hard_pattern,
    hard_class = hard_class,
    classification = classification_input
  )
}

.three_step_safe_inverse <- function(mat, tol = sqrt(.Machine$double.eps)) {
  inv <- tryCatch(solve(mat), error = function(e) NULL)
  if (!is.null(inv) && all(is.finite(inv))) {
    return(inv)
  }
  MASS::ginv(mat, tol = tol)
}

.three_step_coef_table <- function(coefficients,
                                   vcov,
                                   term_names,
                                   conf.level = 0.95,
                                   class_labels = NULL) {
  std.error <- sqrt(pmax(diag(vcov), 0))
  statistic <- coefficients / std.error
  p.value <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  alpha <- 1 - conf.level
  crit <- stats::qnorm(1 - alpha / 2)
  conf.low <- coefficients - crit * std.error
  conf.high <- coefficients + crit * std.error
  odds.ratio <- exp(coefficients)
  or.conf.low <- exp(conf.low)
  or.conf.high <- exp(conf.high)

  if (is.null(class_labels)) {
    out <- data.frame(
      term = term_names,
      estimate = unname(coefficients),
      odds.ratio = unname(odds.ratio),
      std.error = unname(std.error),
      statistic = unname(statistic),
      p.value = unname(p.value),
      conf.low = unname(conf.low),
      conf.high = unname(conf.high),
      or.conf.low = unname(or.conf.low),
      or.conf.high = unname(or.conf.high),
      row.names = NULL
    )
  } else {
    out <- data.frame(
      class = rep(class_labels, each = length(term_names)),
      term = rep(term_names, times = length(class_labels)),
      estimate = unname(coefficients),
      odds.ratio = unname(odds.ratio),
      std.error = unname(std.error),
      statistic = unname(statistic),
      p.value = unname(p.value),
      conf.low = unname(conf.low),
      conf.high = unname(conf.high),
      or.conf.low = unname(or.conf.low),
      or.conf.high = unname(or.conf.high),
      row.names = NULL
    )
  }
  out
}

.three_step_binary_prob <- function(design, coefficients) {
  stats::plogis(drop(design %*% coefficients))
}

.three_step_binary_objective <- function(coefficients, design, observed, misclassification) {
  observed_index <- observed + 1L
  prob <- .three_step_binary_prob(design, coefficients)
  #L = P(W=s|X=0)P(X=0|Z) + P(W=s|X=1)P(X=1|Z)
  #s=0: L = P(W=0|X=0)(1-p) + P(W=0|X=1)p
  #s=1: L = P(W=1|X=0)(1-p) + P(W=1|X=1)p
  like <- misclassification[1, observed_index] * (1 - prob) +
    misclassification[2, observed_index] * prob
  -sum(log(pmax(like, .Machine$double.eps)))
}

.three_step_binary_fit_core <- function(design,
                                        observed,
                                        misclassification,
                                        start = NULL,
                                        conf.level = 0.95) {
  if (!all(observed %in% c(0, 1))) {
    stop("Attribute classifications must be coded as 0/1.", call. = FALSE)
  }

  if (is.null(start)) {
    init <- suppressWarnings(stats::glm.fit(
      x = design,
      y = observed,
      family = stats::binomial()
    ))
    start <- init$coefficients
    start[!is.finite(start)] <- 0
  }

  objective <- function(par) {
    .three_step_binary_objective(par, design, observed, misclassification)
  }

  opt <- stats::optim(
    start,
    objective,
    method = "BFGS",
    control = list(reltol = 1e-10, maxit = 1000)
  )
  hessian <- stats::optimHess(opt$par, objective)
  vcov <- .three_step_safe_inverse(hessian)
  fitted <- .three_step_binary_prob(design, opt$par)
  observed_index <- observed + 1L
  numerator <- misclassification[2, observed_index] * fitted
  denominator <- misclassification[1, observed_index] * (1 - fitted) + numerator
  posterior <- cbind(
    `P(X=0|W,Z)` = 1 - numerator / pmax(denominator, .Machine$double.eps),
    `P(X=1|W,Z)` = numerator / pmax(denominator, .Machine$double.eps)
  )

  names(opt$par) <- colnames(design)
  rownames(vcov) <- colnames(vcov) <- colnames(design)

  list(
    coefficients = opt$par,
    vcov = vcov,
    table = .three_step_coef_table(
      coefficients = opt$par,
      vcov = vcov,
      term_names = colnames(design),
      conf.level = conf.level
    ),
    logLik = -opt$value,
    converged = opt$convergence == 0L,
    counts = opt$counts,
    fitted = fitted,
    posterior = posterior
  )
}

.three_step_multinomial_prob <- function(design, coefficients, n_classes) {
  coefficient_matrix <- matrix(coefficients, nrow = ncol(design), ncol = n_classes - 1L)
  eta <- cbind(design %*% coefficient_matrix, 0)
  eta <- eta - apply(eta, 1, max)
  exp_eta <- exp(eta)
  exp_eta / rowSums(exp_eta)
}

.three_step_multinomial_objective <- function(coefficients,
                                              design,
                                              observed,
                                              misclassification) {
  misclassification <- misclassification / pmax(rowSums(misclassification), .Machine$double.eps)
  probability <- .three_step_multinomial_prob(
    design = design,
    coefficients = coefficients,
    n_classes = ncol(misclassification)
  )
  like <- rowSums(probability * t(misclassification[, observed, drop = FALSE]))
  -sum(log(pmax(like, .Machine$double.eps)))
}

.three_step_multinomial_fit_core <- function(design,
                                             observed,
                                             misclassification,
                                             class_labels = NULL,
                                             start = NULL,
                                             conf.level = 0.95,
                                             reference = c("last", "first")) {
  n_classes <- ncol(misclassification)
  reference <- match.arg(reference)
  if (length(observed) != nrow(design)) {
    stop("observed must have one entry per row of design.", call. = FALSE)
  }
  if (any(observed < 1L) || any(observed > n_classes)) {
    stop("Observed profile classifications are out of range.", call. = FALSE)
  }

  if (is.null(class_labels)) {
    class_labels <- colnames(misclassification)
    if (is.null(class_labels)) {
      class_labels <- paste0("Class ", seq_len(n_classes))
    }
  }

  internal_order <- if (reference == "last") {
    seq_len(n_classes)
  } else {
    c(seq_len(n_classes)[-1L], 1L)
  }
  back_order <- order(internal_order)
  observed_internal <- match(observed, internal_order)
  misclassification_internal <- misclassification[internal_order, internal_order, drop = FALSE]
  class_labels_internal <- class_labels[internal_order]

  n_par <- ncol(design) * (n_classes - 1L)
  if (is.null(start)) {
    start <- rep(0, n_par)
    class_prop <- tabulate(observed_internal, nbins = n_classes) / length(observed_internal)
    class_prop <- pmax(class_prop, .Machine$double.eps)
    start[seq(1, n_par, by = ncol(design))] <- log(class_prop[-n_classes] / class_prop[n_classes])
  }

  objective <- function(par) {
    .three_step_multinomial_objective(par, design, observed_internal, misclassification_internal)
  }

  opt <- stats::optim(
    start,
    objective,
    method = "BFGS",
    control = list(reltol = 1e-10, maxit = 1000)
  )
  hessian <- stats::optimHess(opt$par, objective)
  vcov <- .three_step_safe_inverse(hessian)
  probability_internal <- .three_step_multinomial_prob(design, opt$par, n_classes)
  posterior_internal <- probability_internal * t(misclassification_internal[, observed_internal, drop = FALSE])
  posterior_internal <- posterior_internal / pmax(rowSums(posterior_internal), .Machine$double.eps)

  probability <- probability_internal[, back_order, drop = FALSE]
  posterior <- posterior_internal[, back_order, drop = FALSE]

  coefficient_matrix <- matrix(
    opt$par,
    nrow = ncol(design),
    ncol = n_classes - 1L,
    dimnames = list(colnames(design), class_labels_internal[-n_classes])
  )
  rownames(vcov) <- colnames(vcov) <- as.vector(outer(colnames(design), class_labels_internal[-n_classes], paste, sep = ":"))

  list(
    coefficients = coefficient_matrix,
    vcov = vcov,
    table = .three_step_coef_table(
      coefficients = c(coefficient_matrix),
      vcov = vcov,
      term_names = colnames(design),
      conf.level = conf.level,
      class_labels = class_labels_internal[-n_classes]
    ),
    reference = class_labels_internal[n_classes],
    logLik = -opt$value,
    converged = opt$convergence == 0L,
    counts = opt$counts,
    fitted = stats::setNames(as.data.frame(probability), class_labels),
    posterior = stats::setNames(as.data.frame(posterior), class_labels)
  )
}
