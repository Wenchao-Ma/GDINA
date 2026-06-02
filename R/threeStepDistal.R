#' @title Three-step distal outcome analysis
#'
#' @description
#' Regress distal outcome on attribute classifications or latent
#' profile classifications from a fitted \code{\link{GDINA}} object using a
#' three-step correction. The corrected analysis uses the estimated
#' misclassification matrix \eqn{P(W=s\mid X=t)} from \code{\link{CM}} and can
#' be fit by maximum likelihood (ML) and/or the BCH weighting approach.
#'
#' @param object An estimated \code{GDINA} object returned from
#'   \code{\link{GDINA}}.
#' @param outcome A distal outcome with one value per respondent. Supported
#'   outcomes are binary, ordinal, nominal categorical, and continuous
#'   outcomes.
#' @param formula Optional one-sided or two-sided formula specifying observed
#'   covariates to include additively in the distal outcome regression. Only
#'   the right-hand side is used.
#' @param data Optional data frame containing the covariates in
#'   \code{formula}. It must contain one row per respondent in \code{object}.
#'   Ordered factors are converted to numeric scores in level order so they
#'   enter the regression with a single linear effect.
#' @param level A character string specifying whether the corrected analysis is
#'   carried out at the attribute level (\code{"attribute"}) or the latent
#'   profile level (\code{"profile"}).
#' @param attribute Optional integer vector giving the attributes to include in
#'   the attribute-level regression when \code{level = "attribute"}. Selected
#'   attributes enter the same regression model together. If \code{NULL}, all
#'   attributes are included.
#' @param classification A character string specifying the classification rule.
#'   Supported values are \code{"MAP"}, \code{"MLE"}, and \code{"EAP"}.
#'   Alternatively, a matrix of user-supplied classifications can be provided,
#'   with one row per respondent and one column per attribute.
#' @param method A character vector specifying which corrected distal analyses
#'   to return. Supported values are \code{"ML"} and \code{"BCH"}.
#' @param outcome_type A character string specifying the outcome type.
#'   Supported values are \code{"auto"}, \code{"binary"},
#'   \code{"ordinal"}, \code{"categorical"}, \code{"nominal"}, and
#'   \code{"continuous"}. The value \code{"nominal"} is treated as an
#'   alias for \code{"categorical"}.
#' @param reference A character string specifying the reference class for
#'   profile-level models. Supported values are \code{"last"} (default) and
#'   \code{"first"}. Ignored when \code{level = "attribute"}.
#' @param conf.level Confidence level used for Wald intervals.
#' @param maxit Maximum number of optimization iterations.
#' @param na.action See \code{\link{na.action}} for handling missing values in
#'   observed covariates.
#'
#' @details
#' Let \eqn{W} denote the estimated latent class, \eqn{X} the unobserved true
#' class, \eqn{Y} the distal outcome, and \eqn{Z} the observed covariates. The
#' ML correction is based on the
#' mixture likelihood
#' \deqn{\sum_i \log\left\{\sum_t P(W_i=s_i\mid X_i=t)P(X_i=t)f(Y_i\mid X_i=t, Z_i)\right\}.}
#' Observed covariates enter additively in the distal outcome model together
#' with the selected latent attributes or profiles. For binary distal outcomes,
#' \eqn{f(Y_i\mid X_i=t, Z_i)} is Bernoulli with a logit link; for continuous
#' outcomes it is Gaussian with state-specific means plus additive covariate
#' effects and a common residual standard deviation; for ordinal outcomes it is
#' a cumulative logit model with proportional odds; for nominal categorical
#' outcomes it is a multinomial logit model over the distal categories.
#'
#' The BCH correction replaces the full mixture likelihood with BCH weights
#' derived from the inverse misclassification matrix. Both the naive analysis
#' that ignores classification error and the requested corrected analyses are
#' returned.
#'
#' This function is currently limited to models with binary attributes.
#'
#' @return
#' A list containing the fitted distal outcome analyses. \code{results}
#' contains the requested distal outcome analysis for the selected attribute or
#' profile states. Each fit includes coefficient estimates, Wald summaries,
#' fitted class-specific outcome summaries, threshold estimates for ordinal
#' outcomes, any observed-covariate terms included in the model, and the
#' correction objects used for estimation.
#'
#' @references
#' Bakk, Z., & Kuha, J. (2021). Relating latent class membership to external
#' variables: An overview. \emph{British Journal of Mathematical and Statistical
#' Psychology, 74}(2), 340-362.
#'
#' Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent structure
#' models with categorical variables: One-step versus three-step estimators.
#' \emph{Political Analysis, 12}(1), 3-27.
#'
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved
#' three-step approaches. \emph{Political Analysis, 18}(4), 450-469.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' N <- 1000
#' Q <- sim10GDINA$simQ
#' gs <- data.frame(guess = rep(0.1, nrow(Q)), slip = rep(0.24, nrow(Q)))
#'
#' sim <- simGDINA(N, Q, gs.parm = gs, model = "GDINA")
#' fit <- GDINA(sim$dat, sim$Q, verbose = 0)
#' alpha <- extract(sim, "attribute")
#'
#' # Binary distal outcome at the attribute level
#' y_binary <- rbinom(N, 1, plogis(-0.8 + alpha[, 1] - 0.9 * alpha[, 2] + 0.7 * alpha[, 3]))
#' distal_att <- ThreeStepDistal(fit, y_binary, level = "attribute")
#' print(distal_att)
#'
#' # Continuous distal outcome at the profile level
#' y_continuous <- 0.3 + 0.9 * alpha[, 1] - 1.1 * alpha[, 2] + 0.8 * alpha[, 3] + rnorm(N, sd = 1.1)
#' distal_cont <- ThreeStepDistal(
#'   fit,
#'   y_continuous,
#'   level = "profile",
#'   method = "BCH",
#'   reference = "first"
#' )
#' distal_cont
#'
#' # Binary distal outcome with observed covariates. Ordered factors enter as
#' # numeric scores in level order to mimic a single linear effect.
#' covariates <- data.frame(
#'   gender = factor(sample(c("female", "male"), N, replace = TRUE)),
#'   ses = ordered(sample(c("low", "medium", "high"), N, replace = TRUE),
#'     levels = c("low", "medium", "high")
#'   ),
#'   age = rnorm(N)
#' )
#' y_cov <- rbinom(
#'   N,
#'   1,
#'   plogis(-0.4 + 0.8 * alpha[, 1] - 0.5 * alpha[, 2] + 0.5 * (covariates$gender == "male") +
#'     0.35 * as.numeric(covariates$ses) + 0.2 * covariates$age)
#' )
#' distal_cov <- ThreeStepDistal(
#'   fit,
#'   y_cov,
#'   formula = ~ gender + ses + age,
#'   data = covariates,
#'   level = "attribute",
#'   attribute = 1:3,
#'   method = c("ML", "BCH")
#' )
#' distal_cov
#'
#' # Ordinal distal outcome at the attribute level
#' latent_score <- -0.4 + 1.1 * alpha[, 1] - 0.8 * alpha[, 2] + 0.5 * alpha[, 3] + rnorm(N, sd = 0.8)
#' cuts <- quantile(latent_score, probs = c(1 / 3, 2 / 3))
#' y_ordinal <- ordered(
#'   cut(latent_score, breaks = c(-Inf, cuts, Inf), labels = c("low", "medium", "high")),
#'   levels = c("low", "medium", "high")
#' )
#' distal_ord <- ThreeStepDistal(
#'   fit,
#'   y_ordinal,
#'   level = "attribute",
#'   outcome_type = "ordinal",
#'   method = "BCH"
#' )
#' distal_ord
#' }
#' @export
ThreeStepDistal <- function(object,
                            outcome,
                            formula = NULL,
                            data = NULL,
                            level = c("attribute", "profile"),
                            attribute = NULL,
                            classification = "MAP",
                            method = c("BCH","ML"),
                            outcome_type = c("auto", "binary", "ordinal", "categorical", "nominal", "continuous"),
                            reference = c("last", "first"),
                            conf.level = 0.95,
                            maxit = 1000,
                            na.action = getOption("na.action")) {
  level <- match.arg(level)
  method <- unique(match.arg(method, c("ML", "BCH"), several.ok = TRUE))
  reference <- match.arg(reference)
  inputs <- .distal_prepare_inputs(
    object = object,
    outcome = outcome,
    outcome_type = outcome_type,
    classification = classification,
    formula = formula,
    data = data,
    na.action = na.action
  )

  if (level == "attribute") {
    natt <- ncol(inputs$classification$hard_pattern)
    if (is.null(attribute)) {
      attribute <- seq_len(natt)
    }
    attribute <- unique(as.integer(attribute))
    if (length(attribute) == 0L || any(is.na(attribute)) || any(attribute < 1L) || any(attribute > natt)) {
      stop("attribute must contain valid attribute indices.", call. = FALSE)
    }

    state_info <- .distal_attribute_state_info(
      pattern = inputs$pattern,
      hard_pattern = inputs$classification$hard_pattern,
      posterior = inputs$posterior,
      attribute = attribute
    )
    design_info <- .distal_state_design(
      level = "attribute",
      n_states = nrow(state_info$state_pattern),
      labels = state_info$state_labels,
      attribute = attribute,
      state_pattern = state_info$state_pattern
    )
    fits <- .distal_fit_bundle(
      outcome_info = inputs$outcome,
      state_design = design_info$design,
      covariate_design = inputs$covariates$design,
      observed_state = state_info$hard_class,
      misclassification = state_info$misclassification,
      prior = state_info$prior,
      method = method,
      conf.level = conf.level,
      maxit = maxit
    )

    results <- c(list(
      attribute = attribute,
      observed = state_info$hard_pattern,
      baseline = design_info$baseline,
      prior = stats::setNames(state_info$prior, design_info$state_labels),
      misclassification = state_info$misclassification,
      predictor_labels = .distal_predictor_term_names(design_info$design, inputs$covariates$design),
      state_labels = design_info$state_labels,
      covariate_terms = colnames(inputs$covariates$design)
    ), fits)
  } else {
    misclassification <- CM(
      object,
      classification = inputs$classification$classification,
      matrixtype = "profile"
    )$profile_classification
    state_labels <- colnames(misclassification)
    if (is.null(state_labels)) {
      state_labels <- apply(inputs$pattern, 1L, paste, collapse = "")
    }
    design_info <- .distal_state_design(
      level = "profile",
      n_states = nrow(misclassification),
      labels = state_labels,
      reference = reference
    )
    fits <- .distal_fit_bundle(
      outcome_info = inputs$outcome,
      state_design = design_info$design,
      covariate_design = inputs$covariates$design,
      observed_state = inputs$classification$hard_class,
      misclassification = misclassification,
      prior = colMeans(inputs$posterior),
      method = method,
      conf.level = conf.level,
      maxit = maxit
    )

    results <- c(list(
      observed = inputs$classification$hard_class,
      profile = inputs$classification$hard_pattern,
      baseline = design_info$baseline,
      prior = stats::setNames(colMeans(inputs$posterior), design_info$state_labels),
      misclassification = misclassification,
      predictor_labels = .distal_predictor_term_names(design_info$design, inputs$covariates$design),
      reference = reference,
      covariate_terms = colnames(inputs$covariates$design)
    ), fits)
  }

  ret <- list(
    call = match.call(),
    formula = formula,
    level = level,
    classification = classification,
    method = method,
    outcome_type = inputs$outcome$type,
    outcome_levels = inputs$outcome$levels,
    na.action = inputs$covariates$na.action,
    covariate_design = inputs$covariates$design,
    results = results
  )
  class(ret) <- c("ThreeStepDistal", "threeStepDistal")
  return(ret)
}


.distal_print_method <- function(x) {
  if ("ML" %in% x$method && "ML" %in% names(x$results)) {
    return("ML")
  }
  if ("BCH" %in% x$method && "BCH" %in% names(x$results)) {
    return("BCH")
  }
  if ("naive" %in% names(x$results)) {
    return("naive")
  }
  names(x$results)[vapply(x$results, is.list, logical(1))][1L]
}

#' @export
print.ThreeStepDistal <- function(x, ...) {
  fit_name <- .distal_print_method(x)
  fit <- x$results[[fit_name]]

  cat("\nThree-step distal outcome analysis\n")
  cat("Level:", x$level, "\n")
  cat("Outcome type:", x$outcome_type, "\n")
  cat("Correction methods:", paste(x$method, collapse = ", "), "\n")
  cat("Displayed fit:", fit_name, "\n")

  if (x$level == "attribute") {
    cat("Attributes in regression:", paste0("Attribute ", x$results$attribute, collapse = ", "), "\n")
  } else {
    cat("Reference profile:", x$results$baseline, "\n")
  }

  if (!is.null(x$formula)) {
    cat("Observed covariates:", paste(attr(stats::terms(x$formula), "term.labels"), collapse = ", "), "\n")
  }

  if (!is.null(fit$table)) {
    cat("\nCoefficient table\n")
    print(.three_step_format_print_table(fit$table), row.names = FALSE)
  }

  if (!is.null(fit$threshold_table)) {
    cat("\nThreshold table\n")
    print(.three_step_format_print_table(fit$threshold_table), row.names = FALSE)
  }

  invisible(x)
}


.distal_detect_outcome_type <- function(outcome,
                                        outcome_type = c("auto", "binary", "ordinal", "categorical", "nominal", "continuous")) {
  outcome_type <- match.arg(outcome_type)

  if (outcome_type == "nominal") {
    outcome_type <- "categorical"
  }

  if (outcome_type != "auto") {
    return(outcome_type)
  }

  if (is.ordered(outcome)) {
    n_levels <- nlevels(outcome)
    return(if (n_levels == 2L) "binary" else "ordinal")
  }

  if (is.factor(outcome) || is.character(outcome)) {
    n_levels <- nlevels(factor(outcome))
    return(if (n_levels == 2L) "binary" else "categorical")
  }

  if (is.logical(outcome)) {
    return("binary")
  }

  if (is.numeric(outcome) && all(stats::na.omit(unique(outcome)) %in% c(0, 1))) {
    return("binary")
  }

  "continuous"
}

.distal_ordinal_threshold_names <- function(levels) {
  paste(utils::head(levels, -1L), utils::tail(levels, -1L), sep = "|")
}

.distal_make_strictly_increasing <- function(x, increment = 1e-4) {
  if (length(x) <= 1L) {
    return(x)
  }

  for (idx in 2:length(x)) {
    if (x[idx] <= x[idx - 1L]) {
      x[idx] <- x[idx - 1L] + increment
    }
  }

  x
}

.distal_ordinal_predictor_design <- function(state_design) {
  state_design[, -1L, drop = FALSE]
}

.distal_predictor_term_names <- function(state_design,
                                         covariate_design = NULL,
                                         ordinal = FALSE) {
  state_terms <- colnames(state_design)
  if (ordinal) {
    state_terms <- state_terms[-1L]
  }
  covariate_terms <- if (is.null(covariate_design) || ncol(covariate_design) == 0L) {
    character(0)
  } else {
    colnames(covariate_design)
  }

  c(state_terms, covariate_terms)
}

.distal_observed_predictor_design <- function(state_design,
                                              observed_state,
                                              covariate_design = NULL,
                                              ordinal = FALSE) {
  observed_design <- state_design[observed_state, , drop = FALSE]
  if (ordinal) {
    observed_design <- observed_design[, -1L, drop = FALSE]
  }
  if (!is.null(covariate_design) && ncol(covariate_design) > 0L) {
    observed_design <- cbind(observed_design, covariate_design)
  }
  observed_design
}

.distal_covariate_linear_predictor <- function(covariate_design, beta_covariate, n_obs, n_states) {
  if (is.null(covariate_design) || ncol(covariate_design) == 0L || length(beta_covariate) == 0L) {
    return(matrix(0, nrow = n_obs, ncol = n_states))
  }

  covariate_eta <- drop(covariate_design %*% beta_covariate)
  matrix(covariate_eta, nrow = n_obs, ncol = n_states)
}

.distal_state_linear_predictor <- function(state_design, beta_state, n_obs) {
  state_eta <- drop(state_design %*% beta_state)
  matrix(state_eta, nrow = n_obs, ncol = nrow(state_design), byrow = TRUE)
}

.distal_ordered_to_numeric <- function(data) {
  data_converted <- data
  ordered_columns <- vapply(data_converted, is.ordered, logical(1L))
  if (any(ordered_columns)) {
    data_converted[ordered_columns] <- lapply(data_converted[ordered_columns], as.numeric)
  }
  data_converted
}

.distal_covariate_design <- function(formula,
                                     data,
                                     nobs,
                                     na.action = getOption("na.action")) {
  if (is.null(formula) && is.null(data)) {
    return(list(
      design = matrix(numeric(0), nrow = nobs, ncol = 0L),
      keep = seq_len(nobs),
      na.action = NULL
    ))
  }

  if (is.null(formula) || !inherits(formula, "formula")) {
    stop("formula must be NULL or a one-sided/two-sided formula.", call. = FALSE)
  }
  if (is.null(data) || !is.data.frame(data)) {
    stop("data must be supplied as a data.frame when formula is used.", call. = FALSE)
  }
  if (nrow(data) != nobs) {
    stop("data must have one row per respondent in object.", call. = FALSE)
  }

  data_converted <- .distal_ordered_to_numeric(data)
  terms_obj <- stats::terms(formula, data = data_converted)
  if (attr(terms_obj, "response") > 0L) {
    terms_obj <- stats::delete.response(terms_obj)
  }
  mf <- stats::model.frame(terms_obj, data = data_converted, na.action = na.action)
  design <- stats::model.matrix(terms_obj, data = mf)
  if ("(Intercept)" %in% colnames(design)) {
    design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  }
  omitted <- attr(mf, "na.action")
  keep <- if (is.null(omitted)) seq_len(nobs) else setdiff(seq_len(nobs), as.integer(omitted))

  if (nrow(design) != length(keep)) {
    stop("Failed to align observed covariate rows with retained respondents.", call. = FALSE)
  }

  list(
    design = design,
    keep = keep,
    na.action = omitted
  )
}

.distal_ordinal_decode_thresholds <- function(raw_thresholds) {
  if (length(raw_thresholds) == 0L) {
    return(numeric(0))
  }

  thresholds <- raw_thresholds
  if (length(raw_thresholds) > 1L) {
    thresholds[2:length(raw_thresholds)] <- thresholds[1L] + cumsum(exp(raw_thresholds[2:length(raw_thresholds)]))
  }
  thresholds
}

.distal_ordinal_encode_thresholds <- function(thresholds) {
  thresholds <- .distal_make_strictly_increasing(thresholds)
  c(thresholds[1L], if (length(thresholds) > 1L) log(diff(thresholds)) else numeric(0))
}

.distal_ordinal_parameters <- function(theta,
                                       outcome_info,
                                       state_design,
                                       covariate_design = NULL) {
  n_thresholds <- length(outcome_info$levels) - 1L
  raw_thresholds <- theta[seq_len(n_thresholds)]
  thresholds <- .distal_ordinal_decode_thresholds(raw_thresholds)
  predictor_terms <- .distal_predictor_term_names(state_design, covariate_design, ordinal = TRUE)
  beta <- theta[-seq_len(n_thresholds)]

  names(thresholds) <- .distal_ordinal_threshold_names(outcome_info$levels)
  names(beta) <- predictor_terms

  list(
    thresholds = thresholds,
    raw_thresholds = raw_thresholds,
    beta = beta,
    predictor_terms = predictor_terms
  )
}

.distal_ordinal_start <- function(outcome_info, predictor_design) {
  y <- ordered(outcome_info$levels[outcome_info$y], levels = outcome_info$levels)
  n_thresholds <- length(outcome_info$levels) - 1L
  threshold_names <- .distal_ordinal_threshold_names(outcome_info$levels)
  beta <- rep(0, ncol(predictor_design))

  if (ncol(predictor_design) > 0L) {
    predictor_df <- stats::setNames(as.data.frame(predictor_design), paste0("x", seq_len(ncol(predictor_design))))
    start_fit <- try(suppressWarnings(MASS::polr(
      formula = y ~ .,
      data = cbind(data.frame(y = y), predictor_df),
      method = "logistic",
      Hess = FALSE
    )), silent = TRUE)

    if (!inherits(start_fit, "try-error")) {
      beta_fit <- stats::coef(start_fit)
      if (length(beta_fit) > 0L) {
        beta_index <- as.integer(sub("^x", "", names(beta_fit)))
        beta[beta_index] <- unname(beta_fit)
      }
      thresholds <- unname(start_fit$zeta)
      theta <- c(.distal_ordinal_encode_thresholds(thresholds), beta)
      names(theta) <- c(threshold_names, colnames(predictor_design))
      return(theta)
    }
  }

  prop <- tabulate(outcome_info$y, nbins = length(outcome_info$levels)) / length(outcome_info$y)
  prop <- pmax(prop, .Machine$double.eps)
  prop <- prop / sum(prop)
  thresholds <- stats::qlogis(cumsum(prop)[seq_len(n_thresholds)])
  thresholds <- .distal_make_strictly_increasing(thresholds)
  theta <- c(.distal_ordinal_encode_thresholds(thresholds), beta)
  names(theta) <- c(threshold_names, colnames(predictor_design))
  theta
}

.distal_ordinal_threshold_vcov <- function(raw_thresholds, vcov_threshold_raw) {
  n_thresholds <- length(raw_thresholds)
  jacobian <- matrix(0, nrow = n_thresholds, ncol = n_thresholds)
  jacobian[, 1L] <- 1

  if (n_thresholds > 1L) {
    for (idx in 2:n_thresholds) {
      jacobian[idx:n_thresholds, idx] <- exp(raw_thresholds[idx])
    }
  }

  jacobian %*% vcov_threshold_raw %*% t(jacobian)
}

.distal_ordinal_threshold_table <- function(thresholds, vcov, conf.level = 0.95) {
  std_error <- sqrt(pmax(diag(vcov), 0))
  statistic <- thresholds / pmax(std_error, .Machine$double.eps)
  p_value <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  alpha <- 1 - conf.level
  crit <- stats::qnorm(1 - alpha / 2)
  conf.low <- thresholds - crit * std_error
  conf.high <- thresholds + crit * std_error

  data.frame(
    threshold = names(thresholds),
    estimate = unname(thresholds),
    std.error = unname(std_error),
    statistic = unname(statistic),
    p.value = unname(p_value),
    conf.low = unname(conf.low),
    conf.high = unname(conf.high),
    row.names = NULL
  )
}

.distal_normalize_outcome <- function(outcome, outcome_type) {
  if (outcome_type == "nominal") {
    outcome_type <- "categorical"
  }

  if (length(outcome) == 0L) {
    stop("outcome must contain at least one observation.", call. = FALSE)
  }
  if (anyNA(outcome)) {
    stop("outcome cannot contain missing values.", call. = FALSE)
  }

  if (outcome_type == "binary") {
    if (is.factor(outcome) || is.character(outcome)) {
      outcome <- factor(outcome)
      if (nlevels(outcome) != 2L) {
        stop("Binary outcomes must have exactly two levels.", call. = FALSE)
      }
      y <- as.integer(outcome == levels(outcome)[2L])
      return(list(
        type = "binary",
        y = y,
        levels = levels(outcome),
        baseline = levels(outcome)[1L],
        event = levels(outcome)[2L]
      ))
    }

    y <- if (is.logical(outcome)) as.integer(outcome) else as.numeric(outcome)
    if (!all(y %in% c(0, 1))) {
      stop("Binary outcomes must be coded as 0/1, logical, or a two-level factor.", call. = FALSE)
    }

    return(list(
      type = "binary",
      y = y,
      levels = c("0", "1"),
      baseline = "0",
      event = "1"
    ))
  }

  if (outcome_type == "ordinal") {
    y_factor <- if (is.ordered(outcome)) {
      outcome
    } else if (is.factor(outcome)) {
      ordered(outcome, levels = levels(outcome))
    } else {
      ordered(outcome)
    }

    if (nlevels(y_factor) < 3L) {
      stop("Ordinal outcomes must have at least three ordered levels.", call. = FALSE)
    }

    return(list(
      type = "ordinal",
      y = as.integer(y_factor),
      levels = levels(y_factor)
    ))
  }

  if (outcome_type == "categorical") {
    y_factor <- factor(outcome)
    if (nlevels(y_factor) < 3L) {
      stop("Categorical outcomes must have at least three levels.", call. = FALSE)
    }

    return(list(
      type = "categorical",
      y = as.integer(y_factor),
      levels = levels(y_factor),
      baseline = levels(y_factor)[nlevels(y_factor)]
    ))
  }

  y <- as.numeric(outcome)
  if (!all(is.finite(y))) {
    stop("Continuous outcomes must be numeric and finite.", call. = FALSE)
  }

  list(
    type = "continuous",
    y = y,
    mean = mean(y),
    sd = stats::sd(y)
  )
}

.distal_coef_table <- function(coefficients,
                               vcov,
                               term_names,
                               outcome_info,
                               conf.level = 0.95) {
  if (outcome_info$type != "categorical") {
    return(.three_step_coef_table(
      coefficients = coefficients,
      vcov = vcov,
      term_names = term_names,
      conf.level = conf.level
    ))
  }

  std_error <- sqrt(pmax(diag(vcov), 0))
  statistic <- coefficients / pmax(std_error, .Machine$double.eps)
  p_value <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  alpha <- 1 - conf.level
  crit <- stats::qnorm(1 - alpha / 2)
  conf.low <- coefficients - crit * std_error
  conf.high <- coefficients + crit * std_error
  odds.ratio <- exp(coefficients)
  or.conf.low <- exp(conf.low)
  or.conf.high <- exp(conf.high)
  target_levels <- outcome_info$levels[-length(outcome_info$levels)]

  data.frame(
    outcome_level = rep(target_levels, each = length(term_names)),
    term = rep(term_names, times = length(target_levels)),
    estimate = unname(coefficients),
    odds.ratio = unname(odds.ratio),
    std.error = unname(std_error),
    statistic = unname(statistic),
    p.value = unname(p_value),
    conf.low = unname(conf.low),
    conf.high = unname(conf.high),
    or.conf.low = unname(or.conf.low),
    or.conf.high = unname(or.conf.high),
    row.names = NULL
  )
}

.distal_state_design <- function(level,
                                 n_states,
                                 labels = NULL,
                                 reference = c("last", "first"),
                                 attribute = NULL,
                                 state_pattern = NULL) {
  reference <- match.arg(reference)

  if (level == "attribute") {
    if (is.null(attribute) || is.null(state_pattern)) {
      stop("attribute and state_pattern are required for attribute-level designs.", call. = FALSE)
    }
    design <- cbind(`(Intercept)` = 1, state_pattern)
    colnames(design) <- c("(Intercept)", paste0("Attribute ", attribute))
    if (is.null(labels)) {
      labels <- apply(state_pattern, 1L, paste, collapse = "")
    }
    baseline_index <- which(rowSums(state_pattern) == 0)
    if (length(baseline_index) == 0L) {
      baseline_index <- 1L
    } else {
      baseline_index <- baseline_index[1L]
    }
    return(list(
      design = design,
      baseline = labels[baseline_index],
      state_labels = labels,
      nonreference = colnames(design)[-1L]
    ))
  }

  if (is.null(labels)) {
    labels <- paste0("Class ", seq_len(n_states))
  }

  baseline_index <- if (reference == "last") n_states else 1L
  nonreference <- setdiff(seq_len(n_states), baseline_index)
  design <- cbind(`(Intercept)` = 1, diag(n_states)[, nonreference, drop = FALSE])
  colnames(design) <- c("(Intercept)", labels[nonreference])

  list(
    design = design,
    baseline = labels[baseline_index],
    state_labels = labels,
    nonreference = labels[nonreference]
  )
}

.distal_match_rows <- function(x, table) {
  x_key <- apply(as.matrix(x), 1L, paste, collapse = "\r")
  table_key <- apply(as.matrix(table), 1L, paste, collapse = "\r")
  match(x_key, table_key)
}

.distal_attribute_state_info <- function(pattern,
                                         hard_pattern,
                                         posterior,
                                         attribute) {
  state_pattern <- pattern[, attribute, drop = FALSE]
  state_pattern <- state_pattern[!duplicated(state_pattern), , drop = FALSE]
  state_labels <- apply(state_pattern, 1L, paste, collapse = "")

  profile_to_state <- .distal_match_rows(pattern[, attribute, drop = FALSE], state_pattern)
  posterior_state <- vapply(seq_len(nrow(state_pattern)), function(idx) {
    rowSums(posterior[, profile_to_state == idx, drop = FALSE])
  }, numeric(nrow(posterior)))
  if (!is.matrix(posterior_state)) {
    posterior_state <- matrix(posterior_state, ncol = 1L)
  }

  hard_pattern <- hard_pattern[, attribute, drop = FALSE]
  hard_class <- .distal_match_rows(hard_pattern, state_pattern)
  observed_assignment <- diag(nrow(state_pattern))[hard_class, , drop = FALSE]
  misclassification <- crossprod(posterior_state, observed_assignment)
  misclassification <- misclassification / pmax(rowSums(misclassification), .Machine$double.eps)
  colnames(misclassification) <- rownames(misclassification) <- state_labels

  list(
    state_pattern = state_pattern,
    state_labels = state_labels,
    posterior = posterior_state,
    hard_pattern = hard_pattern,
    hard_class = hard_class,
    prior = colMeans(posterior_state),
    misclassification = misclassification
  )
}

.distal_prepare_inputs <- function(object,
                                   outcome,
                                   outcome_type = c("auto", "binary", "ordinal", "categorical", "nominal", "continuous"),
                                   classification = "MAP",
                                   formula = NULL,
                                   data = NULL,
                                   na.action = getOption("na.action")) {
  if (!inherits(object, "GDINA")) {
    stop("object must be a GDINA object.", call. = FALSE)
  }

  pattern <- extract(object, "attributepattern")
  if (any(pattern > 1)) {
    stop("threeStepDistal currently supports only binary attributes.", call. = FALSE)
  }

  nobs <- extract(object, "nobs")
  covariates <- .distal_covariate_design(
    formula = formula,
    data = data,
    nobs = nobs,
    na.action = na.action
  )
  keep <- covariates$keep

  outcome_type <- .distal_detect_outcome_type(outcome, outcome_type)
  if (length(outcome) != nobs) {
    stop("outcome must have one value per respondent in object.", call. = FALSE)
  }
  outcome <- outcome[keep]
  outcome_info <- .distal_normalize_outcome(outcome, outcome_type)

  posterior <- exp(indlogPost(object))
  posterior <- posterior / rowSums(posterior)
  posterior <- posterior[keep, , drop = FALSE]
  classification_info <- .three_step_classification(object, classification)
  classification_info$hard_pattern <- classification_info$hard_pattern[keep, , drop = FALSE]
  classification_info$hard_class <- classification_info$hard_class[keep]

  list(
    outcome = outcome_info,
    posterior = posterior,
    classification = classification_info,
    pattern = pattern,
    mp = personparm(object, what = "mp")[keep, , drop = FALSE],
    covariates = covariates
  )
}

.distal_initial_theta <- function(outcome_info,
                                  state_design,
                                  observed_state,
                                  covariate_design = NULL) {
  observed_design <- .distal_observed_predictor_design(state_design, observed_state, covariate_design)
  p <- ncol(observed_design)

  if (outcome_info$type == "binary") {
    fit <- suppressWarnings(stats::glm.fit(
      x = observed_design,
      y = outcome_info$y,
      family = stats::binomial()
    ))
    theta <- fit$coefficients
    theta[!is.finite(theta)] <- 0
    names(theta) <- colnames(state_design)
    return(theta)
  }

  if (outcome_info$type == "continuous") {
    fit <- stats::lm.wfit(x = observed_design, y = outcome_info$y, w = rep(1, length(outcome_info$y)))
    beta <- fit$coefficients
    beta[!is.finite(beta)] <- 0
    residual <- outcome_info$y - drop(observed_design %*% beta)
    sigma <- sqrt(mean(residual^2))
    sigma <- max(sigma, .Machine$double.eps)
    theta <- c(beta, log_sigma = log(sigma))
    names(theta) <- c(colnames(state_design), "log_sigma")
    return(theta)
  }

  if (outcome_info$type == "ordinal") {
    observed_ordinal_design <- .distal_observed_predictor_design(
      state_design = state_design,
      observed_state = observed_state,
      covariate_design = covariate_design,
      ordinal = TRUE
    )
    return(.distal_ordinal_start(outcome_info, observed_ordinal_design))
  }

  n_outcome <- length(outcome_info$levels)
  prop <- tabulate(outcome_info$y, nbins = n_outcome) / length(outcome_info$y)
  prop <- pmax(prop, .Machine$double.eps)
  coef_mat <- matrix(0, nrow = p, ncol = n_outcome - 1L)
  coef_mat[1L, ] <- log(prop[-n_outcome] / prop[n_outcome])
  theta <- c(coef_mat)
  names(theta) <- as.vector(outer(colnames(state_design), outcome_info$levels[-n_outcome], paste, sep = ":"))
  theta
}

.distal_likelihood_matrix <- function(theta,
                                      outcome_info,
                                      state_design,
                                      covariate_design = NULL) {
  n_obs <- length(outcome_info$y)
  n_states <- nrow(state_design)
  eps <- .Machine$double.eps
  n_state_terms <- ncol(state_design)
  n_covariates <- if (is.null(covariate_design)) 0L else ncol(covariate_design)

  if (outcome_info$type == "binary") {
    beta_state <- theta[seq_len(n_state_terms)]
    beta_covariate <- if (n_covariates > 0L) theta[n_state_terms + seq_len(n_covariates)] else numeric(0)
    eta <- .distal_state_linear_predictor(state_design, beta_state, n_obs) +
      .distal_covariate_linear_predictor(covariate_design, beta_covariate, n_obs, n_states)
    prob_mat <- pmin(pmax(stats::plogis(eta), eps), 1 - eps)
    like <- ifelse(matrix(outcome_info$y, nrow = n_obs, ncol = n_states) == 1, prob_mat, 1 - prob_mat)

    return(list(
      like = pmax(like, eps),
      state_fitted = if (n_covariates == 0L) prob_mat[1L, ] else prob_mat,
      state_summary = data.frame(state = seq_len(n_states), event_probability = colMeans(prob_mat))
    ))
  }

  if (outcome_info$type == "continuous") {
    beta_state <- theta[seq_len(n_state_terms)]
    beta_covariate <- if (n_covariates > 0L) theta[n_state_terms + seq_len(n_covariates)] else numeric(0)
    sigma <- exp(theta[length(theta)])
    sigma <- max(sigma, eps)
    mu <- .distal_state_linear_predictor(state_design, beta_state, n_obs) +
      .distal_covariate_linear_predictor(covariate_design, beta_covariate, n_obs, n_states)
    like <- matrix(
      stats::dnorm(
        x = rep(outcome_info$y, each = n_states),
        mean = c(t(mu)),
        sd = sigma
      ),
      nrow = n_obs,
      byrow = TRUE
    )

    return(list(
      like = pmax(like, eps),
      state_fitted = if (n_covariates == 0L) mu[1L, ] else mu,
      state_summary = data.frame(state = seq_len(n_states), mean = colMeans(mu), sigma = sigma)
    ))
  }

  if (outcome_info$type == "ordinal") {
    ordinal_par <- .distal_ordinal_parameters(theta, outcome_info, state_design, covariate_design)
    n_state_predictors <- ncol(.distal_ordinal_predictor_design(state_design))
    beta_state <- if (n_state_predictors > 0L) ordinal_par$beta[seq_len(n_state_predictors)] else numeric(0)
    beta_covariate <- if (n_covariates > 0L) ordinal_par$beta[n_state_predictors + seq_len(n_covariates)] else numeric(0)
    state_predictor <- if (n_state_predictors > 0L) {
      .distal_state_linear_predictor(.distal_ordinal_predictor_design(state_design), beta_state, n_obs)
    } else {
      matrix(0, nrow = n_obs, ncol = n_states)
    }
    eta <- state_predictor + .distal_covariate_linear_predictor(covariate_design, beta_covariate, n_obs, n_states)
    prob_array <- array(NA_real_, dim = c(n_obs, n_states, length(outcome_info$levels)))
    like <- matrix(NA_real_, nrow = n_obs, ncol = n_states)
    for (state_index in seq_len(n_states)) {
      cumulative <- stats::plogis(outer(eta[, state_index], ordinal_par$thresholds, function(e, z) z - e))
      cumulative <- cbind(0, cumulative, 1)
      prob <- cumulative[, -1L, drop = FALSE] - cumulative[, -ncol(cumulative), drop = FALSE]
      prob <- pmin(pmax(prob, eps), 1 - eps)
      prob <- prob / rowSums(prob)
      colnames(prob) <- outcome_info$levels
      prob_array[, state_index, ] <- prob
      like[, state_index] <- prob[cbind(seq_len(n_obs), outcome_info$y)]
    }

    return(list(
      like = pmax(like, eps),
      state_fitted = if (n_covariates == 0L) prob_array[1L, , , drop = FALSE][1L, , ] else prob_array,
      state_summary = stats::setNames(as.data.frame(apply(prob_array, c(2, 3), mean)), outcome_info$levels),
      thresholds = ordinal_par$thresholds,
      raw_thresholds = ordinal_par$raw_thresholds
    ))
  }

  n_outcome <- length(outcome_info$levels)
  coef_mat <- matrix(theta, nrow = n_state_terms + n_covariates, ncol = n_outcome - 1L)
  eta <- array(0, dim = c(n_obs, n_states, n_outcome))
  for (outcome_index in seq_len(n_outcome - 1L)) {
    beta_state <- coef_mat[seq_len(n_state_terms), outcome_index]
    beta_covariate <- if (n_covariates > 0L) coef_mat[n_state_terms + seq_len(n_covariates), outcome_index] else numeric(0)
    eta[, , outcome_index] <- .distal_state_linear_predictor(state_design, beta_state, n_obs) +
      .distal_covariate_linear_predictor(covariate_design, beta_covariate, n_obs, n_states)
  }
  eta <- sweep(eta, c(1, 2), apply(eta, c(1, 2), max), FUN = "-")
  prob_array <- exp(eta)
  prob_array <- sweep(prob_array, c(1, 2), apply(prob_array, c(1, 2), sum), FUN = "/")
  like <- matrix(NA_real_, nrow = n_obs, ncol = n_states)
  for (state_index in seq_len(n_states)) {
    like[, state_index] <- prob_array[cbind(seq_len(n_obs), rep(state_index, n_obs), outcome_info$y)]
  }

  list(
    like = pmax(like, eps),
    state_fitted = if (n_covariates == 0L) prob_array[1L, , , drop = FALSE][1L, , ] else prob_array,
    state_summary = stats::setNames(as.data.frame(apply(prob_array, c(2, 3), mean)), outcome_info$levels)
  )
}

.distal_unpack_coefficients <- function(theta,
                                        outcome_info,
                                        state_design,
                                        covariate_design = NULL) {
  term_names <- .distal_predictor_term_names(state_design, covariate_design)
  if (outcome_info$type == "categorical") {
    coef_mat <- matrix(theta, nrow = length(term_names), ncol = length(outcome_info$levels) - 1L)
    rownames(coef_mat) <- term_names
    colnames(coef_mat) <- outcome_info$levels[-length(outcome_info$levels)]
    return(coef_mat)
  }

  if (outcome_info$type == "ordinal") {
    return(.distal_ordinal_parameters(theta, outcome_info, state_design, covariate_design)$beta)
  }

  theta[seq_len(length(term_names))]
}

.distal_table_components <- function(theta,
                                     vcov,
                                     outcome_info,
                                     state_design,
                                     covariate_design = NULL) {
  term_names <- .distal_predictor_term_names(state_design, covariate_design)
  if (outcome_info$type == "categorical") {
    return(list(
      coefficients = theta,
      vcov = vcov,
      term_names = term_names
    ))
  }

  if (outcome_info$type == "ordinal") {
    n_thresholds <- length(outcome_info$levels) - 1L
    beta_index <- seq.int(n_thresholds + 1L, length(theta))
    return(list(
      coefficients = theta[beta_index],
      vcov = vcov[beta_index, beta_index, drop = FALSE],
      term_names = .distal_predictor_term_names(state_design, covariate_design, ordinal = TRUE)
    ))
  }

  beta_index <- seq_len(length(term_names))
  list(
    coefficients = theta[beta_index],
    vcov = vcov[beta_index, beta_index, drop = FALSE],
    term_names = term_names
  )
}

.distal_posterior <- function(prior, misclassification, observed_state, like) {
  state_weight <- t(misclassification[, observed_state, drop = FALSE])
  state_weight <- state_weight * matrix(prior, nrow = nrow(state_weight), ncol = ncol(state_weight), byrow = TRUE)
  numerator <- state_weight * like
  denominator <- rowSums(numerator)
  numerator / pmax(denominator, .Machine$double.eps)
}

.distal_ml_fit_core <- function(outcome_info,
                                state_design,
                                covariate_design,
                                observed_state,
                                misclassification,
                                prior,
                                conf.level = 0.95,
                                maxit = 1000,
                                start = NULL) {
  misclassification <- misclassification / pmax(rowSums(misclassification), .Machine$double.eps)
  prior <- as.numeric(prior)
  prior <- pmax(prior, .Machine$double.eps)
  prior <- prior / sum(prior)

  if (is.null(start)) {
    start <- .distal_initial_theta(outcome_info, state_design, observed_state, covariate_design)
  }

  objective <- function(par) {
    like <- .distal_likelihood_matrix(par, outcome_info, state_design, covariate_design)$like
    state_weight <- t(misclassification[, observed_state, drop = FALSE])
    state_weight <- state_weight * matrix(prior, nrow = nrow(like), ncol = ncol(like), byrow = TRUE)
    mixture <- rowSums(state_weight * like)
    -sum(log(pmax(mixture, .Machine$double.eps)))
  }

  opt <- stats::optim(
    par = start,
    fn = objective,
    method = "BFGS",
    control = list(maxit = maxit, reltol = 1e-10)
  )

  hessian <- stats::optimHess(opt$par, objective)
  vcov <- .three_step_safe_inverse(hessian)
  like <- .distal_likelihood_matrix(opt$par, outcome_info, state_design, covariate_design)
  posterior <- .distal_posterior(prior, misclassification, observed_state, like$like)
  coefficient_object <- .distal_unpack_coefficients(opt$par, outcome_info, state_design, covariate_design)
  table_components <- .distal_table_components(opt$par, vcov, outcome_info, state_design, covariate_design)
  threshold_vcov <- NULL

  if (outcome_info$type == "ordinal") {
    n_thresholds <- length(outcome_info$levels) - 1L
    threshold_vcov <- .distal_ordinal_threshold_vcov(
      raw_thresholds = like$raw_thresholds,
      vcov_threshold_raw = vcov[seq_len(n_thresholds), seq_len(n_thresholds), drop = FALSE]
    )
  }

  list(
    coefficients = coefficient_object,
    raw_coefficients = opt$par,
    vcov = vcov,
    table = .distal_coef_table(
      coefficients = table_components$coefficients,
      vcov = table_components$vcov,
      term_names = table_components$term_names,
      outcome_info = outcome_info,
      conf.level = conf.level
    ),
    logLik = -opt$value,
    converged = opt$convergence == 0L,
    counts = opt$counts,
    posterior = posterior,
    state_fitted = like$state_fitted,
    state_summary = like$state_summary,
    thresholds = if (outcome_info$type == "ordinal") like$thresholds else NULL,
    threshold_vcov = threshold_vcov,
    threshold_table = if (outcome_info$type == "ordinal") {
      .distal_ordinal_threshold_table(like$thresholds, threshold_vcov, conf.level = conf.level)
    } else {
      NULL
    },
    residual_sd = if (outcome_info$type == "continuous") exp(opt$par[length(opt$par)]) else NULL
  )
}

.distal_bch_weights <- function(misclassification, observed_state) {
  misclassification <- misclassification / pmax(rowSums(misclassification), .Machine$double.eps)
  assignment <- diag(nrow(misclassification))[observed_state, , drop = FALSE]
  inverse_misclassification <- .three_step_safe_inverse(misclassification)

  list(
    weights = assignment %*% inverse_misclassification,
    inverse_misclassification = inverse_misclassification,
    misclassification = misclassification
  )
}

.distal_bch_fit_core <- function(outcome_info,
                                 state_design,
                                 covariate_design,
                                 bch_weights,
                                 conf.level = 0.95,
                                 maxit = 1000,
                                 start = NULL,
                                 observed_state = NULL) {
  if (is.null(start)) {
    if (is.null(observed_state)) {
      stop("observed_state is required when start is not supplied.", call. = FALSE)
    }
    start <- .distal_initial_theta(outcome_info, state_design, observed_state, covariate_design)
  }

  objective <- function(par) {
    like <- .distal_likelihood_matrix(par, outcome_info, state_design, covariate_design)$like
    -sum(bch_weights * log(pmax(like, .Machine$double.eps)))
  }

  opt <- stats::optim(
    par = start,
    fn = objective,
    method = "BFGS",
    control = list(maxit = maxit, reltol = 1e-10)
  )

  hessian <- stats::optimHess(opt$par, objective)
  vcov <- .three_step_safe_inverse(hessian)
  like <- .distal_likelihood_matrix(opt$par, outcome_info, state_design, covariate_design)
  coefficient_object <- .distal_unpack_coefficients(opt$par, outcome_info, state_design, covariate_design)
  table_components <- .distal_table_components(opt$par, vcov, outcome_info, state_design, covariate_design)
  threshold_vcov <- NULL

  if (outcome_info$type == "ordinal") {
    n_thresholds <- length(outcome_info$levels) - 1L
    threshold_vcov <- .distal_ordinal_threshold_vcov(
      raw_thresholds = like$raw_thresholds,
      vcov_threshold_raw = vcov[seq_len(n_thresholds), seq_len(n_thresholds), drop = FALSE]
    )
  }

  list(
    coefficients = coefficient_object,
    raw_coefficients = opt$par,
    vcov = vcov,
    table = .distal_coef_table(
      coefficients = table_components$coefficients,
      vcov = table_components$vcov,
      term_names = table_components$term_names,
      outcome_info = outcome_info,
      conf.level = conf.level
    ),
    pseudo_logLik = -opt$value,
    converged = opt$convergence == 0L,
    counts = opt$counts,
    state_fitted = like$state_fitted,
    state_summary = like$state_summary,
    thresholds = if (outcome_info$type == "ordinal") like$thresholds else NULL,
    threshold_vcov = threshold_vcov,
    threshold_table = if (outcome_info$type == "ordinal") {
      .distal_ordinal_threshold_table(like$thresholds, threshold_vcov, conf.level = conf.level)
    } else {
      NULL
    },
    residual_sd = if (outcome_info$type == "continuous") exp(opt$par[length(opt$par)]) else NULL
  )
}

.distal_fit_bundle <- function(outcome_info,
                               state_design,
                               covariate_design,
                               observed_state,
                               misclassification,
                               prior,
                               method,
                               conf.level,
                               maxit) {
  naive_fit <- .distal_ml_fit_core(
    outcome_info = outcome_info,
    state_design = state_design,
    covariate_design = covariate_design,
    observed_state = observed_state,
    misclassification = diag(nrow(misclassification)),
    prior = rep(1 / nrow(misclassification), nrow(misclassification)),
    conf.level = conf.level,
    maxit = maxit
  )

  out <- list(naive = naive_fit)

  if ("ML" %in% method) {
    out$ML <- .distal_ml_fit_core(
      outcome_info = outcome_info,
      state_design = state_design,
      covariate_design = covariate_design,
      observed_state = observed_state,
      misclassification = misclassification,
      prior = prior,
      start = naive_fit$raw_coefficients,
      conf.level = conf.level,
      maxit = maxit
    )
  }

  if ("BCH" %in% method) {
    bch <- .distal_bch_weights(misclassification, observed_state)
    out$BCH <- .distal_bch_fit_core(
      outcome_info = outcome_info,
      state_design = state_design,
      covariate_design = covariate_design,
      bch_weights = bch$weights,
      start = naive_fit$raw_coefficients,
      observed_state = observed_state,
      conf.level = conf.level,
      maxit = maxit
    )
    out$bch_weights <- bch
  }

  out
}
