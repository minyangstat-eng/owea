# ===========================================================================
# model.R -- formula-style model builder for standard GLMs.
#
# For the common models (linear / logistic / log-linear count) a user can
# describe the model with a compact formula instead of hand-writing
# info_vector / info_matrix:
#
#   link  : "identity" (linear), "logit" (logistic), "loglinear" (Poisson log
#           link), "multinomial" (baseline-category logit) or "cumulative"
#           (proportional-odds ordinal logit); the last two need ncat = J and
#           produce a k x k per-point information matrix (no info-vector form)
#   f     : main FACTOR effects            -- factor indices, e.g. c(1, 3)
#   x     : main CONTINUOUS linear effects -- continuous indices, e.g. c(1, 2, 3)
#   fx    : factor x continuous interactions -- two-digit codes c(11, 23) OR a
#           list of index pairs list(c(1,1), c(2,3))  (first = factor idx,
#           second = continuous idx)
#   ff    : factor x factor interactions -- two-digit codes c(12, 13) OR list of
#           pairs list(c(1,2), c(1,3)); the two factors must differ
#   xx    : continuous quadratic / continuous x continuous -- two-digit codes
#           c(11, 13) (equal digits = quadratic) OR list of pairs
#   intercept : include an intercept column (default TRUE)
#   coding    : factor contrast coding, "zero-sum" (default) or "baseline"
#
# Covariates are indexed WITHIN their own kind, in design_box order: factor
# index j = the j-th factor; continuous index j = the j-th continuous covariate.
# Factors come from design_box (a single integer L >= 2 = a factor with levels
# 1..L, see factors.R) or, for the raw-design helpers, from factor_levels /
# auto-detection.  The assembled model row q has the term order
# [intercept], f, x, fx, xx, ff, each factor term expanded to its contrast cols.
# ===========================================================================

# Normalize an fx / xx / ff specification (two-digit numeric codes OR a list of
# length-2 index pairs) into a list of integer c(i, j) pairs.
.norm_pairs <- function(codes, what) {
  if (is.null(codes)) return(list())
  if (is.list(codes)) {
    return(lapply(codes, function(p) {
      p <- as.integer(round(as.numeric(p)))
      if (length(p) != 2L || any(is.na(p)) || any(p < 1L))
        stop(sprintf("each '%s' term must be a pair of positive indices c(i, j).",
                     what), call. = FALSE)
      p
    }))
  }
  codes <- as.integer(round(as.numeric(codes)))
  lapply(codes, function(cd) {
    if (is.na(cd) || cd < 11L || cd > 99L || (cd %% 10L) == 0L)
      stop(sprintf(paste0("two-digit '%s' code %s is invalid; both digits must ",
                          "be 1..9, or pass an explicit list of index pairs ",
                          "c(i, j) for indices > 9."), what, cd), call. = FALSE)
    c(cd %/% 10L, cd %% 10L)
  })
}

# The length-(L-1) coded vector for a factor at integer `level` (1..L).
.contrast_code <- function(level, L, coding) {
  level <- as.integer(round(level))
  v <- numeric(max(L - 1L, 0L))
  if (L >= 2L) {
    if (level < L)      v[level] <- 1
    else if (coding == "zero-sum") v[] <- -1   # level == L
    # baseline: top level stays all zeros
  }
  v
}

# Covariate metadata for the builder from EITHER a design_box (preferred) OR a
# factor_levels vector (one entry per covariate; NA/0/1 = continuous, L>=2 =
# factor).  Returns the same shape as .parse_design_box (minus lo/hi when built
# from factor_levels).
.model_meta <- function(design_box = NULL, factor_levels = NULL, ncov = NULL) {
  if (!is.null(design_box)) return(.parse_design_box(design_box))
  if (!is.null(factor_levels)) {
    n <- length(factor_levels)
    m <- .factor_levels_to_meta(factor_levels, n)
    m$lo <- rep(NA_real_, n); m$hi <- rep(NA_real_, n)
    return(m)
  }
  stop("Supply 'design_box' or 'factor_levels' to describe the covariates.",
       call. = FALSE)
}

# Validate the model spec against the covariate metadata and precompute a "plan"
# (index maps, factor levels, total number of parameters k).
.build_model_terms <- function(meta, f, x, fx, xx, ff = NULL,
                               intercept = TRUE, coding = "zero-sum",
                               link = "identity", cov_names = NULL,
                               ncat = NULL) {
  coding <- match.arg(coding, c("zero-sum", "baseline"))
  link   <- match.arg(link, c("identity", "logit", "loglinear",
                              "multinomial", "cumulative"))
  multicat <- link %in% c("multinomial", "cumulative")
  J <- NA_integer_
  if (multicat) {
    if (is.null(ncat) || length(ncat) != 1L || is.na(ncat) ||
        ncat < 2 || abs(ncat - round(ncat)) > 1e-8)
      stop("'ncat' (number of response categories, an integer >= 2) is required ",
           "for the '", link, "' link.", call. = FALSE)
    J <- as.integer(round(ncat))
  }
  # The proportional-odds (cumulative) model carries J-1 threshold intercepts, so
  # its slope block uses the NON-intercept terms; any 'intercept' request is
  # folded into the thresholds.
  base_intercept <- if (link == "cumulative") FALSE else isTRUE(intercept)
  is_factor <- meta$is_factor
  fac_pos   <- which(is_factor)
  cont_pos  <- which(!is_factor)
  nfac  <- length(fac_pos)
  ncont <- length(cont_pos)
  nlev_fac <- as.integer(meta$nlevels[fac_pos])

  f  <- if (is.null(f)) integer(0) else as.integer(round(as.numeric(f)))
  x  <- if (is.null(x)) integer(0) else as.integer(round(as.numeric(x)))
  fx <- .norm_pairs(fx, "fx")
  xx <- .norm_pairs(xx, "xx")
  ff <- .norm_pairs(ff, "ff")

  if (length(f) && (any(f < 1L) || any(f > nfac)))
    stop("f indices must be within 1..", nfac, " (number of factors).",
         call. = FALSE)
  if (length(x) && (any(x < 1L) || any(x > ncont)))
    stop("x indices must be within 1..", ncont,
         " (number of continuous covariates).", call. = FALSE)
  for (p in fx) {
    if (p[1] > nfac)
      stop("fx factor index ", p[1], " exceeds number of factors (", nfac, ").",
           call. = FALSE)
    if (p[2] > ncont)
      stop("fx continuous index ", p[2],
           " exceeds number of continuous covariates (", ncont, ").",
           call. = FALSE)
  }
  for (p in xx)
    if (any(p > ncont))
      stop("xx indices must be within 1..", ncont,
           " (number of continuous covariates).", call. = FALSE)
  for (p in ff) {
    if (any(p > nfac))
      stop("ff indices must be within 1..", nfac, " (number of factors).",
           call. = FALSE)
    if (p[1] == p[2])
      stop("ff interaction must be between two DIFFERENT factors (got ",
           p[1], " x ", p[2], ").", call. = FALSE)
  }
  if ((length(f) || length(fx) || length(ff)) && nfac == 0L)
    stop("factor terms (f/fx/ff) requested but the design has no factor ",
         "covariates.", call. = FALSE)
  ref_fac <- unique(c(f,
                      if (length(fx)) vapply(fx, `[`, integer(1), 1L),
                      if (length(ff)) unlist(ff)))
  if (length(ref_fac) && any(nlev_fac[ref_fac] < 2L))
    stop("a referenced factor has fewer than 2 levels; supply 'factor_levels' ",
         "or a valid design_box.", call. = FALSE)

  # covariate display names (from design_box / colnames, else generic F#/X#)
  fac_nm  <- function(i) {
    nm <- if (!is.null(cov_names) && !is.na(cov_names[fac_pos[i]]) &&
              nzchar(cov_names[fac_pos[i]])) cov_names[fac_pos[i]]
          else paste0("F", i)
    nm
  }
  cont_nm <- function(j) {
    if (!is.null(cov_names) && !is.na(cov_names[cont_pos[j]]) &&
        nzchar(cov_names[cont_pos[j]])) cov_names[cont_pos[j]]
    else paste0("X", j)
  }
  fac_cols <- function(i) {                       # contrast-column labels
    L <- nlev_fac[i]; nm <- fac_nm(i)
    if (L == 2L) nm else paste0(nm, ".", seq_len(L - 1L))
  }

  # labels and length of the BASE model row f(x) (the per-point regressors)
  nm <- character(0)
  if (base_intercept) nm <- c(nm, "(Intercept)")
  for (i in f) nm <- c(nm, fac_cols(i))
  for (j in x) nm <- c(nm, cont_nm(j))
  for (p in fx) nm <- c(nm, paste0(fac_cols(p[1]), ":", cont_nm(p[2])))
  for (p in xx)
    nm <- c(nm, if (p[1] == p[2]) paste0(cont_nm(p[1]), "^2")
                else paste0(cont_nm(p[1]), ":", cont_nm(p[2])))
  for (p in ff)
    nm <- c(nm, as.vector(t(outer(fac_cols(p[1]), fac_cols(p[2]),
                                  paste, sep = ":"))))   # first factor slow

  base_k <- (if (base_intercept) 1L else 0L) +
       (if (length(f)) sum(nlev_fac[f] - 1L) else 0L) +
       length(x) +
       (if (length(fx))
          sum(vapply(fx, function(p) nlev_fac[p[1]] - 1L, integer(1))) else 0L) +
       length(xx) +
       (if (length(ff))
          sum(vapply(ff, function(p) (nlev_fac[p[1]] - 1L) *
                                     (nlev_fac[p[2]] - 1L), integer(1))) else 0L)
  if (length(nm) != base_k)
    stop("internal error: base coefficient label count (", length(nm),
         ") != base_k (", base_k, ").", call. = FALSE)

  # full parameter count k and coefficient labels:
  #  * single-response links: k = base_k
  #  * multinomial: k = (J-1) * base_k, theta stacked by category (baseline = J)
  #  * cumulative:  k = (J-1) thresholds + base_k shared slopes
  if (!multicat) {
    k <- base_k; coef_names <- nm
  } else if (link == "multinomial") {
    if (base_k < 1L)
      stop("the multinomial model has no per-category terms.", call. = FALSE)
    k <- (J - 1L) * base_k
    coef_names <- unlist(lapply(seq_len(J - 1L),
                                function(j) paste0("cat", j, ":", nm)))
  } else {                                            # cumulative
    k <- (J - 1L) + base_k
    coef_names <- c(sprintf("(threshold %d)", seq_len(J - 1L)), nm)
  }
  if (k < 1L)
    stop("the model has no terms (need an intercept or at least one effect).",
         call. = FALSE)

  list(intercept = base_intercept, f = f, x = x, fx = fx, xx = xx,
       ff = ff, fac_pos = fac_pos, cont_pos = cont_pos, nlev_fac = nlev_fac,
       coding = coding, link = link, ncat = J, base_k = as.integer(base_k),
       k = as.integer(k), coef_names = coef_names)
}

# Assemble the base model row f(x) (length plan$base_k) for a covariate point x.
.model_row <- function(x, plan) {
  x   <- as.numeric(x)
  q   <- numeric(plan$base_k)
  pos <- 0L
  if (plan$intercept) { pos <- pos + 1L; q[pos] <- 1 }
  for (j in plan$f) {
    cc <- .contrast_code(x[plan$fac_pos[j]], plan$nlev_fac[j], plan$coding)
    if (length(cc)) { q[(pos + 1L):(pos + length(cc))] <- cc; pos <- pos + length(cc) }
  }
  for (j in plan$x) { pos <- pos + 1L; q[pos] <- x[plan$cont_pos[j]] }
  for (p in plan$fx) {
    cc <- .contrast_code(x[plan$fac_pos[p[1]]], plan$nlev_fac[p[1]], plan$coding) *
          x[plan$cont_pos[p[2]]]
    if (length(cc)) { q[(pos + 1L):(pos + length(cc))] <- cc; pos <- pos + length(cc) }
  }
  for (p in plan$xx) {
    pos <- pos + 1L
    q[pos] <- x[plan$cont_pos[p[1]]] * x[plan$cont_pos[p[2]]]
  }
  for (p in plan$ff) {
    ai <- .contrast_code(x[plan$fac_pos[p[1]]], plan$nlev_fac[p[1]], plan$coding)
    bj <- .contrast_code(x[plan$fac_pos[p[2]]], plan$nlev_fac[p[2]], plan$coding)
    cc <- kronecker(ai, bj)          # first factor slow (matches coef_names)
    if (length(cc)) { q[(pos + 1L):(pos + length(cc))] <- cc; pos <- pos + length(cc) }
  }
  q
}

# Build the per-point info_vector function for a plan.  identity -> function(x)
# (theta unused); logit / loglinear -> function(x, theta).
.make_model_info_vector <- function(plan) {
  f <- if (plan$link == "identity") {
    function(x) .model_row(x, plan)
  } else if (plan$link == "logit") {
    function(x, theta) {
      q   <- .model_row(x, plan)
      eta <- min(max(sum(q * theta), -500), 500)   # guard exp() overflow
      (exp(eta / 2) / (1 + exp(eta))) * q
    }
  } else {                                          # loglinear (Poisson log link)
    function(x, theta) {
      q   <- .model_row(x, plan)
      eta <- min(max(sum(q * theta), -500), 500)
      exp(eta / 2) * q
    }
  }
  attr(f, "coef_names") <- plan$coef_names
  attr(f, "link")       <- plan$link
  attr(f, "coding")     <- plan$coding
  f
}

# Per-observation k x k information at a BASE model row frow (length p) for the
# baseline-category multinomial: M = W (x) frow frow', W = diag(pi) - pi pi'
# over categories 1..J-1, pi = softmax, theta = (beta_1,...,beta_{J-1}) stacked.
.info_multinomial <- function(frow, theta, J) {
  p   <- length(frow)
  B   <- matrix(as.numeric(theta), nrow = p, ncol = J - 1L)  # col j = beta_j
  eta <- pmin(pmax(as.numeric(crossprod(frow, B)), -500), 500)
  ex  <- exp(eta)
  pri <- ex / (1 + sum(ex))                          # pi_1..pi_{J-1}
  W   <- diag(pri, J - 1L) - tcrossprod(pri)
  kronecker(W, tcrossprod(frow))
}

# Per-observation k x k information at a base (non-intercept) row frow (length q)
# for the proportional-odds cumulative model: sum_c (1/pi_c)(g_c-g_{c-1})(.)',
# g_j = lambda_j (e_j ; frow), gamma_j = plogis(alpha_j + frow'beta),
# theta = (alpha_1,...,alpha_{J-1}, beta) with alpha strictly increasing.
.info_cumulative <- function(frow, theta, J) {
  q     <- length(frow); kk <- (J - 1L) + q
  theta <- as.numeric(theta)
  alpha <- theta[seq_len(J - 1L)]
  beta  <- if (q > 0L) theta[(J - 1L) + seq_len(q)] else numeric(0)
  fb    <- if (q > 0L) sum(frow * beta) else 0
  eta   <- pmin(pmax(alpha + fb, -500), 500)         # length J-1
  gamma <- plogis(eta)                               # cumulative probs
  lam   <- gamma * (1 - gamma)                       # densities lambda_j
  pri   <- diff(c(0, gamma, 1))                      # pi_1..pi_J (>0 iff alpha up)
  G <- matrix(0, kk, J - 1L)                         # g_j as columns
  for (j in seq_len(J - 1L)) {
    G[j, j] <- lam[j]
    if (q > 0L) G[(J - 1L) + seq_len(q), j] <- lam[j] * frow
  }
  Gext <- cbind(0, G, 0)                             # g_0, g_1..g_{J-1}, g_J
  M <- matrix(0, kk, kk)
  for (cc in seq_len(J)) {
    d  <- Gext[, cc + 1L] - Gext[, cc]               # g_c - g_{c-1}
    pc <- pri[cc]
    if (pc > 1e-12) M <- M + tcrossprod(d) / pc
  }
  M
}

# Build the per-point k x k information matrix function for a multi-category
# model.  Returns function(x, theta) -> k x k matrix (rank J-1 per point).
.make_model_info_matrix_multicat <- function(plan) {
  J <- plan$ncat
  im <- if (plan$link == "multinomial")
          function(x, theta) .info_multinomial(.model_row(x, plan), theta, J)
        else
          function(x, theta) .info_cumulative(.model_row(x, plan), theta, J)
  attr(im, "coef_names") <- plan$coef_names
  attr(im, "link")       <- plan$link
  attr(im, "coding")     <- plan$coding
  attr(im, "ncat")       <- J
  im
}

# Covariate metadata for the RAW-DESIGN helpers (infor_matrix /
# design_information): use factor_levels if given, else auto-detect (a column
# whose values are all positive integers is a factor with L = column maximum).
.factor_meta_from_design <- function(X, factor_levels = NULL, quiet = FALSE) {
  X  <- as.matrix(X)
  nc <- ncol(X)
  if (!is.null(factor_levels)) return(.factor_levels_to_meta(factor_levels, nc))
  is_factor <- logical(nc); nlevels <- rep(NA_integer_, nc)
  for (j in seq_len(nc)) {
    v <- X[, j]
    if (all(is.finite(v)) && all(v > 0) && all(abs(v - round(v)) < 1e-8)) {
      is_factor[j] <- TRUE
      nlevels[j]   <- as.integer(max(round(v)))
    }
  }
  if (any(is_factor) && !quiet) {
    fj <- which(is_factor)
    message(sprintf(paste0("model+link: auto-detected factor column(s) %s with ",
                           "levels %s (L = column max). Pass 'factor_levels' to ",
                           "override."),
                    paste(fj, collapse = ", "),
                    paste(nlevels[fj], collapse = ", ")))
  }
  list(is_factor = is_factor, nlevels = nlevels)
}

# Draw a random theta for a local design when none is supplied (with a warning).
# For the cumulative model the threshold block is sorted so that pi_c > 0 at
# every x.
.random_theta <- function(plan) {
  if (plan$link == "cumulative") {
    J <- plan$ncat
    c(sort(stats::rnorm(J - 1L)), stats::rnorm(plan$base_k))
  } else {
    stats::rnorm(plan$k)
  }
}

# Resolve a model spec into (info_vector, info_matrix, theta) for the solvers.
# When `link` is NULL the caller supplied its own model -- return untouched.
# For single-response links (logit/loglinear) an info_vector is built; the
# multi-category links (multinomial/cumulative) build a k x k info_matrix.  A
# missing theta for a local model is drawn ~ N(0,1) (with a warning).
.resolve_model_spec <- function(link, f, x, fx, xx, intercept, coding,
                                design_box, candidate_set, factor_levels,
                                info_vector, info_matrix, theta, ff = NULL,
                                ncat = NULL) {
  if (is.null(link))
    return(list(info_vector = info_vector, info_matrix = info_matrix,
                theta = theta, factor_levels = factor_levels,
                coef_names = NULL, link = NULL, spec_given = FALSE))
  if (!is.null(info_vector) || !is.null(info_matrix))
    stop("Supply a model spec ('link' with f/x/fx/ff/xx) OR an ",
         "'info_vector'/'info_matrix', not both.", call. = FALSE)

  det_fl <- factor_levels
  cov_names <- names(design_box)
  meta <- if (!is.null(design_box)) {
            .parse_design_box(design_box)
          } else if (!is.null(candidate_set)) {
            if (is.null(cov_names))
              cov_names <- colnames(as.matrix(candidate_set))
            m <- .factor_meta_from_design(candidate_set, factor_levels)
            if (is.null(factor_levels))
              det_fl <- ifelse(m$is_factor, m$nlevels, NA_integer_)
            m
          } else {
            stop("Supply 'design_box' (or 'candidate_set') with a model spec.",
                 call. = FALSE)
          }
  plan     <- .build_model_terms(meta, f, x, fx, xx, ff, intercept, coding, link,
                                 cov_names, ncat = ncat)
  multicat <- plan$link %in% c("multinomial", "cumulative")

  # local model (everything but identity) needs theta; draw it if missing.
  if (plan$link != "identity" && is.null(theta)) {
    theta <- .random_theta(plan)
    warning(sprintf(paste0("no 'theta' supplied for the %s model; each of the ",
                           "%d parameter(s) was drawn from N(0,1)%s: %s. The ",
                           "design is only locally optimal at this theta -- ",
                           "supply 'theta' for a specific parameter value."),
                    plan$link, plan$k,
                    if (plan$link == "cumulative") " (thresholds sorted)" else "",
                    paste(sprintf("%.4f", theta), collapse = ", ")),
            call. = FALSE)
  }
  if (!is.null(theta)) {
    theta <- as.numeric(theta)
    if (length(theta) != plan$k)
      stop(sprintf("the %s model has %d parameters but 'theta' has length %d.",
                   plan$link, plan$k, length(theta)), call. = FALSE)
    if (plan$link == "cumulative") {
      a <- theta[seq_len(plan$ncat - 1L)]
      if (length(a) > 1L && any(diff(a) <= 0))
        stop("cumulative model: the first ", plan$ncat - 1L, " elements of ",
             "'theta' are the thresholds and must be strictly increasing.",
             call. = FALSE)
    }
  }

  if (multicat) {
    im <- .make_model_info_matrix_multicat(plan)
    return(list(info_vector = NULL, info_matrix = im, theta = theta,
                factor_levels = det_fl, coef_names = plan$coef_names,
                link = plan$link, spec_given = TRUE))
  }
  iv <- .make_model_info_vector(plan)
  list(info_vector = iv, info_matrix = NULL, theta = theta,
       factor_levels = det_fl, coef_names = plan$coef_names, link = plan$link,
       spec_given = TRUE)
}

#' Build an information-vector model from a formula-style spec.
#'
#' Constructs an \code{info_vector(x, theta)} (or \code{info_vector(x)} for the
#' identity link) for a standard linear / logistic / log-linear model, generating
#' factor contrast coding automatically, so it can be passed to
#' \code{\link{optimal_design}}, \code{\link{exact_design}} or
#' \code{\link{infor_matrix}}. (Those solvers also accept the same
#' \code{link}/\code{f}/\code{x}/\code{fx}/\code{ff}/\code{xx} arguments
#' directly.) The returned function carries the coefficient labels in its
#' \code{"coef_names"} attribute (see \code{\link{model_summary}}).
#'
#' @param design_box list describing the covariates, as in
#'   \code{\link{optimal_design}} (a \code{c(lo, hi)} pair is continuous; a
#'   single integer \code{L >= 2} is a factor with levels \code{1..L}). Supply
#'   this OR \code{factor_levels}. A named list labels the coefficients.
#' @param factor_levels integer vector with one entry per covariate marking
#'   factor columns (\code{L >= 2}) versus continuous (\code{NA}/\code{0}/
#'   \code{1}); an alternative to \code{design_box}.
#' @param link \code{"identity"} (linear, default), \code{"logit"} (logistic),
#'   \code{"loglinear"} (Poisson log-link count model), \code{"multinomial"}
#'   (baseline-category logit) or \code{"cumulative"} (proportional-odds ordinal
#'   logit). The last two need \code{ncat} and have no information-vector form
#'   (use \code{\link{model_info_matrix}} / \code{\link{optimal_design}}).
#' @param ncat number of response categories \eqn{J \ge 2} for the
#'   \code{"multinomial"} / \code{"cumulative"} links (ignored otherwise).
#' @param f main factor-effect indices (which factors), e.g. \code{c(1, 3)}.
#' @param x main continuous linear-effect indices, e.g. \code{c(1, 2, 3)}.
#' @param fx factor \eqn{\times} continuous interactions: two-digit codes
#'   \code{c(11, 23)} (first digit = factor index, second = continuous index) OR
#'   a list of index pairs \code{list(c(1, 1), c(2, 3))}.
#' @param xx continuous quadratic / continuous \eqn{\times} continuous terms:
#'   two-digit codes \code{c(11, 13)} (equal digits = quadratic) OR a list of
#'   index pairs.
#' @param ff factor \eqn{\times} factor interactions: two-digit codes
#'   \code{c(12, 13)} (both digits are factor indices) OR a list of index pairs
#'   \code{list(c(1, 2), c(1, 3))}. Each such term expands to
#'   \eqn{(L_i-1)(L_j-1)} contrast columns; the two factors must differ.
#' @param intercept include an intercept column (default \code{TRUE}).
#' @param coding factor contrast coding: \code{"zero-sum"} (default) or
#'   \code{"baseline"}.
#' @return a function \code{info_vector(x, theta)} (or \code{info_vector(x)} for
#'   the identity link) returning the length-\eqn{k} information vector, with a
#'   \code{"coef_names"} attribute.
#' @seealso \code{\link{model_info_matrix}}, \code{\link{model_summary}},
#'   \code{\link{optimal_design}}.
#' @export
model_info_vector <- function(design_box = NULL, factor_levels = NULL,
                              link = "identity", f = NULL, x = NULL,
                              fx = NULL, xx = NULL, ff = NULL, intercept = TRUE,
                              coding = "zero-sum", ncat = NULL) {
  meta <- .model_meta(design_box, factor_levels)
  cov_names <- if (!is.null(design_box)) names(design_box) else names(factor_levels)
  plan <- .build_model_terms(meta, f, x, fx, xx, ff, intercept, coding, link,
                             cov_names, ncat = ncat)
  if (plan$link %in% c("multinomial", "cumulative"))
    stop("the '", plan$link, "' model has no information-vector form; use ",
         "model_info_matrix() or optimal_design(link = \"", plan$link, "\").",
         call. = FALSE)
  .make_model_info_vector(plan)
}

#' Build an information-matrix model from a formula-style spec.
#'
#' The information-matrix counterpart of \code{\link{model_info_vector}}: returns
#' a function giving the per-point \eqn{k \times k} information matrix -- the
#' rank-1 \eqn{f f'} (identity) or GLM-weighted \eqn{f f'} (logit / loglinear),
#' or the full rank-\eqn{(J-1)} matrix of a \code{"multinomial"} /
#' \code{"cumulative"} model.
#'
#' @inheritParams model_info_vector
#' @return a function \code{info_matrix(x, theta)} (or \code{info_matrix(x)} for
#'   the identity link) returning the \eqn{k \times k} per-point information
#'   matrix.
#' @seealso \code{\link{model_info_vector}}, \code{\link{design_information}}.
#' @export
model_info_matrix <- function(design_box = NULL, factor_levels = NULL,
                              link = "identity", f = NULL, x = NULL,
                              fx = NULL, xx = NULL, ff = NULL, intercept = TRUE,
                              coding = "zero-sum", ncat = NULL) {
  meta <- .model_meta(design_box, factor_levels)
  cov_names <- if (!is.null(design_box)) names(design_box) else names(factor_levels)
  plan <- .build_model_terms(meta, f, x, fx, xx, ff, intercept, coding, link,
                             cov_names, ncat = ncat)
  if (plan$link %in% c("multinomial", "cumulative"))
    return(.make_model_info_matrix_multicat(plan))
  iv <- .make_model_info_vector(plan)
  im <- if (.info_vec_needs_theta(iv))
          function(x, theta) tcrossprod(as.numeric(iv(x, theta)))
        else
          function(x) tcrossprod(as.numeric(iv(x)))
  attr(im, "coef_names") <- attr(iv, "coef_names")
  attr(im, "link")       <- attr(iv, "link")
  attr(im, "coding")     <- attr(iv, "coding")
  im
}

#' Summarise a formula-style model: label each coefficient by its term.
#'
#' Prints the model's coefficient/term labels (intercept, factor contrasts,
#' continuous effects and interactions), and, when available, the parameter
#' value \code{theta} used and the design criterion. Works on the result of
#' \code{\link{optimal_design}} / \code{\link{exact_design}} built from a model
#' spec, or on a function returned by \code{\link{model_info_vector}} /
#' \code{\link{model_info_matrix}}.
#'
#' @param object a design result (from \code{\link{optimal_design}} or
#'   \code{\link{exact_design}} with a model spec), or a model function from
#'   \code{\link{model_info_vector}} / \code{\link{model_info_matrix}}.
#' @param ... ignored.
#' @return invisibly, a \code{data.frame} with one row per coefficient
#'   (\code{term}, and \code{theta} when known).
#' @seealso \code{\link{model_info_vector}}, \code{\link{optimal_design}}.
#' @export
model_summary <- function(object, ...) {
  cn <- NULL; th <- NULL; link <- NULL; crit <- NULL; pp <- NULL
  if (is.function(object)) {
    cn <- attr(object, "coef_names"); link <- attr(object, "link")
  } else if (is.list(object)) {
    cn <- object$coef_names; th <- object$theta; link <- object$link
    crit <- object$criterion; pp <- object$p
  }
  if (is.null(cn))
    stop("no model-spec coefficient labels found: this object was not built ",
         "from a formula-style model spec ('link' + terms).", call. = FALSE)

  has_theta <- !is.null(th) && length(th) == length(cn)
  crit_lab  <- if (!is.null(pp))
                 switch(as.character(pp), "0" = "D", "1" = "A",
                        paste0("Phi_", pp)) else NULL

  cat(sprintf("Model: %s link, %d coefficient(s)\n",
              if (is.null(link)) "identity" else link, length(cn)))
  if (!is.null(crit))
    cat(sprintf("Criterion%s: %.6f\n",
                if (is.null(crit_lab)) "" else sprintf(" (%s-optimality)", crit_lab),
                crit))
  wid <- max(nchar(cn), nchar("term"))
  if (has_theta) {
    cat(sprintf("  %-3s  %-*s  %s\n", "#", wid, "term", "theta"))
    for (i in seq_along(cn))
      cat(sprintf("  %-3d  %-*s  % .4f\n", i, wid, cn[i], th[i]))
  } else {
    cat(sprintf("  %-3s  %s\n", "#", "term"))
    for (i in seq_along(cn))
      cat(sprintf("  %-3d  %s\n", i, cn[i]))
  }
  invisible(data.frame(term = cn,
                       theta = if (has_theta) th else NA_real_,
                       stringsAsFactors = FALSE))
}
