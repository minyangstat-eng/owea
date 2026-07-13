# ===========================================================================
# model_fit.R -- fit a formula-style model to an OBSERVED data set (the
# "analyse" half of simulate_design(), with no simulation and no true theta).
# The maximum-likelihood fitters themselves live in model_sim.R (.fit_model).
# ===========================================================================

# Split a user data set into a numeric covariate matrix and the response.
# Returns list(X, y, cov_names, fac_cols, fac_levels).  Factor / character /
# logical covariate columns are recoded to integer levels 1..L (the coding the
# design side uses); numeric columns are passed through untouched.
.split_fit_data <- function(data, response) {
  if (is.matrix(data)) data <- as.data.frame(data)
  if (!is.data.frame(data))
    stop("'data' must be a data.frame or a matrix.", call. = FALSE)
  if (!nrow(data)) stop("'data' has no rows.", call. = FALSE)

  if (is.null(response)) {
    # default: the column named "y" if there is one (what simulate_design
    # produces), otherwise the LAST column.
    j <- match("y", names(data))
    if (is.na(j)) j <- ncol(data)
  } else if (is.numeric(response) && length(response) == 1L) {
    j <- as.integer(response)
    if (j < 1L || j > ncol(data))
      stop("'response' column index is out of range.", call. = FALSE)
  } else {
    if (!is.character(response) || length(response) != 1L)
      stop("'response' must be a single column name or index.", call. = FALSE)
    j <- match(response, names(data))
    if (is.na(j))
      stop("response column '", response, "' not found in 'data' (columns: ",
           paste(names(data), collapse = ", "), ").", call. = FALSE)
  }

  y   <- data[[j]]
  cov <- data[, -j, drop = FALSE]
  if (!ncol(cov))
    stop("'data' must have at least one covariate column besides the response.",
         call. = FALSE)

  d   <- ncol(cov)
  X   <- matrix(0, nrow(cov), d)
  fac <- logical(d)
  lev <- vector("list", d)
  for (j2 in seq_len(d)) {
    v <- cov[[j2]]
    if (is.factor(v) || is.character(v) || is.logical(v)) {
      fv <- factor(v)
      fac[j2]   <- TRUE
      lev[[j2]] <- levels(fv)
      X[, j2]   <- as.integer(fv)
    } else if (is.numeric(v)) {
      X[, j2] <- as.numeric(v)
    } else {
      stop("covariate column '", names(cov)[j2],
           "' has unsupported type '", class(v)[1], "'.", call. = FALSE)
    }
  }
  colnames(X) <- names(cov)
  list(X = X, y = y, resp_name = names(data)[j], cov_names = names(cov),
       fac_cols = fac, fac_levels = lev)
}

# Covariate metadata for a fitted data set: design_box / factor_levels win; then
# columns that arrived as R factors; then the usual auto-detection.
.fit_meta <- function(sp, design_box, factor_levels) {
  if (!is.null(design_box)) return(.parse_design_box(design_box))
  if (!is.null(factor_levels))
    return(.factor_meta_from_design(sp$X, factor_levels))
  meta <- .factor_meta_from_design(sp$X, NULL, quiet = any(sp$fac_cols))
  for (j in which(sp$fac_cols)) {
    meta$is_factor[j] <- TRUE
    meta$nlevels[j]   <- length(sp$fac_levels[[j]])
  }
  meta
}

# Coerce / validate the response for a given link; returns the numeric y the
# fitters expect (0/1 for logit, counts for loglinear, 1..J for multi-category).
.check_response <- function(y, plan, nm) {
  link <- plan$link
  if (link == "identity") {
    if (!is.numeric(y))
      stop("the '", nm, "' column must be numeric for the identity link.",
           call. = FALSE)
    return(as.numeric(y))
  }
  if (link == "logit") {
    if (is.logical(y)) return(as.numeric(y))
    if (is.factor(y) || is.character(y)) {
      fy <- factor(y)
      if (nlevels(fy) != 2L)
        stop("the '", nm, "' column must have exactly 2 levels for the logit ",
             "link; got ", nlevels(fy), ".", call. = FALSE)
      message("logit: treating '", levels(fy)[2], "' as the success (y = 1).")
      return(as.numeric(fy) - 1)
    }
    y <- as.numeric(y)
    if (!all(y %in% c(0, 1)))
      stop("the '", nm, "' column must be 0/1 (or logical, or a 2-level ",
           "factor) for the logit link.", call. = FALSE)
    return(y)
  }
  if (link == "loglinear") {
    y <- as.numeric(y)
    if (any(y < 0) || any(abs(y - round(y)) > 1e-8))
      stop("the '", nm, "' column must be nonnegative integer counts for the ",
           "loglinear link.", call. = FALSE)
    return(round(y))
  }
  # multinomial / cumulative: integer category codes 1..J
  J <- plan$ncat
  if (is.factor(y) || is.character(y)) {
    fy <- factor(y)
    if (nlevels(fy) > J)
      stop("the '", nm, "' column has ", nlevels(fy), " levels but ncat = ", J,
           ".", call. = FALSE)
    if (is.character(y))
      message(plan$link, ": category codes taken from the sorted levels (",
              paste(levels(fy), collapse = " < "), "). Supply '", nm,
              "' as a factor with the levels in the intended order, or as ",
              "integer codes 1..", J, ", to control the ordering.")
    return(as.integer(fy))
  }
  y <- as.integer(round(as.numeric(y)))
  if (any(y < 1L) || any(y > J))
    stop("the '", nm, "' column must be integer category codes in 1..", J,
         " (ncat = ", J, ") for the '", plan$link, "' link.", call. = FALSE)
  y
}

#' Fit a formula-style model to an observed data set (maximum likelihood).
#'
#' The estimation half of \code{\link{simulate_design}}, on its own: given a data
#' set of runs (covariate columns plus a response) and the same formula-style
#' model specification used elsewhere in the package (\code{link} plus
#' \code{f}/\code{x}/\code{fx}/\code{ff}/\code{xx}), obtain the maximum-likelihood
#' estimate of \eqn{\theta} together with its standard errors and covariance
#' matrix. No design, no true \code{theta} and no simulation are involved -- the
#' rows of \code{data} are the observed runs, replicates included.
#'
#' The returned \code{theta_hat} uses exactly the parameterisation and coefficient
#' ordering of \code{\link{optimal_design}} and \code{\link{simulate_design}}
#' (stacked by category for \code{"multinomial"}, thresholds-then-slopes for
#' \code{"cumulative"}), so it can be fed straight back in as a locally-optimal
#' design's \code{theta}.
#'
#' @param data a data frame (or matrix) with one row per run: the covariate
#'   columns, in the same order as the design's columns, followed by the response
#'   column. The \code{data} element returned by \code{\link{simulate_design}} has
#'   exactly this layout. Covariate columns may be numeric, or \code{factor} /
#'   \code{character} / \code{logical} (recoded to integer levels 1..L).
#' @param link,f,x,fx,ff,xx,intercept,coding,ncat the model, exactly as in
#'   \code{\link{optimal_design}}.
#' @param design_box,factor_levels how to interpret the covariate columns (which
#'   are factors and their levels). Supply one; otherwise columns that arrive as
#'   \code{factor}/\code{character}/\code{logical} are treated as factors and the
#'   remaining columns are auto-detected (all-positive-integer columns), which may
#'   misclassify integer-valued continuous covariates -- and note that
#'   auto-detection can only see the levels that actually occur in \code{data}.
#' @param response which column of \code{data} holds the response. By default the
#'   column named \code{"y"} if there is one, else the LAST column; every other
#'   column is a covariate. Give a name or a column index to override (e.g. when
#'   the response sits in the middle of the data frame).
#' @return a list with \code{theta_hat} (named), \code{se}, \code{vcov},
#'   \code{loglik}, \code{converged}, \code{sigma_hat} (identity link only),
#'   \code{coef_names}, \code{link}, \code{ncat} and \code{N} (number of runs).
#' @seealso \code{\link{simulate_design}}, \code{\link{optimal_design}},
#'   \code{\link{design_information}}.
#' @examples
#' \dontrun{
#' sim <- simulate_design(support = cbind(c(-1, 0, 1)), counts = c(10, 10, 10),
#'                        link = "logit", x = 1, xx = list(c(1, 1)),
#'                        theta = c(0.5, 1, -0.5), seed = 1, fit = FALSE)
#' fit_design(sim$data, link = "logit", x = 1, xx = list(c(1, 1)))
#' }
#' @export
fit_design <- function(data, link = NULL,
                       f = NULL, x = NULL, fx = NULL, xx = NULL, ff = NULL,
                       intercept = TRUE, coding = "zero-sum", ncat = NULL,
                       design_box = NULL, factor_levels = NULL,
                       response = NULL) {
  if (is.null(link)) stop("'link' is required.", call. = FALSE)
  sp <- .split_fit_data(data, response)

  meta      <- .fit_meta(sp, design_box, factor_levels)
  cov_names <- if (!is.null(design_box)) names(design_box) else sp$cov_names
  if (length(meta$is_factor) != ncol(sp$X))
    stop(sprintf(paste0("the model spec describes %d covariate(s) but 'data' ",
                        "has %d covariate column(s)."),
                 length(meta$is_factor), ncol(sp$X)), call. = FALSE)

  plan <- .build_model_terms(meta, f, x, fx, xx, ff, intercept, coding, link,
                             cov_names, ncat = ncat)

  y <- .check_response(sp$y, plan, sp$resp_name)
  N <- nrow(sp$X)
  if (N < plan$k)
    warning(sprintf(paste0("only %d run(s) for %d parameter(s); the fit may ",
                           "be unidentifiable."), N, plan$k), call. = FALSE)

  # run-level base model matrix f(x), one row per observation
  Xrun <- matrix(0, N, plan$base_k)
  for (i in seq_len(N)) Xrun[i, ] <- .model_row(sp$X[i, ], plan)

  fitr <- .fit_model(plan, Xrun, y, sigma = 1)

  vcov <- as.matrix(fitr$vcov)
  dimnames(vcov) <- list(plan$coef_names, plan$coef_names)
  theta_hat <- stats::setNames(fitr$theta_hat, plan$coef_names)
  se        <- stats::setNames(sqrt(diag(vcov)), plan$coef_names)

  out <- list(theta_hat = theta_hat, se = se, vcov = vcov,
              loglik = fitr$loglik, converged = fitr$converged,
              coef_names = plan$coef_names, link = plan$link,
              ncat = plan$ncat, N = N)
  if (!is.null(fitr$sigma_hat)) out$sigma_hat <- fitr$sigma_hat
  if (!isTRUE(fitr$converged))
    warning("the maximum-likelihood fit did not converge.", call. = FALSE)
  out
}
