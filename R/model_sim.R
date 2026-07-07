# ===========================================================================
# model_sim.R -- simulate responses for a design under a formula-style model and
# fit the model back (maximum likelihood), optionally over many Monte Carlo
# replications to compare the empirical covariance of the estimates with the
# design-based large-sample covariance M(xi, theta)^{-1} / N.
# ===========================================================================

# ---- vectorised categorical sampling --------------------------------------
# One draw per row of the N x J probability matrix P (rows sum to 1).
.rcat <- function(P) {
  J  <- ncol(P)
  cs <- P %*% upper.tri(matrix(0, J, J), diag = TRUE)   # cumulative sums by row
  1L + rowSums(cs < runif(nrow(P)))
}

# ---- response simulation for a run-level model matrix ---------------------
# Xrun: N x base_k model matrix (base row f(x) per run; NON-intercept for
# cumulative).  Returns the length-N response (integer category for multi-cat).
.simulate_y <- function(plan, Xrun, theta, sigma) {
  N <- nrow(Xrun); link <- plan$link
  if (link == "identity") {
    stats::rnorm(N, as.numeric(Xrun %*% theta), sigma)
  } else if (link == "logit") {
    stats::rbinom(N, 1L, stats::plogis(as.numeric(Xrun %*% theta)))
  } else if (link == "loglinear") {
    stats::rpois(N, exp(pmin(as.numeric(Xrun %*% theta), 700)))
  } else if (link == "multinomial") {
    J <- plan$ncat; p <- plan$base_k
    B   <- matrix(theta, p, J - 1L)
    Eta <- pmin(pmax(Xrun %*% B, -500), 500)
    EX  <- exp(Eta); denom <- 1 + rowSums(EX)
    .rcat(cbind(EX / denom, 1 / denom))               # baseline = category J
  } else {                                            # cumulative
    J <- plan$ncat; q <- plan$base_k
    alpha <- theta[seq_len(J - 1L)]
    beta  <- if (q > 0L) theta[(J - 1L) + seq_len(q)] else numeric(0)
    lin   <- if (q > 0L) as.numeric(Xrun %*% beta) else rep(0, N)
    G  <- stats::plogis(outer(lin, alpha, `+`))       # N x (J-1) cumulative
    Gf <- cbind(0, G, 1)
    .rcat(Gf[, 2:(J + 1L), drop = FALSE] - Gf[, 1:J, drop = FALSE])
  }
}

# ---- maximum-likelihood fitters -------------------------------------------
# Each returns list(theta_hat, vcov, loglik, converged) in the SAME
# parameterisation the design/spec uses.

.fit_gaussian <- function(X, y) {
  XtX <- crossprod(X); b <- solve(XtX, crossprod(X, y))
  N <- nrow(X); p <- ncol(X)
  s2 <- sum((y - X %*% b)^2) / max(N - p, 1L)
  list(theta_hat = as.numeric(b), vcov = s2 * solve(XtX),
       loglik = sum(stats::dnorm(y, as.numeric(X %*% b), sqrt(s2), log = TRUE)),
       converged = TRUE, sigma_hat = sqrt(s2))
}

.fit_glm <- function(X, y, link) {
  fam <- if (link == "logit") stats::binomial() else stats::poisson()
  fit <- suppressWarnings(stats::glm.fit(X, y, family = fam))
  b   <- as.numeric(fit$coefficients)
  mu  <- fam$linkinv(as.numeric(X %*% b))
  w   <- if (link == "logit") mu * (1 - mu) else mu
  info <- crossprod(X * w, X)
  ll  <- if (link == "logit") sum(stats::dbinom(y, 1L, mu, log = TRUE))
         else sum(stats::dpois(y, mu, log = TRUE))
  list(theta_hat = b, vcov = solve(info), loglik = ll,
       converged = isTRUE(fit$converged))
}

# Baseline-category multinomial via Newton-Raphson (canonical -> observed =
# expected information).  b is (beta_1,...,beta_{J-1}) stacked.
.fit_multinomial <- function(X, y, J, maxit = 100L, tol = 1e-8) {
  N <- nrow(X); p <- ncol(X); k <- p * (J - 1L)
  probs <- function(b) {
    Eta <- pmin(pmax(X %*% matrix(b, p, J - 1L), -500), 500)
    EX  <- exp(Eta); denom <- 1 + rowSums(EX)
    cbind(EX / denom, 1 / denom)                      # N x J
  }
  Y  <- outer(y, seq_len(J - 1L), `==`) + 0           # N x (J-1) indicators
  b  <- numeric(k); conv <- FALSE
  H  <- diag(k)
  for (it in seq_len(maxit)) {
    P    <- probs(b); Pnb <- P[, seq_len(J - 1L), drop = FALSE]
    grad <- as.vector(crossprod(X, Y - Pnb))
    H    <- matrix(0, k, k)
    for (j in seq_len(J - 1L)) for (l in j:(J - 1L)) {
      w   <- if (j == l) Pnb[, j] * (1 - Pnb[, j]) else -Pnb[, j] * Pnb[, l]
      blk <- crossprod(X * w, X)
      rj  <- (j - 1L) * p + seq_len(p); rl <- (l - 1L) * p + seq_len(p)
      H[rj, rl] <- blk; if (j != l) H[rl, rj] <- t(blk)
    }
    step <- tryCatch(solve(H, grad), error = function(e) NULL)
    if (is.null(step)) break
    b <- b + step
    if (max(abs(grad)) < tol) { conv <- TRUE; break }
  }
  P  <- probs(b); ll <- sum(log(pmax(P[cbind(seq_len(N), y)], 1e-12)))
  list(theta_hat = b,
       vcov = tryCatch(solve(H), error = function(e) matrix(NA_real_, k, k)),
       loglik = ll, converged = conv)
}

# Proportional-odds cumulative model via optim with an ordered-threshold
# reparameterisation (alpha_1, log-increments, beta); vcov from the Fisher
# information at the MLE in the ORIGINAL (alpha, beta) parameterisation.
.fit_cumulative <- function(X, y, J, maxit = 300L, want_vcov = TRUE) {
  N <- nrow(X); q <- ncol(X)
  to_ab <- function(psi) {
    incr  <- if (J > 2L) exp(psi[2:(J - 1L)]) else numeric(0)
    alpha <- psi[1] + c(0, cumsum(incr))
    list(alpha = alpha, beta = if (q > 0L) psi[(J - 1L) + seq_len(q)] else numeric(0))
  }
  negll <- function(psi) {
    ab  <- to_ab(psi)
    lin <- if (q > 0L) as.numeric(X %*% ab$beta) else rep(0, N)
    G   <- stats::plogis(outer(lin, ab$alpha, `+`))
    Gf  <- cbind(0, G, 1)
    P   <- Gf[, 2:(J + 1L), drop = FALSE] - Gf[, 1:J, drop = FALSE]
    -sum(log(pmax(P[cbind(seq_len(N), y)], 1e-12)))
  }
  a0   <- stats::qlogis(seq_len(J - 1L) / J)
  psi0 <- c(a0[1], if (J > 2L) log(pmax(diff(a0), 1e-3)) else numeric(0), rep(0, q))
  opt  <- tryCatch(stats::optim(psi0, negll, method = "BFGS",
                                control = list(maxit = maxit)),
                   error = function(e) NULL)
  if (is.null(opt))
    return(list(theta_hat = rep(NA_real_, (J - 1L) + q),
                vcov = matrix(NA_real_, (J - 1L) + q, (J - 1L) + q),
                loglik = NA_real_, converged = FALSE))
  ab <- to_ab(opt$par); theta_hat <- c(ab$alpha, ab$beta)
  kk <- (J - 1L) + q
  vcov <- matrix(NA_real_, kk, kk)
  if (isTRUE(want_vcov)) {
    info <- matrix(0, kk, kk)
    for (i in seq_len(N)) info <- info + .info_cumulative(X[i, ], theta_hat, J)
    vcov <- tryCatch(solve(info), error = function(e) matrix(NA_real_, kk, kk))
  }
  list(theta_hat = theta_hat, vcov = vcov,
       loglik = -opt$value, converged = (opt$convergence == 0L))
}

.fit_model <- function(plan, Xrun, y, sigma, want_vcov = TRUE) {
  switch(plan$link,
         identity    = .fit_gaussian(Xrun, y),
         logit       = .fit_glm(Xrun, y, "logit"),
         loglinear   = .fit_glm(Xrun, y, "loglinear"),
         multinomial = .fit_multinomial(Xrun, y, plan$ncat),
         cumulative  = .fit_cumulative(Xrun, y, plan$ncat, want_vcov = want_vcov))
}

#' Simulate responses for a design under a formula-style model and fit it back.
#'
#' Given a design (support points and integer run \code{counts}), a formula-style
#' model (\code{link} plus \code{f}/\code{x}/\code{fx}/\code{ff}/\code{xx}, as in
#' \code{\link{optimal_design}}) and the TRUE parameters \code{theta}, simulate
#' the responses and (by default) obtain the maximum-likelihood estimate. With
#' \code{nsim > 1} the simulation-and-fit is repeated and the empirical
#' covariance of \eqn{\hat\theta} is compared to the design-based large-sample
#' covariance \eqn{M(\xi,\theta)^{-1}/N}.
#'
#' @param support \eqn{m \times d} matrix of distinct support points (one per
#'   row; factors as integer levels), OR an \code{"exact_design"} object (its
#'   \code{support}/\code{counts} are used).
#' @param counts integer vector of run counts, one per support point (ignored
#'   when \code{support} is an \code{"exact_design"}).
#' @param link,f,x,fx,ff,xx,intercept,coding,ncat the model, exactly as in
#'   \code{\link{optimal_design}}.
#' @param theta the TRUE parameter vector to simulate from (required); its length
#'   and ordering match the model (see \code{\link{optimal_design}}: stacked by
#'   category for \code{"multinomial"}, thresholds-then-slopes for
#'   \code{"cumulative"}).
#' @param design_box,factor_levels how to interpret the columns of \code{support}
#'   (which are factors and their levels). Supply one; if neither, factor columns
#'   are auto-detected (all-positive-integer columns), which may misclassify
#'   integer-valued continuous covariates.
#' @param sigma residual standard deviation for the \code{"identity"} (Gaussian)
#'   link (default 1; ignored otherwise).
#' @param nsim number of Monte Carlo replications (default 1).
#' @param seed optional integer seed for reproducibility.
#' @param fit if \code{TRUE} (default) also fit the model; if \code{FALSE} (only
#'   valid for \code{nsim = 1}) just return the simulated data.
#' @return For \code{nsim = 1}: a list with \code{data} (the simulated runs and
#'   response), and -- when \code{fit = TRUE} -- \code{theta_hat}, \code{se},
#'   \code{vcov}, \code{loglik}, \code{converged}. For \code{nsim > 1}: a list
#'   with \code{estimates} (\code{nsim} \eqn{\times} k), \code{theta_hat_mean},
#'   \code{bias}, \code{mse} (per-parameter mean squared error
#'   \eqn{\mathrm{E}[(\hat\theta-\theta)^2]} over the converged replicates,
#'   \eqn{= \mathrm{bias}^2 + \mathrm{variance}}), \code{cov_empirical},
#'   \code{cov_design} (\eqn{M(\xi,\theta)^{-1}/N}, times \eqn{\sigma^2} for
#'   identity), \code{se_empirical}, \code{se_design} and \code{n_converged}. Both carry
#'   \code{theta}, \code{coef_names}, \code{link} and \code{N}.
#' @seealso \code{\link{optimal_design}}, \code{\link{exact_design}},
#'   \code{\link{design_information}}.
#' @export
simulate_design <- function(support, counts = NULL, link = NULL,
                            f = NULL, x = NULL, fx = NULL, xx = NULL, ff = NULL,
                            intercept = TRUE, coding = "zero-sum", ncat = NULL,
                            theta = NULL, design_box = NULL, factor_levels = NULL,
                            sigma = 1, nsim = 1L, seed = NULL, fit = TRUE) {
  if (inherits(support, "exact_design")) {
    if (is.null(counts)) counts <- support$counts
    support <- support$support
  }
  support <- as.matrix(support); storage.mode(support) <- "double"
  m <- nrow(support)
  if (is.null(counts)) stop("'counts' (runs per support point) is required.",
                            call. = FALSE)
  counts <- as.integer(round(counts))
  if (length(counts) != m || any(counts < 0L) || sum(counts) < 1L)
    stop("'counts' must be one nonnegative integer per support point, summing ",
         "to at least 1.", call. = FALSE)
  if (is.null(link)) stop("'link' is required.", call. = FALSE)
  if (is.null(theta)) stop("'theta' (the true parameters) is required.",
                           call. = FALSE)
  nsim <- as.integer(nsim)
  if (nsim < 1L) stop("'nsim' must be >= 1.", call. = FALSE)
  if (!isTRUE(fit) && nsim > 1L)
    stop("fit = FALSE is only allowed with nsim = 1.", call. = FALSE)

  meta <- if (!is.null(design_box)) .parse_design_box(design_box)
          else .factor_meta_from_design(support, factor_levels)
  cov_names <- if (!is.null(design_box)) names(design_box) else colnames(support)
  plan <- .build_model_terms(meta, f, x, fx, xx, ff, intercept, coding, link,
                             cov_names, ncat = ncat)
  k <- plan$k
  theta <- as.numeric(theta)
  if (length(theta) != k)
    stop(sprintf("'theta' must have length %d for this model; got %d.",
                 k, length(theta)), call. = FALSE)
  if (plan$link == "cumulative") {
    a <- theta[seq_len(plan$ncat - 1L)]
    if (length(a) > 1L && any(diff(a) <= 0))
      stop("cumulative model: the first ncat-1 elements of 'theta' are the ",
           "thresholds and must be strictly increasing.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(as.integer(seed))

  # base model matrix per support point, then expanded to runs
  Fmat <- matrix(0, m, plan$base_k)
  for (i in seq_len(m)) Fmat[i, ] <- .model_row(support[i, ], plan)
  runidx <- rep(seq_len(m), counts)
  Xrun   <- Fmat[runidx, , drop = FALSE]
  N      <- sum(counts)

  col_nm <- if (!is.null(cov_names) && length(cov_names) == ncol(support))
              cov_names else paste0("x", seq_len(ncol(support)))

  # ---- single dataset --------------------------------------------------
  if (nsim == 1L) {
    y    <- .simulate_y(plan, Xrun, theta, sigma)
    dat  <- data.frame(support[runidx, , drop = FALSE], y = y)
    names(dat) <- c(col_nm, "y")
    out  <- list(data = dat, theta = theta, coef_names = plan$coef_names,
                 link = plan$link, ncat = plan$ncat, N = N)
    if (!isTRUE(fit)) return(out)
    fitr <- .fit_model(plan, Xrun, y, sigma)
    out$theta_hat <- fitr$theta_hat
    out$se        <- sqrt(diag(as.matrix(fitr$vcov)))
    out$vcov      <- fitr$vcov
    out$loglik    <- fitr$loglik
    out$converged <- fitr$converged
    if (!is.null(fitr$sigma_hat)) out$sigma_hat <- fitr$sigma_hat
    names(out$theta_hat) <- plan$coef_names
    names(out$se)        <- plan$coef_names
    return(out)
  }

  # ---- Monte Carlo -----------------------------------------------------
  E    <- matrix(NA_real_, nsim, k, dimnames = list(NULL, plan$coef_names))
  conv <- logical(nsim)
  for (s in seq_len(nsim)) {
    y  <- .simulate_y(plan, Xrun, theta, sigma)
    fr <- tryCatch(.fit_model(plan, Xrun, y, sigma, want_vcov = FALSE),
                   error = function(e) NULL)
    if (!is.null(fr) && isTRUE(fr$converged)) { E[s, ] <- fr$theta_hat; conv[s] <- TRUE }
  }

  # design-based large-sample covariance M(xi, theta)^{-1} / N
  fl  <- ifelse(meta$is_factor, meta$nlevels, NA_integer_)
  Mbar <- design_information(support, counts / N, link = link, f = f, x = x,
                             fx = fx, xx = xx, ff = ff, intercept = intercept,
                             coding = coding, ncat = ncat, theta = theta,
                             factor_levels = fl)
  scale <- if (plan$link == "identity") sigma^2 else 1
  cov_design <- tryCatch(scale * solve(N * Mbar),
                         error = function(e) matrix(NA_real_, k, k))
  dimnames(cov_design) <- list(plan$coef_names, plan$coef_names)

  ok <- which(conv)
  list(estimates = E, theta = theta, coef_names = plan$coef_names,
       link = plan$link, ncat = plan$ncat, N = N, nsim = nsim,
       n_converged = length(ok),
       theta_hat_mean = colMeans(E[ok, , drop = FALSE]),
       bias = colMeans(E[ok, , drop = FALSE]) - theta,
       mse = colMeans(sweep(E[ok, , drop = FALSE], 2, theta, "-")^2),
       cov_empirical = stats::cov(E[ok, , drop = FALSE]),
       cov_design = cov_design,
       se_empirical = apply(E[ok, , drop = FALSE], 2, stats::sd),
       se_design = sqrt(diag(cov_design)))
}
