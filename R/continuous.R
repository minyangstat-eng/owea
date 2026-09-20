# ===========================================================================
# continuous.R -- grid-free search over CONTINUOUS covariates.
#
# optimal_design(continuous = TRUE) replaces the grid discretisation of the
# continuous covariates by a search in the continuous design region.  Factor
# covariates keep their levels; only continuous coordinates move.
#
#   1. initial design : the OWEA engine solved on a SMALL coarse grid
#                       (init_levels equally spaced levels per continuous
#                       covariate x all levels of every factor), or the
#                       user-supplied init_points;
#   2. polishing      : with the weights held fixed, the criterion is minimised
#                       over the continuous coordinates of ALL support points by
#                       L-BFGS-B (stats::optim), using the closed-form gradient
#                          d crit / d x_i = -(w_i / v) d tr(part I(x)) / dx |_{x_i}
#                       with `part` the matrix behind the directional derivative
#                       (part_coeff() in owea_engine.cpp) and v the number of
#                       parameters of interest;
#   3. weights        : the package's own Newton weight step (the active-set
#                       engine) on the current support;
#   4. merging        : support points closer than merge_tol (a fraction of each
#                       covariate's range, factor levels identical) are fused;
#   5. new point      : the maximiser of the directional derivative d(x, xi) over
#                       the region by multi-start L-BFGS-B (starts: the support
#                       points, the box vertices, n_starts random points) for
#                       every level combination of the factor covariates;
#   6. stopping       : max d <= eps0 AND a random audit of n_audit points (its
#                       best points polished by L-BFGS-B) finds no violation.
#                       Otherwise the violating point is added and the loop
#                       continues.
#
# Derivatives of the per-point information with respect to the continuous
# coordinates come from, in this order, the user's `info_jacobian`, the analytic
# Jacobian a formula-style model carries (attribute "jacobian", see model.R),
# or finite differences of the model function (central inside the box,
# one-sided at its boundary).
#
# The certificate of optimality is the largest directional derivative found over
# every point examined (multi-start maximisations + audit), NOT an exhaustive
# scan; efficiency_lower_bound is computed from it and documented as such.
# ===========================================================================

# ---- inverse of a symmetric positive (semi)definite matrix -----------------
# Cholesky first, then an increasing ridge, then a pseudo-inverse (mirrors
# spd_inv() of owea_engine.cpp).
.cont_inv <- function(A) {
  A <- 0.5 * (A + t(A))
  R <- tryCatch(chol2inv(chol(A)), error = function(e) NULL)
  if (!is.null(R)) return(R)
  n <- nrow(A); reg <- 1e-12 * max(1, sum(diag(A)) / n)
  for (t in 1:12) {
    R <- tryCatch(chol2inv(chol(A + reg * diag(n))), error = function(e) NULL)
    if (!is.null(R)) return(R)
    reg <- reg * 10
  }
  s <- svd(A); pos <- s$d > 1e-12 * max(s$d)
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

# ---- `part` matrix and `coeff` of the directional derivative --------------
# pp = 0 : part = T' (T wb')^{-1} T, coeff = 1 ;  pp = 1 : part = T'T, coeff = 1/v
# with T = wb M^{-1}  (part_coeff() of owea_engine.cpp).
.cont_part <- function(M, wb, pp) {
  Tm <- wb %*% .cont_inv(M)
  if (as.integer(pp) == 0L)
    list(part = crossprod(Tm, .cont_inv(Tm %*% t(wb)) %*% Tm), coeff = 1)
  else
    list(part = crossprod(Tm), coeff = 1 / nrow(wb))
}

# Candidate-design information from b-scaled storage columns C (k x m or k^2 x m).
.cont_M <- function(C, w, info_mode, k) {
  M <- if (info_mode == 0L) C %*% (w * t(C)) else matrix(C %*% w, k, k)
  0.5 * (M + t(M))
}

# Directional derivatives d(x_n) for the storage columns C, given part/coeff and
# diff1 = tr(part * opt_infor).  Entry n is coeff * (tr(part I_n) - diff1).
.cont_dir <- function(C, part, coeff, diff1, info_mode) {
  q <- if (info_mode == 0L) colSums(C * (part %*% C))
       else as.numeric(as.numeric(part) %*% C)
  coeff * (q - diff1)
}

# Gradient of tr(part I(x)) with respect to the continuous coordinates, from
# the storage column f at x and the Jacobian J of that column (len x ncont).
.cont_tr_grad <- function(f, J, part, info_mode) {
  if (info_mode == 0L) 2 * as.numeric(crossprod(J, part %*% f))
  else as.numeric(as.numeric(part) %*% J)
}

# ---- model evaluator --------------------------------------------------------
# Bundles the per-point information, its b-scaled storage column, a vectorised
# evaluation over many points, and the Jacobian of the storage column with
# respect to the CONTINUOUS coordinates (columns in design_box order of the
# continuous covariates).  info_vector / info_matrix are the NORMALISED model
# functions (info_vector(x, theta), info_matrix(x)).
.cont_evaluator <- function(info_mode, info_vector, info_matrix, theta, k, b,
                            is_factor, lo, hi, info_jacobian = NULL) {
  cont  <- which(!is_factor); ncont <- length(cont)
  len   <- if (info_mode == 0L) k else k * k
  sc    <- if (b == 1) 1 else if (info_mode == 0L) sqrt(b) else b
  raw_col <- function(x)
    as.numeric(.info_col_at(x, info_mode, info_vector, info_matrix, theta))
  col <- function(x) sc * raw_col(x)

  fn_vec <- if (info_mode == 0L) attr(info_vector, "vectorized")
            else attr(info_matrix, "vectorized")
  cols <- function(X) {
    X <- as.matrix(X); storage.mode(X) <- "double"
    if (nrow(X) == 0L) return(matrix(0, len, 0L))
    out <- NULL
    if (is.function(fn_vec)) {
      out <- tryCatch(as.matrix(fn_vec(X, theta)), error = function(e) NULL)
      if (!is.null(out) && !identical(dim(out), c(len, nrow(X)))) out <- NULL
    }
    if (is.null(out))
      out <- .build_info_data(X, info_mode, info_vector, info_matrix, theta)$info_data
    sc * out
  }

  jac_fun <- NULL
  if (!is.null(info_jacobian)) {
    if (!is.function(info_jacobian))
      stop("'info_jacobian' must be a function(x, theta) (or function(x)).",
           call. = FALSE)
    np <- tryCatch(length(formals(info_jacobian)), error = function(e) 1L)
    jac_fun <- if (np >= 2L) info_jacobian else function(x, theta) info_jacobian(x)
    # validate the shape once, up front: inside the optimisers a failure would
    # only be caught and fallen back from silently
    probe <- .factor_probe_point(lo, hi, is_factor)
    J0 <- tryCatch(as.matrix(jac_fun(probe, theta)), error = function(e)
      stop("'info_jacobian' failed at a probe point: ", conditionMessage(e), call. = FALSE))
    if (nrow(J0) != len || ncol(J0) != ncont)
      stop(sprintf(paste0("'info_jacobian' must return a %d x %d matrix (rows: the ",
                          "information %s, columns: the continuous covariates in ",
                          "design_box order); got %d x %d."),
                   len, ncont, if (info_mode == 0L) "vector" else "matrix (column-major)",
                   nrow(J0), ncol(J0)), call. = FALSE)
  } else {
    a <- if (info_mode == 0L) attr(info_vector, "jacobian") else attr(info_matrix, "jacobian")
    if (is.function(a)) jac_fun <- a
  }
  jac <- if (is.function(jac_fun)) {
    function(x) {
      J <- as.matrix(jac_fun(x, theta))
      if (nrow(J) != len || ncol(J) != ncont)
        stop(sprintf(paste0("the information Jacobian must be a %d x %d matrix ",
                            "(rows: the information %s, columns: the continuous ",
                            "covariates in design_box order); got %d x %d."),
                     len, ncont, if (info_mode == 0L) "vector" else "matrix (column-major)",
                     nrow(J), ncol(J)), call. = FALSE)
      sc * J
    }
  } else {
    function(x) {                              # finite differences
      J  <- matrix(0, len, ncont)
      f0 <- NULL
      for (jj in seq_len(ncont)) {
        j <- cont[jj]; h <- 1e-6 * max(1, abs(x[j]))
        xp <- x; xm <- x
        if (x[j] + h > hi[j]) {                # one-sided at the upper bound
          if (is.null(f0)) f0 <- raw_col(x)
          xm[j] <- x[j] - h; J[, jj] <- (f0 - raw_col(xm)) / h
        } else if (x[j] - h < lo[j]) {         # one-sided at the lower bound
          if (is.null(f0)) f0 <- raw_col(x)
          xp[j] <- x[j] + h; J[, jj] <- (raw_col(xp) - f0) / h
        } else {
          xp[j] <- x[j] + h; xm[j] <- x[j] - h
          J[, jj] <- (raw_col(xp) - raw_col(xm)) / (2 * h)
        }
      }
      sc * J
    }
  }
  list(col = col, cols = cols, jac = jac, len = len, k = k, info_mode = info_mode,
       cont = cont, ncont = ncont, analytic = is.function(jac_fun))
}

# ---- the small coarse grid of the initial design -----------------------------
.cont_initial_grid <- function(lo, hi, is_factor, nlevels, init_levels = 3L) {
  L <- max(2L, as.integer(init_levels))
  axes <- lapply(seq_along(lo), function(j) {
    if (is_factor[j]) seq_len(nlevels[j])
    else if (hi[j] > lo[j]) seq(lo[j], hi[j], length.out = L)
    else lo[j]
  })
  X <- unname(as.matrix(expand.grid(axes)))
  storage.mode(X) <- "double"
  X
}

# ---- weights on a fixed (off-grid) support ---------------------------------
# The package's Newton weight step in its active-set form, warm-started from
# equal weights, with the support itself as the candidate set (the exchange
# loop then has nothing to add).  Returns the kept indices, their weights and
# the criterion value.
.cont_weights <- function(pp, wb, info_mode, C, infor0, msup, eps0) {
  m <- ncol(C)
  if (m == 1L)
    return(list(index = 1L, weight = 1,
                value = criterion_cpp(as.integer(pp), 1L, 1.0, as.integer(info_mode),
                                      C, wb, infor0)))
  r <- appro_opt_cpp(as.integer(pp), wb, as.integer(info_mode), C, infor0,
                     seq_len(m), as.integer(msup), 100L, as.numeric(eps0), FALSE,
                     1L, 1L)
  list(index = as.integer(r$index), weight = as.numeric(r$weight), value = r$value)
}

# ---- polishing of the support-point locations ------------------------------
# Weights fixed; every continuous coordinate of every support point is free
# within the box.  Never returns a design worse than the input.
.cont_polish <- function(ev, X, w, pp, wb, infor0, lo, hi, is_factor) {
  cont <- ev$cont; ncont <- ev$ncont; m <- nrow(X); v <- nrow(wb)
  if (m == 0L || ncont == 0L) return(X)
  info_mode <- ev$info_mode; k <- ev$k
  par0  <- as.numeric(X[, cont])
  lower <- rep(lo[cont], each = m); upper <- rep(hi[cont], each = m)
  put <- function(par) { Xm <- X; Xm[, cont] <- par; Xm }
  fn <- function(par) {
    val <- criterion_cpp(as.integer(pp), seq_len(m), w, as.integer(info_mode),
                         ev$cols(put(par)), wb, infor0)
    if (!is.finite(val)) 1e300 else val
  }
  gr <- function(par) {
    Xm <- put(par); C <- ev$cols(Xm)
    pc <- .cont_part(infor0 + .cont_M(C, w, info_mode, k), wb, pp)
    G  <- matrix(0, m, ncont)
    for (i in seq_len(m))
      if (w[i] > 0)
        G[i, ] <- -(w[i] / v) * .cont_tr_grad(C[, i], ev$jac(Xm[i, ]), pc$part, info_mode)
    as.numeric(G)
  }
  f0 <- fn(par0)
  o  <- tryCatch(stats::optim(par0, fn, gr, method = "L-BFGS-B", lower = lower,
                              upper = upper, control = list(factr = 1e3, maxit = 200L)),
                 error = function(e) NULL)
  if (is.null(o) || !is.finite(o$value) || o$value > f0) return(X)
  put(pmin(pmax(o$par, lower), upper))
}

# ---- merging of nearly coincident support points ----------------------------
# Distance = Euclidean distance of the continuous coordinates scaled by the
# range of each covariate; only points with identical factor levels merge, and
# the levels are kept exact.
.cont_merge <- function(X, w, is_factor, lo, hi, tol) {
  cont <- which(!is_factor)
  if (!length(cont) || !is.finite(tol) || tol <= 0) return(list(X = X, w = w))
  rng  <- pmax(hi[cont] - lo[cont], .Machine$double.eps)
  repeat {
    m <- nrow(X)
    if (m < 2L) break
    S <- sweep(X[, cont, drop = FALSE], 2, rng, "/")
    D <- as.matrix(stats::dist(S)); D[lower.tri(D, diag = TRUE)] <- Inf
    if (any(is_factor)) {
      same <- as.matrix(stats::dist(X[, is_factor, drop = FALSE])) < 0.5
      D[!same] <- Inf
    }
    kk <- which.min(D)
    if (!is.finite(D[kk]) || D[kk] >= tol) break
    ij <- arrayInd(kk, dim(D)); i <- ij[1]; j <- ij[2]
    xn <- X[i, ]
    xn[cont] <- (w[i] * X[i, cont] + w[j] * X[j, cont]) / (w[i] + w[j])
    X <- rbind(X[-c(i, j), , drop = FALSE], xn); w <- c(w[-c(i, j)], w[i] + w[j])
  }
  list(X = unname(X), w = w)
}

# ---- local maximisation of d(x) from one start (continuous coords only) ----
.cont_local_max <- function(ev, x_start, part, coeff, diff1, lo, hi) {
  cont <- ev$cont; info_mode <- ev$info_mode; k <- ev$k
  dval <- function(x) {
    f <- ev$col(x)
    q <- if (info_mode == 0L) sum(f * (part %*% f)) else sum(part * matrix(f, k, k))
    coeff * (q - diff1)
  }
  full <- function(x1) { x <- x_start; x[cont] <- x1; x }
  fn <- function(x1) { v <- -dval(full(x1)); if (is.finite(v)) v else 1e300 }
  gr <- function(x1) { x <- full(x1); -coeff * .cont_tr_grad(ev$col(x), ev$jac(x), part, info_mode) }
  x1 <- x_start[cont]
  o  <- tryCatch(stats::optim(x1, fn, gr, method = "L-BFGS-B", lower = lo[cont],
                              upper = hi[cont]),
                 error = function(e) NULL)
  if (is.null(o) || !is.finite(o$value)) return(list(x = x_start, d = dval(x_start)))
  x <- full(pmin(pmax(o$par, lo[cont]), hi[cont]))
  list(x = x, d = -o$value)
}

# ---- Step 5: the maximiser of the directional derivative over the region ---
# Multi-start L-BFGS-B for every level combination of the factor covariates.
.cont_new_point <- function(ev, X, part, coeff, diff1, lo, hi, is_factor, nlevels,
                            n_starts = 30L) {
  d <- length(lo); cont <- ev$cont; fac <- which(is_factor)
  combos <- if (length(fac))
              as.matrix(expand.grid(lapply(fac, function(j) seq_len(nlevels[j]))))
            else matrix(numeric(0), 1L, 0L)
  ncombo <- nrow(combos)
  verts  <- if (length(cont) <= 4L)
              as.matrix(expand.grid(lapply(cont, function(j) c(lo[j], hi[j]))))
            else NULL
  n_rand <- max(3L, ceiling(as.integer(n_starts) / ncombo))
  best <- list(d = -Inf, x = NULL)
  for (r in seq_len(ncombo)) {
    x2 <- if (length(fac)) combos[r, ] else numeric(0)
    base <- numeric(d); base[fac] <- x2
    starts <- list()
    if (nrow(X)) {
      own <- if (length(fac))
               which(apply(X[, fac, drop = FALSE], 1, function(z) all(abs(z - x2) < 0.5)))
             else seq_len(nrow(X))
      for (i in own) starts[[length(starts) + 1L]] <- X[i, cont]
    }
    if (!is.null(verts)) for (i in seq_len(nrow(verts))) starts[[length(starts) + 1L]] <- verts[i, ]
    for (i in seq_len(n_rand))
      starts[[length(starts) + 1L]] <- stats::runif(length(cont), lo[cont], hi[cont])
    for (s in starts) {
      x0 <- base; x0[cont] <- s
      res <- .cont_local_max(ev, x0, part, coeff, diff1, lo, hi)
      if (is.finite(res$d) && res$d > best$d) best <- res
    }
  }
  best
}

# ---- the random audit --------------------------------------------------------
# d(x) at n random points of the region; the n_top best are polished locally.
.cont_audit <- function(ev, n, part, coeff, diff1, lo, hi, is_factor, nlevels,
                        n_top = 5L) {
  d <- length(lo); n <- as.integer(n)
  U <- matrix(0, n, d)
  for (j in seq_len(d))
    U[, j] <- if (is_factor[j]) sample.int(nlevels[j], n, replace = TRUE)
              else stats::runif(n, lo[j], hi[j])
  dd <- .cont_dir(ev$cols(U), part, coeff, diff1, ev$info_mode)
  dd[!is.finite(dd)] <- -Inf
  raw_max <- max(dd)
  best <- list(d = raw_max, x = U[which.max(dd), ])
  for (t in order(dd, decreasing = TRUE)[seq_len(min(n_top, n))]) {
    res <- .cont_local_max(ev, U[t, ], part, coeff, diff1, lo, hi)
    if (is.finite(res$d) && res$d > best$d) best <- res
  }
  list(max_d = best$d, x = best$x, raw_max_d = raw_max, n = n)
}

# ---- maximum directional derivative of a GIVEN design over the region -------
# Used by verify_optimality(continuous = TRUE) and by the solver's final check.
.cont_max_sensitivity <- function(ev, X, w, pp, wb, infor0, lo, hi, is_factor,
                                  nlevels, n_starts, n_audit) {
  C <- ev$cols(X)
  opt_infor <- .cont_M(C, w, ev$info_mode, ev$k)
  pc <- .cont_part(infor0 + opt_infor, wb, pp)
  diff1 <- sum(pc$part * opt_infor)
  np <- .cont_new_point(ev, X, pc$part, pc$coeff, diff1, lo, hi, is_factor, nlevels,
                        n_starts)
  au <- NULL
  if (n_audit > 0) {
    au <- .cont_audit(ev, n_audit, pc$part, pc$coeff, diff1, lo, hi, is_factor, nlevels)
    if (au$max_d > np$d) np <- list(d = au$max_d, x = au$x)
  }
  list(max_d = np$d, x = np$x, audit = au)
}

# ---- the solver ---------------------------------------------------------------
.cont_solve <- function(ev, pp, wb, infor0, msup, lo, hi, is_factor, nlevels,
                        init_support, init_weights, eps0, max_iter, n_starts,
                        n_audit, merge_tol, verbose = FALSE) {
  t0 <- proc.time()[3]
  info_mode <- ev$info_mode; k <- ev$k
  X <- as.matrix(init_support); storage.mode(X) <- "double"
  ws <- .cont_weights(pp, wb, info_mode, ev$cols(X), infor0, msup, eps0)
  X <- X[ws$index, , drop = FALSE]; w <- ws$weight
  hist <- list(); converged <- FALSE; audit <- NULL
  x_new <- NULL; max_d <- NA_real_; crit <- ws$value
  best_crit <- Inf; stall <- 0L
  for (it in seq_len(max_iter)) {
    ## 2. polish the locations (weights fixed), 3. re-optimise the weights
    X  <- .cont_polish(ev, X, w, pp, wb, infor0, lo, hi, is_factor)
    ws <- .cont_weights(pp, wb, info_mode, ev$cols(X), infor0, msup, eps0)
    X  <- X[ws$index, , drop = FALSE]; w <- ws$weight
    ## 4. merge
    mg <- .cont_merge(X, w, is_factor, lo, hi, merge_tol)
    if (nrow(mg$X) < nrow(X)) {
      ws <- .cont_weights(pp, wb, info_mode, ev$cols(mg$X), infor0, msup, eps0)
      X  <- mg$X[ws$index, , drop = FALSE]; w <- ws$weight
    }
    crit <- ws$value
    ## 5. the maximiser of the directional derivative
    C <- ev$cols(X); opt_infor <- .cont_M(C, w, info_mode, k)
    pc <- .cont_part(infor0 + opt_infor, wb, pp)
    diff1 <- sum(pc$part * opt_infor)
    np <- .cont_new_point(ev, X, pc$part, pc$coeff, diff1, lo, hi, is_factor, nlevels,
                          n_starts)
    max_d <- np$d; x_new <- np$x
    ## 6. the audit, only when the search itself finds no violation
    if (is.finite(max_d) && max_d <= eps0 && n_audit > 0) {
      audit <- .cont_audit(ev, n_audit, pc$part, pc$coeff, diff1, lo, hi, is_factor,
                           nlevels)
      if (audit$max_d > max_d) { max_d <- audit$max_d; x_new <- audit$x }
    }
    hist[[it]] <- data.frame(iter = it, time = proc.time()[3] - t0, support = nrow(X),
                             criterion = crit, max_d = max_d)
    if (verbose)
      cat(sprintf("  continuous iter %3d: |S| = %2d  crit = %.8f  max d = %.3e  (%.2f s)\n",
                  it, nrow(X), crit, max_d, proc.time()[3] - t0))
    if (is.finite(max_d) && max_d <= eps0) { converged <- TRUE; break }
    ## stall guard: no criterion progress for several rounds
    if (!is.finite(best_crit) || crit < best_crit - 1e-12 * max(1, abs(best_crit))) {
      best_crit <- crit; stall <- 0L
    } else stall <- stall + 1L
    if (stall >= 6L) break
    ## add the violating point and re-optimise the weights
    X  <- rbind(X, x_new)
    ws <- .cont_weights(pp, wb, info_mode, ev$cols(X), infor0, msup, eps0)
    X  <- X[ws$index, , drop = FALSE]; w <- ws$weight
  }
  list(support = unname(X), weights = w, criterion = crit, max_d = max_d,
       converged = converged, iterations = it, maximiser = x_new, audit = audit,
       history = do.call(rbind, hist), time = proc.time()[3] - t0)
}
