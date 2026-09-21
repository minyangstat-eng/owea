# ===========================================================================
# elfving.R -- certificate of c-optimality (a single linear combination).
#
# When the quantity of interest is one linear combination c'theta (wb = c', one
# row; a one-row grad_g; a one-parameter subset), the criterion is
#     V(w) = c' M(w)^- c ,   M(w) = infor0 + sum_i w_i A_i ,
# with A_i the (b-scaled) per-point information.  The c-optimal design is often
# SINGULAR (fewer support points than parameters), and the sensitivity function
# the package uses elsewhere -- built from a (pseudo-)inverse of M -- cannot
# certify a singular design: the general equivalence theorem then needs a
# specific generalized inverse (Pukelsheim, 1993, Thm 2.16 / 7.19).
#
# Elfving's theorem, in its dual (Lagrangian) form, gives a certificate that
# needs no inverse at all.  For EVERY h with h'c = 1 and every design w,
#     c' M(w)^- c  >=  1 / ( h' infor0 h + max_i h' A_i h ) ,
# because (h'c)^2 <= (h'Mh)(c'M^-c) and sum_i w_i = 1.  So the right-hand side
# is a lower bound on the c-optimal value V*, for any h; the best h attains V*
# (strong duality: the criterion is convex in w, and c'M^-c = max_h 2h'c - h'Mh
# is linear in w and concave in h).  A design of value V is certified optimal
# when V - bound <= eps0 * max(1, V), and efficiency >= bound / V always holds.
#
# Finding h: with information VECTORS and no existing design the problem is a
# linear program (min max_i |h'f_i| s.t. h'c = 1), solved exactly with lpSolve
# when it is installed; otherwise (matrix-mode information, an existing design,
# or no lpSolve) Nelder-Mead minimises the convex piecewise-quadratic
# h' infor0 h + max_i h' A_i h on the affine set h'c = 1, started from c/c'c and
# from the pseudo-inverse direction M^+ c.  Any h gives a VALID bound, so the
# fallback can only be looser than the LP, never wrong.
# ===========================================================================

# The LP: variables h = hp - hm (hp, hm >= 0) and t >= 0; minimise t subject to
# -t <= f_i' h <= t for every candidate i and c' h = 1.  Returns h, or NULL.
.c_lp <- function(cvec, F) {
  k <- nrow(F); N <- ncol(F); Ft <- t(F)
  A   <- rbind(cbind(Ft, -Ft, -1), cbind(-Ft, Ft, -1), c(cvec, -cvec, 0))
  sol <- lpSolve::lp("min", c(rep(0, 2 * k), 1), A,
                     c(rep("<=", 2 * N), "="), c(rep(0, 2 * N), 1))
  if (!identical(as.integer(sol$status), 0L)) return(NULL)
  sol$solution[seq_len(k)] - sol$solution[k + seq_len(k)]
}

# The certificate on a FINITE candidate set.
#   cvec      : the vector c (length k)
#   info_data : the b-scaled candidate information (k x N vectors, or k^2 x N)
#   infor0    : a * I_xi0 (k x k; zeros for a single-stage design)
#   value     : the design's value V = c' M^- c
#   M         : optional combined information of the design (for the M^+ c start)
# Returns list(bound, h, gap = V - bound, certified, efficiency = min(1, bound/V)).
.c_certificate <- function(cvec, info_mode, info_data, infor0, k, value,
                           eps0 = 1e-6, M = NULL) {
  cvec <- as.numeric(cvec); info_data <- as.matrix(info_data)
  infor0 <- as.matrix(infor0)
  quad_max <- function(h) {                          # max_i h' A_i h
    h <- as.numeric(h)
    if (info_mode == 0L) max(as.numeric(crossprod(info_data, h))^2)
    else max(as.numeric(crossprod(info_data, as.numeric(tcrossprod(h)))))
  }
  Dfun <- function(h) { h <- as.numeric(h); as.numeric(crossprod(h, infor0 %*% h)) + quad_max(h) }

  hs <- list()
  if (info_mode == 0L && all(infor0 == 0) &&
      requireNamespace("lpSolve", quietly = TRUE)) {
    h_lp <- tryCatch(.c_lp(cvec, info_data), error = function(e) NULL)
    if (!is.null(h_lp) && all(is.finite(h_lp))) hs[[length(hs) + 1L]] <- h_lp
  }
  h0 <- cvec / sum(cvec^2)                           # h'c = 1
  starts <- list(h0)
  if (!is.null(M)) {
    hp <- tryCatch(as.numeric(.cont_inv(as.matrix(M)) %*% cvec), error = function(e) NULL)
    if (!is.null(hp) && all(is.finite(hp)) && abs(sum(hp * cvec)) > 1e-12)
      starts[[length(starts) + 1L]] <- hp / sum(hp * cvec)
  }
  if (k > 1L) {                                      # minimise D on h'c = 1
    Nb <- qr.Q(qr(matrix(cvec, k, 1L)), complete = TRUE)[, -1L, drop = FALSE]
    for (s in starts) {
      z0 <- as.numeric(crossprod(Nb, s - h0))
      if (k == 2L) {                                 # one free direction: Brent
        sc <- max(1, abs(z0), sqrt(sum(h0^2)))
        o  <- tryCatch(stats::optimize(function(z) Dfun(h0 + Nb %*% z),
                                       c(z0 - 1e3 * sc, z0 + 1e3 * sc), tol = 1e-12),
                       error = function(e) NULL)
        if (!is.null(o) && is.finite(o$objective)) {
          o2 <- stats::optimize(function(z) Dfun(h0 + Nb %*% z),
                                c(o$minimum - sc, o$minimum + sc), tol = 1e-14)
          hs[[length(hs) + 1L]] <- as.numeric(h0 + Nb %*% o2$minimum)
        }
      } else {
        o <- tryCatch(stats::optim(z0, function(z) Dfun(h0 + Nb %*% z),
                                   method = "Nelder-Mead",
                                   control = list(maxit = 5000L, reltol = 1e-12)),
                      error = function(e) NULL)
        if (!is.null(o) && is.finite(o$value))
          hs[[length(hs) + 1L]] <- as.numeric(h0 + Nb %*% o$par)
      }
    }
  }
  hs <- c(hs, starts)
  Ds <- vapply(hs, Dfun, numeric(1))
  ok <- which(is.finite(Ds) & Ds > 0)
  if (!length(ok))
    return(list(bound = NA_real_, h = h0, gap = NA_real_, certified = FALSE,
                efficiency = NA_real_))
  best <- ok[which.min(Ds[ok])]
  h <- hs[[best]]; bound <- 1 / Ds[best]
  gap <- value - bound
  list(bound = bound, h = h, gap = gap,
       certified = is.finite(gap) && gap <= eps0 * max(1, abs(value)),
       efficiency = if (is.finite(bound) && is.finite(value) && value > 0)
                      min(1, bound / value) else NA_real_)
}

# max over the continuous region of h' A(x) h, from one start (L-BFGS-B on the
# continuous coordinates; factor levels fixed).  Used by the continuous path.
.cont_quad_max <- function(ev, h, x0, lo, hi) {
  cont <- ev$cont; k <- ev$k; hh <- as.numeric(tcrossprod(h))
  qf <- function(x) {
    f <- ev$col(x)
    if (ev$info_mode == 0L) sum(h * f)^2 else sum(hh * f)
  }
  full <- function(x1) { x <- x0; x[cont] <- x1; x }
  gr <- function(x1) {
    x <- full(x1); J <- ev$jac(x)
    if (ev$info_mode == 0L) -2 * sum(h * ev$col(x)) * as.numeric(crossprod(J, h))
    else -as.numeric(crossprod(J, hh))
  }
  o <- tryCatch(stats::optim(x0[cont], function(x1) -qf(full(x1)), gr,
                             method = "L-BFGS-B", lower = lo[cont], upper = hi[cont]),
                error = function(e) NULL)
  q0 <- qf(x0)
  if (is.null(o) || !is.finite(o$value)) q0 else max(q0, -o$value)
}

# The certificate for the CONTINUOUS design region: h is chosen on a finite set
# (the support, a 3-level grid and n_points random points), then the maximum
# of h' A(x) h is refined over the region by multi-start L-BFGS-B.  As every
# other check on the continuous path, the maximum is over the points examined.
.c_certificate_cont <- function(ev, cvec, infor0, lo, hi, is_factor, nlevels,
                                support, value, eps0 = 1e-6, n_points = 2000L) {
  k <- ev$k; d <- length(lo)
  G <- .cont_initial_grid(lo, hi, is_factor, nlevels, 3L)
  U <- matrix(0, n_points, d)
  for (j in seq_len(d))
    U[, j] <- if (is_factor[j]) sample.int(nlevels[j], n_points, replace = TRUE)
              else stats::runif(n_points, lo[j], hi[j])
  S <- rbind(as.matrix(support), G, U)
  C <- ev$cols(S)
  cert <- .c_certificate(cvec, ev$info_mode, C, infor0, k, value, eps0)
  h <- cert$h
  q <- if (ev$info_mode == 0L) as.numeric(crossprod(C, h))^2
       else as.numeric(crossprod(C, as.numeric(tcrossprod(h))))
  top <- order(q, decreasing = TRUE)[seq_len(min(10L, nrow(S)))]
  starts <- unique(rbind(as.matrix(support), S[top, , drop = FALSE]))
  qmax <- max(q)
  for (i in seq_len(nrow(starts))) {
    r <- .cont_quad_max(ev, h, starts[i, ], lo, hi)
    if (is.finite(r) && r > qmax) qmax <- r
  }
  D <- as.numeric(crossprod(h, as.matrix(infor0) %*% h)) + qmax
  bound <- if (is.finite(D) && D > 0) 1 / D else NA_real_
  gap <- value - bound
  list(bound = bound, h = h, gap = gap,
       certified = is.finite(gap) && gap <= eps0 * max(1, abs(value)),
       efficiency = if (is.finite(bound) && is.finite(value) && value > 0)
                      min(1, bound / value) else NA_real_,
       heuristic = TRUE)
}
