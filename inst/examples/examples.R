# =====================================================================
#  examples.R -- worked examples for the information-VECTOR (fast) API.
#
#  Run with:  source(system.file("examples", "examples.R", package = "owea"))
# =====================================================================

library(owea)

cat("\n")
cat("=====================================================================\n")
cat(" Approximate optimal design -- information-vector API\n")
cat("=====================================================================\n")


# ---------------------------------------------------------------------
#  Example 0 -- sanity check: simple linear regression, D-optimality
#   model     : E[y] = theta1 + theta2 * x ,   x in [-1, 1]
#   D-optimal : weight 1/2 at x = -1 and weight 1/2 at x = +1
# ---------------------------------------------------------------------
cat("\n--- Example 0: simple linear regression, D-optimality ---\n")

x_grid        <- seq(-1, 1, by = 0.01)        # 201 candidate points
infor_vec_all <- rbind(1, x_grid)             # k x N  (k = 2)
wb            <- diag(2)

res0 <- appro_opt(pp = 0, wb = wb, infor_vec_all = infor_vec_all,
                  max_iter = 100, tol = 1e-8)

cat(sprintf("support points x = %s\n",
            paste(sprintf("%.2f", x_grid[res0$index]), collapse = ", ")))
cat(sprintf("weights          = %s\n",
            paste(sprintf("%.4f", res0$weight), collapse = ", ")))
cat(sprintf("max directional derivative = %.2e   (iter = %d)\n",
            res0$sensitivity, res0$iter))
cat("expected: x = -1 and x = +1, each with weight 0.5\n")


# ---------------------------------------------------------------------
#  Example 1 -- ESD logistic experiment (7 parameters), D-optimality
# ---------------------------------------------------------------------
cat("\n--- Example 1: ESD logistic experiment, D-optimality ---\n")

build_logistic <- function(theta) {
    lev  <- c(-1, 1)
    vols <- seq(25, 45, by = 0.001)                 # 20001 voltage levels
    g    <- expand.grid(vol = vols, Pulse = lev, ESD = lev,
                        B = lev, A = lev)
    Fmat <- rbind(1, g$A, g$B, g$ESD, g$Pulse, g$vol, g$ESD * g$Pulse)
    xb   <- as.numeric(crossprod(theta, Fmat))      # theta' F
    xb   <- pmin(pmax(xb, -500), 500)               # guard exp() overflow
    scal <- exp(xb / 2) / (1 + exp(xb))
    list(infor_vec_all = Fmat * rep(scal, each = 7L),
         settings      = as.matrix(g[, c("A", "B", "ESD", "Pulse", "vol")]))
}

theta1 <- c(-7.5, 1.50, -0.2, -0.15, 0.25, 0.35, 0.4)
lg     <- build_logistic(theta1)
cat(sprintf("candidate grid: %d points (k = 7)\n", ncol(lg$infor_vec_all)))

wb1 <- diag(7)
t1  <- system.time(
    res1 <- appro_opt(pp = 0, wb = wb1, infor_vec_all = lg$infor_vec_all,
                      max_iter = 300, tol = 1e-6, verbose = FALSE)
)
cat(sprintf("solved in %.3f s,  %d outer iterations\n", t1["elapsed"], res1$iter))
cat(sprintf("max directional derivative = %.3e\n", res1$sensitivity))
cat(sprintf("criterion  log|Sigma| / 7  = %.6f\n", res1$value))
cat(sprintf("number of support points   = %d\n", length(res1$index)))

design_pts <- cbind(lg$settings[res1$index, , drop = FALSE], weight = res1$weight)
cat("optimal design (settings + weight):\n")
print(round(design_pts, 4))


# ---------------------------------------------------------------------
#  Example 2 -- linear main-effects model, A-optimality (slopes only)
# ---------------------------------------------------------------------
cat("\n--- Example 2: linear main-effects model, A-optimality ---\n")

set.seed(1)
Xlin          <- matrix(runif(6 * 5000), nrow = 6)   # 6 x 5000
infor_vec_lin <- rbind(1, Xlin)                      # k x N  (k = 7)
wb2           <- diag(7)[2:7, ]                      # exclude intercept

t2 <- system.time(
    res2 <- appro_opt(pp = 1, wb = wb2, infor_vec_all = infor_vec_lin,
                      max_iter = 1000, tol = 1e-6)
)
cat(sprintf("solved in %.3f s,  %d outer iterations\n", t2["elapsed"], res2$iter))
cat(sprintf("max directional derivative = %.3e\n", res2$sensitivity))
cat(sprintf("criterion  (tr Sigma / 6)  = %.6f\n", res2$value))
cat(sprintf("number of support points   = %d\n", length(res2$index)))


# ---------------------------------------------------------------------
#  Example 3 -- paper verification: Yang, Biedermann & Tang (2013),
#               Example 4, model (10).
#   Published optimal design (paper, p.1418):
#     (0, 0.3509), (0.3011, 0.4438), (0.7926, 0.1491), (1, 0.0562)
# ---------------------------------------------------------------------
cat("\n--- Example 3: paper Example 4 / model (10) -- verification ---\n")

th  <- c(1, 0.5, 1, 1)
xg  <- seq(0, 1, by = 1e-4)
infor_vec_m10 <- rbind(exp(th[2] * xg),
                       th[1] * xg * exp(th[2] * xg),
                       exp(th[4] * xg),
                       th[3] * xg * exp(th[4] * xg))
wb3 <- matrix(c(0.5, 1, 1, 1), nrow = 1)          # dg/dtheta  (v = 1)

res3 <- appro_opt(pp = 0, wb = wb3, infor_vec_all = infor_vec_m10,
                  max_iter = 200, tol = 1e-8)

# merge support points whose x-values are within 2.5e-3 of each other
ord <- order(xg[res3$index])
sx  <- xg[res3$index][ord]
sw  <- res3$weight[ord]
mx  <- numeric(0); mw <- numeric(0)
for (i in seq_along(sx)) {
    if (length(mx) && abs(sx[i] - mx[length(mx)]) < 2.5e-3) {
        nw <- mw[length(mw)] + sw[i]
        mx[length(mx)] <- (mx[length(mx)] * mw[length(mw)] + sx[i] * sw[i]) / nw
        mw[length(mw)] <- nw
    } else { mx <- c(mx, sx[i]); mw <- c(mw, sw[i]) }
}
cat(sprintf("max directional derivative = %.2e   (iter = %d)\n",
            res3$sensitivity, res3$iter))
cat("computed optimal design (x, weight):\n")
for (i in seq_along(mx))
    cat(sprintf("   x = %.4f   weight = %.4f\n", mx[i], mw[i]))
cat("paper target:\n")
cat("   x = 0.0000   weight = 0.3509\n")
cat("   x = 0.3011   weight = 0.4438\n")
cat("   x = 0.7926   weight = 0.1491\n")
cat("   x = 1.0000   weight = 0.0562\n")


# ---------------------------------------------------------------------
#  Example 4 -- 3-covariate logistic model on a fine 3-D grid
#   design region  x1 in [-2,2], x2 in [-1,1], x3 in [-3,3]  (step 0.05)
#   => 401,841 candidate points; D-optimal for the full parameter vector.
# ---------------------------------------------------------------------
cat("\n--- Example 4: 3-covariate logistic model, D-optimality ---\n")

build_logit3 <- function(theta) {
    grid <- 0.05
    x1 <- seq(-2, 2, by = grid)
    x2 <- seq(-1, 1, by = grid)
    x3 <- seq(-3, 3, by = grid)
    g  <- expand.grid(x1 = x1, x2 = x2, x3 = x3)
    Fmat <- rbind(1, g$x1, g$x2, g$x3)
    xb   <- as.numeric(crossprod(theta, Fmat))
    xb   <- pmin(pmax(xb, -500), 500)
    scal <- exp(xb / 2) / (1 + exp(xb))
    list(infor_vec_all = Fmat * rep(scal, each = 4L), settings = as.matrix(g))
}

theta4 <- c(1, -0.5, 0.5, 1)
lg3    <- build_logit3(theta4)
cat(sprintf("candidate grid: %d points (k = 4)\n", ncol(lg3$infor_vec_all)))

t4  <- system.time(
    res4 <- appro_opt(pp = 0, wb = diag(4), infor_vec_all = lg3$infor_vec_all,
                      max_iter = 100, tol = 1e-6)
)
cat(sprintf("solved in %.3f s,  %d outer iterations\n", t4["elapsed"], res4$iter))
cat(sprintf("max directional derivative = %.3e\n", res4$sensitivity))
cat(sprintf("criterion  log|Sigma| / 4  = %.6f\n", res4$value))
cat(sprintf("number of support points   = %d\n", length(res4$index)))

design_pts4 <- cbind(lg3$settings[res4$index, , drop = FALSE], weight = res4$weight)
cat("optimal design (settings + weight):\n")
print(round(design_pts4, 4))


# ---------------------------------------------------------------------
#  Example 5 -- convenience features: make_grid(), appro_opt(X=, ...),
#               and appro_opt_seq() (sequential multi-resolution).
# ---------------------------------------------------------------------
cat("\n--- Example 5: make_grid(), X interface, and sequential search ---\n")

infor_vec_logit3 <- function(x, theta) {
    out <- c(1, x[1], x[2], x[3])
    xb  <- sum(out * theta)
    xb  <- min(max(xb, -500), 500)
    (exp(xb / 2) / (1 + exp(xb))) * out
}
theta5 <- c(1, -0.5, 0.5, 1)

Xc <- make_grid(lower = c(-2, -1, -3), upper = c(2, 1, 3), by = 0.2)
cat(sprintf("(a) make_grid() at step 0.2  ->  %d candidate points\n", nrow(Xc)))
ta <- system.time(
    resa <- appro_opt(pp = 0, wb = diag(4),
                      X = Xc, infor_vec = infor_vec_logit3, theta = theta5)
)
cat(sprintf("    appro_opt(X=, infor_vec=, theta=): %.3f s, criterion %.6f\n",
            ta["elapsed"], resa$value))

tb <- system.time(
    resb <- appro_opt_seq(pp = 0, wb = diag(4),
                          lower  = c(-2, -1, -3),
                          upper  = c( 2,  1,  3),
                          by_seq = list(0.5, 0.1, 0.02),
                          infor_vec = infor_vec_logit3, theta = theta5,
                          verbose = TRUE)
)
cat(sprintf("(b) sequential search: %.3f s total\n", tb["elapsed"]))
cat(sprintf("    final criterion log|Sigma|/4 = %.6f  (grid resolution 0.02)\n",
            resb$value))
cat(sprintf("    number of support points = %d\n", nrow(resb$points)))
cat("    optimal design (settings + weight):\n")
print(round(resb$design, 4))


cat("\n=====================================================================\n")
cat(" All examples finished.\n")
cat("=====================================================================\n")
