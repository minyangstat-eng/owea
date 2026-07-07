# =====================================================================
#  example_logistic3d.R
#
#  A worked tutorial: approximate optimal design for a 3-covariate
#  logistic regression model on a 3-D design region, using the
#  information-VECTOR (fast) API.  Demonstrates make_grid(),
#  appro_opt(X=, infor_vec=, theta=), and appro_opt_seq().
#
#  Run with:
#    source(system.file("examples", "example_logistic3d.R", package = "owea"))
# =====================================================================

library(owea)

# ---------------------------------------------------------------------
#  The model: logistic regression with linear predictor
#        eta = theta1 + theta2 x1 + theta3 x2 + theta4 x3
#  Information vector  f(x, theta) = scale(x) * (1, x1, x2, x3),
#  scale(x) = exp(eta/2) / (1 + exp(eta)).
# ---------------------------------------------------------------------
infor_vec <- function(x, theta) {
    out <- c(1, x[1], x[2], x[3])
    eta <- sum(out * theta)
    eta <- min(max(eta, -500), 500)
    (exp(eta / 2) / (1 + exp(eta))) * out
}

theta <- c(1, -0.5, 0.5, 1)
wb    <- diag(4)                              # full parameter vector
lower <- c(x1 = -2, x2 = -1, x3 = -3)
upper <- c(x1 =  2, x2 =  1, x3 =  3)

cat("=====================================================================\n")
cat(" 3-covariate logistic model -- D-optimal design\n")
cat(" theta =", theta, "\n")
cat("=====================================================================\n")


# =====================================================================
#  1.  make_grid()  +  appro_opt(X = , infor_vec = , theta = )
# =====================================================================
cat("\n[1] make_grid() + appro_opt(X = , infor_vec = , theta = )\n")
cat("    -----------------------------------------------------------\n")

X <- make_grid(lower, upper, by = 0.1)
cat(sprintf("    make_grid(step = 0.1)  ->  %d candidate points\n", nrow(X)))
cat("    first few grid points:\n")
print(head(X, 3))

tm1 <- system.time(
    res1 <- appro_opt(pp = 0, wb = wb, X = X, infor_vec = infor_vec,
                      theta = theta, tol = 1e-6)
)
cat(sprintf("\n    solved in %.2f s  (%d outer iterations)\n",
            tm1["elapsed"], res1$iter))
cat(sprintf("    criterion  log|Sigma|/4 = %.6f      max d = %.2e\n",
            res1$value, res1$sensitivity))
cat(sprintf("    %d support points (settings + weight):\n", nrow(res1$points)))
print(round(cbind(res1$points, weight = res1$weight), 4))


# =====================================================================
#  2.  appro_opt_seq()  --  sequential multi-resolution search
# =====================================================================
cat("\n[2] appro_opt_seq() : sequential coarse -> fine search\n")
cat("    -----------------------------------------------------------\n")

tm2 <- system.time(
    res2 <- appro_opt_seq(pp = 0, wb = wb, lower = lower, upper = upper,
                          by_seq = list(0.5, 0.1, 0.02),
                          infor_vec = infor_vec, theta = theta,
                          tol = 1e-6, verbose = TRUE)
)
cat(sprintf("\n    sequential search solved in %.2f s\n", tm2["elapsed"]))
cat(sprintf("    criterion  log|Sigma|/4 = %.6f   (final resolution 0.02)\n",
            res2$value))
cat(sprintf("    %d support points (settings + weight):\n", nrow(res2$points)))
print(round(res2$design, 4))


# =====================================================================
#  3.  infor_matrix()  --  information matrix of a design
# =====================================================================
cat("\n[3] infor_matrix() : information matrix of the optimal design\n")
cat("    -----------------------------------------------------------\n")

M <- infor_matrix(res2$points, res2$weight, infor_vec, theta)
cat("    M = sum_i weight_i * f(x_i,theta) f(x_i,theta)' :\n")
print(round(M, 5))

# Cross-check: for D-optimality of the full theta the reported criterion equals
# -log|M| / k, so exp(-criterion) recovers det(M)^(1/k).
cat(sprintf("\n    det(M)^(1/4)            = %.6f\n", det(M)^(1 / 4)))
cat(sprintf("    exp(-criterion)         = %.6f   (should match)\n",
            exp(-res2$value)))

cat("\n=====================================================================\n")
cat(" Done.\n")
cat("=====================================================================\n")
