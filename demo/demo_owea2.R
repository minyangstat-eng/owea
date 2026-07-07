###############################################################################
# demo_owea2.R
#
# Demo: owea() on a FIXED candidate grid built with candidate_grid()
# (information MATRIX input).
#
# Model: a logistic GLM with three continuous factors (Yang/Chen/Yang 2026,
# Example 5), solved on a single pre-discretized grid instead of with the
# multistage refinement of optimal_design().
#
# NOTE ON GRID SIZE
#   This demo uses step 0.05 (~0.4 million points), which runs in a few
#   seconds.  For a genuinely fine solution use the multistage optimal_design()
#   (see demo_general.R): it reaches very fine steps cheaply because the fine
#   stages only scan small neighbourhoods.
###############################################################################

library(owea)

# ---- 1. User supplies the information matrix ------------------------------
beta <- c(1.0, -0.5, 0.5, 1.0)

logistic_info <- function(x) {
  q   <- c(1.0, x[1], x[2], x[3])
  eta <- sum(q * beta)
  e   <- exp(eta)
  nu  <- e / (1 + e)^2
  nu * tcrossprod(q)
}

# ---- 2. User supplies the fixed candidate set -----------------------------
design_box <- list(c(-2.0, 2.0), c(-1.0, 1.0), c(-3.0, 3.0))
X <- candidate_grid(design_box, 0.05)            # 81 x 41 x 121 points

cat(strrep("=", 72), "\n")
cat("Demo of owea() on a fixed candidate grid -- logistic GLM (Example 5)\n")
cat(sprintf("candidate set: %d points (step 0.05)\n", nrow(X)))
cat(strrep("=", 72), "\n")

# ---- 3. D-optimal design --------------------------------------------------
prob_D <- DesignProblem(X = X, theta = beta, p = 0,
                        info_matrix = logistic_info)
t_D    <- system.time(res_D <- owea(prob_D))[3]
print_design(sprintf("D-optimal design  (time = %.4f s)", t_D), res_D)

# ---- 4. A-optimal design --------------------------------------------------
prob_A <- DesignProblem(X = X, theta = beta, p = 1,
                        info_matrix = logistic_info)
t_A    <- system.time(res_A <- owea(prob_A))[3]
print_design(sprintf("A-optimal design  (time = %.4f s)", t_A), res_A)

cat("\nEquivalence check (max d_p over X, should be <= eps0 = 1e-6):\n")
cat(sprintf("  D-optimal : max d = %+.3e\n", res_D$max_d))
cat(sprintf("  A-optimal : max d = %+.3e\n", res_A$max_d))
