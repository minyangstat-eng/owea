###############################################################################
# demo_general.R
#
# Demo: the general-purpose optimal_design() on Example 5 of Yang/Chen/Yang
# (2026) -- an A-optimal design under a logistic GLM with three continuous
# factors (information MATRIX input).
#
# To solve a different problem the user only changes:
#   - the info_matrix(x) function (the Fisher information matrix at x)
#   - the design_box
#   - the step_sequence
###############################################################################

library(owea)

# ---- 1. User supplies the information matrix ------------------------------
# Logistic GLM:  I_x = nu(x) * q(x) q(x)^T  with  nu(eta) = e^eta / (1+e^eta)^2.
beta <- c(1.0, -0.5, 0.5, 1.0)

logistic_info <- function(x) {
  q   <- c(1.0, x[1], x[2], x[3])     # predictor vector q(x)
  eta <- sum(q * beta)
  e   <- exp(eta)
  nu  <- e / (1 + e)^2                # GLM variance weight
  nu * tcrossprod(q)                  # k x k information matrix
}

# ---- 2. User supplies the design box --------------------------------------
design_box <- list(c(-2.0, 2.0), c(-1.0, 1.0), c(-3.0, 3.0))

# ---- 3. User supplies the grid sequence (coarse -> fine) ------------------
step_sequence <- c(0.1, 0.05, 0.02, 0.0125, 0.01)

cat(strrep("=", 72), "\n")
cat("Demo of general optimal_design() on Yang/Chen/Yang 2026 Example 5\n")
cat(strrep("=", 72), "\n")

res <- optimal_design(
  info_matrix   = logistic_info,
  design_box    = design_box,
  step_sequence = step_sequence,
  p             = 1,             # A-optimality
  verbose       = TRUE
)

print_result(res, title = "Example 5 A-optimal design")
