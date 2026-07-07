###############################################################################
# demo_owea.R
#
# Demo: using owea() directly on a FIXED, finite candidate set (information
# MATRIX input).  Use owea() when the set of admissible settings is already a
# fixed finite list -- no continuous design region, no grid refinement.
#
# Model: the bi-exponential nonlinear regression of Yang/Biedermann/Tang
# (2013), Example 1,
#     Y ~ theta1*exp(-theta2*x) + theta3*exp(-theta4*x) + N(0, sigma^2),
# with x in [0, 3] and assumed (theta1, theta2, theta3, theta4) = (1, 1, 1, 2).
###############################################################################

library(owea)

# ---- 1. User supplies the information matrix ------------------------------
theta_ex1 <- c(1.0, 1.0, 1.0, 2.0)        # (theta1, theta2, theta3, theta4)

biexp_info <- function(x) {
  th1 <- theta_ex1[1]; th2 <- theta_ex1[2]
  th3 <- theta_ex1[3]; th4 <- theta_ex1[4]
  xx  <- x[1]
  f <- c(            exp(-th2 * xx),
         -th1 * xx * exp(-th2 * xx),
                     exp(-th4 * xx),
         -th3 * xx * exp(-th4 * xx))
  tcrossprod(f)                           # k x k information matrix
}

# ---- 2. User supplies the fixed candidate set -----------------------------
X <- matrix(seq(0.0, 3.0, length.out = 601L), ncol = 1L)

cat(strrep("=", 72), "\n")
cat("Demo of owea() on a fixed candidate set -- Yang/Biedermann/Tang Ex.1\n")
cat(strrep("=", 72), "\n")

# ---- 3. D-optimal design for all four parameters --------------------------
prob_D <- DesignProblem(X = X, theta = theta_ex1, p = 0,
                        info_matrix = biexp_info)
t_D    <- system.time(res_D <- owea(prob_D))[3]
print_design(sprintf("D-optimal design for theta  (time = %.4f s)", t_D), res_D)

# ---- 4. Variation: A-optimality for a SUBSET of parameters ----------------
prob_A <- DesignProblem(X = X, theta = theta_ex1, p = 1,
                        info_matrix = biexp_info, subset = c(2, 4))
res_A  <- owea(prob_A)
print_design("A-optimal design for (theta2, theta4)", res_A)

# ---- 5. Variation: a two-stage (multistage) design ------------------------
# An initial design xi0 with n0 runs has already been carried out; allocate n1
# further runs so the COMBINED design is D-optimal.
prob_2stage <- DesignProblem(
  X           = X,
  theta       = theta_ex1,
  p           = 0,
  info_matrix = biexp_info,
  xi0_points  = matrix(c(0.0, 1.0, 2.0, 3.0), ncol = 1L),  # already-run support
  xi0_weights = c(0.25, 0.25, 0.25, 0.25),                 # already-run weights
  n0          = 40,                                        # runs already done
  n1          = 80                                         # runs to allocate
)
res_2 <- owea(prob_2stage)
print_design("Two-stage D-optimal design (n0=40, n1=80)", res_2)

# ---- 6. Equivalence-theorem check -----------------------------------------
cat("\nEquivalence check (max d_p over X, should be <= eps0 = 1e-6):\n")
cat(sprintf("  D-optimal        : max d = %+.3e\n", res_D$max_d))
cat(sprintf("  A-optimal subset : max d = %+.3e\n", res_A$max_d))
cat(sprintf("  Two-stage        : max d = %+.3e\n", res_2$max_d))
