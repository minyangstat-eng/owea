###############################################################################
# demo_subset.R
#
# Two ways to target only PART of the parameter vector:
#   * subset = c(i, j, ...)        design optimal for theta[i], theta[j], ...
#   * grad_g = function(theta)     design optimal for a general differentiable
#                                  g(theta); returns the v x k Jacobian dg/dtheta'
#
# Both are accepted by DesignProblem() (with owea()) and by optimal_design().
# Supply AT MOST ONE; omitting both targets ALL parameters.
###############################################################################

library(owea)

# ---- Model: bi-exponential nonlinear regression (4 parameters) ------------
theta0 <- c(1.0, 1.0, 1.0, 2.0)

biexp_info <- function(x) {
  xx <- x[1]
  f <- c(            exp(-theta0[2] * xx),
         -theta0[1] * xx * exp(-theta0[2] * xx),
                     exp(-theta0[4] * xx),
         -theta0[3] * xx * exp(-theta0[4] * xx))
  tcrossprod(f)
}

###############################################################################
# A.  DesignProblem() + owea()  --  fixed candidate set
###############################################################################
cat(strrep("=", 72), "\n")
cat("A.  DesignProblem() + owea()  -- fixed candidate set\n")
cat(strrep("=", 72), "\n")

X <- matrix(seq(0, 3, length.out = 601), ncol = 1)

# --- A1. subset: A-optimal design for the subset (theta2, theta4) ----------
prob_A1 <- DesignProblem(X = X, theta = theta0, p = 1,
                         info_matrix = biexp_info, subset = c(2, 4))
res_A1  <- owea(prob_A1)
print_design("A1.  DesignProblem(subset = c(2, 4))", res_A1)

# --- A2. grad_g: the SAME design written with an explicit selector ---------
selector_24 <- function(theta) matrix(c(0, 1, 0, 0,
                                        0, 0, 0, 1), nrow = 2, byrow = TRUE)
prob_A2 <- DesignProblem(X = X, theta = theta0, p = 1,
                         info_matrix = biexp_info, grad_g = selector_24)
res_A2  <- owea(prob_A2)
print_design("A2.  DesignProblem(grad_g = selector for theta2, theta4)", res_A2)
cat(sprintf("\n   -> A1 (subset) and A2 (grad_g) agree: criterion diff = %.2e\n",
            abs(res_A1$criterion - res_A2$criterion)))

# --- A3. grad_g: a GENERAL function of the parameters ----------------------
# Interest in  g(theta) = ( theta1 ,  theta2 - theta4 ).
grad_g_general <- function(theta) matrix(c(1, 0, 0,  0,
                                           0, 1, 0, -1), nrow = 2, byrow = TRUE)
prob_A3 <- DesignProblem(X = X, theta = theta0, p = 1,
                         info_matrix = biexp_info, grad_g = grad_g_general)
res_A3  <- owea(prob_A3)
print_design("A3.  DesignProblem(grad_g for g = (theta1, theta2 - theta4))",
             res_A3)

###############################################################################
# B.  optimal_design()  --  continuous design space (multistage refinement)
###############################################################################
cat("\n")
cat(strrep("=", 72), "\n")
cat("B.  optimal_design()  -- continuous design space\n")
cat(strrep("=", 72), "\n")

design_box    <- list(c(0.0, 3.0))
step_sequence <- c(0.05, 0.02, 0.01, 0.005, 0.002, 0.001)

# --- B1. subset -----------------------------------------------------------
res_B1 <- optimal_design(info_matrix = biexp_info, design_box = design_box,
                         step_sequence = step_sequence, p = 1,
                         theta = theta0, subset = c(2, 4))
print_result(res_B1, title = "B1.  optimal_design(subset = c(2, 4))")

# --- B2. grad_g selector --------------------------------------------------
res_B2 <- optimal_design(info_matrix = biexp_info, design_box = design_box,
                         step_sequence = step_sequence, p = 1,
                         theta = theta0, grad_g = selector_24)
print_result(res_B2, title = "B2.  optimal_design(grad_g = selector)")
cat(sprintf("\n   -> B1 (subset) and B2 (grad_g) agree: criterion diff = %.2e\n",
            abs(res_B1$criterion - res_B2$criterion)))

# --- B3. grad_g general ---------------------------------------------------
res_B3 <- optimal_design(info_matrix = biexp_info, design_box = design_box,
                         step_sequence = step_sequence, p = 1,
                         theta = theta0, grad_g = grad_g_general)
print_result(res_B3,
             title = "B3.  optimal_design(grad_g for g = (theta1, theta2 - theta4))")

cat("\nSummary: subset = c(i, j) is the convenient form of a row-selector",
    "grad_g;\n         use grad_g directly for any other function g(theta).\n")
