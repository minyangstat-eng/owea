###############################################################################
# demo_exact.R
#
# Demo: exact_design() turns the approximate Phi_p-optimal design into a
# highly efficient EXACT design for a user-specified sample size n.
#
#   1. round the approximate weights to integer counts (apportionment),
#   2. repair the total to exactly n using the sensitivity function,
#   3. improve by random exchanges,
#   4. report the efficiency relative to the approximate optimum.
#
# To solve a different problem the user changes only the model
# (info_matrix / info_vector), the design region, and the sample size n.
###############################################################################

library(owea)

# ---- A. Simple linear regression on [-1, 1] (D-optimal) -------------------
# The approximate D-optimal design puts weight 1/2 at each of -1 and +1, so for
# an even n the exact design splits evenly and is exactly optimal (efficiency 1).
lin_info <- function(x) tcrossprod(c(1, x[1]))

cat(strrep("=", 72), "\n")
cat("Exact D-optimal design for simple linear regression, n = 10\n")
cat(strrep("=", 72), "\n")
r1 <- exact_design(n = 10, info_matrix = lin_info,
                   design_box = list(c(-1, 1)), step_sequence = c(0.1, 0.05),
                   p = 0, seed = 1)
print(r1)

# ---- B. Bi-exponential model on a fixed candidate set (D-optimal) ----------
# I_x = f f' with the bi-exponential gradient f(x, theta).
theta_biexp <- c(1.0, 1.0, 1.0, 2.0)
biexp_vec <- function(x, th)
  c(           exp(-th[2] * x[1]),
    -th[1] * x[1] * exp(-th[2] * x[1]),
                 exp(-th[4] * x[1]),
    -th[3] * x[1] * exp(-th[4] * x[1]))
biexp_info <- function(x) tcrossprod(biexp_vec(x, theta_biexp))

X <- matrix(seq(0, 3, length.out = 301), ncol = 1)

cat("\n", strrep("=", 72), "\n", sep = "")
cat("Exact D-optimal bi-exponential designs: efficiency rises with n\n")
cat(strrep("=", 72), "\n")
for (n in c(6, 10, 25, 100)) {
  e <- exact_design(n = n, candidate_set = X, theta = theta_biexp,
                    info_matrix = biexp_info, p = 0, seed = 42)
  cat(sprintf("n = %3d  |support| = %d  efficiency = %.4f\n",
              n, nrow(e$support), e$efficiency))
}

cat("\nFull printout for n = 25:\n")
print(exact_design(n = 25, candidate_set = X, theta = theta_biexp,
                   info_matrix = biexp_info, p = 0, seed = 42))
