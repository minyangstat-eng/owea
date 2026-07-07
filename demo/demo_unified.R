###############################################################################
# demo_unified.R
#
# Shows the features unique to the merged package:
#   (1) the SAME model solved via an information VECTOR (fast) and an
#       information MATRIX (general) -- identical optimal design;
#   (2) an EXISTING (two-stage) design via the fast information-vector path;
#   (3) OPTIONAL neighbourhood merging (default off).
###############################################################################

library(owea)

# Bi-exponential nonlinear model, theta = (1, 1, 1, 2).
theta <- c(1, 1, 1, 2)

# Information VECTOR  f(x, theta) = d eta / d theta   (per-point info = f f').
biexp_vec <- function(x, th) {
  xx <- x[1]
  c(            exp(-th[2] * xx),
    -th[1] * xx * exp(-th[2] * xx),
                exp(-th[4] * xx),
    -th[3] * xx * exp(-th[4] * xx))
}
# Information MATRIX  I_x = f f'  for the SAME model.
biexp_mat <- function(x) tcrossprod(biexp_vec(x, theta))

design_box    <- list(c(0.0, 3.0))
step_sequence <- c(0.05, 0.02, 0.01)

cat(strrep("=", 72), "\n")
cat("(1) Same model via information VECTOR vs information MATRIX\n")
cat(strrep("=", 72), "\n")

res_vec <- optimal_design(info_vector = biexp_vec, theta = theta,
                          design_box = design_box, step_sequence = step_sequence,
                          p = 0)
res_mat <- optimal_design(info_matrix = biexp_mat,
                          design_box = design_box, step_sequence = step_sequence,
                          p = 0)
print_result(res_vec, "Vector input (fast path)")
print_result(res_mat, "Matrix input (general path)")
cat(sprintf("\n   -> criterion difference (vector vs matrix) = %.2e\n",
            abs(res_vec$criterion - res_mat$criterion)))

cat("\n", strrep("=", 72), "\n", sep = "")
cat("(2) Two-stage design via the FAST information-vector path\n")
cat(strrep("=", 72), "\n")

res_2stage <- optimal_design(
  info_vector = biexp_vec, theta = theta,
  design_box  = design_box, step_sequence = step_sequence, p = 0,
  xi0_points  = matrix(c(0, 1, 2, 3), ncol = 1),
  xi0_weights = rep(0.25, 4), n0 = 40, n1 = 80)
print_result(res_2stage, "Two-stage D-optimal (n0=40, n1=80), vector input")

cat("\n", strrep("=", 72), "\n", sep = "")
cat("(3) Optional neighbourhood merging (default off)\n")
cat(strrep("=", 72), "\n")

res_nomerge <- optimal_design(info_vector = biexp_vec, theta = theta,
                              design_box = design_box,
                              step_sequence = step_sequence, p = 0,
                              merge = FALSE)
res_merge   <- optimal_design(info_vector = biexp_vec, theta = theta,
                              design_box = design_box,
                              step_sequence = step_sequence, p = 0,
                              merge = TRUE, merge_factor = 1.5)
cat(sprintf("support points  : merge=FALSE -> %d   merge=TRUE -> %d\n",
            nrow(res_nomerge$support), nrow(res_merge$support)))
cat(sprintf("criterion       : merge=FALSE -> %.6f   merge=TRUE -> %.6f\n",
            res_nomerge$criterion, res_merge$criterion))
