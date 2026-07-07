###############################################################################
# demo_model.R
#
# FORMULA-STYLE model specification.  Instead of hand-writing an info_vector /
# info_matrix function (with factor contrast coding and the GLM variance weight),
# describe the model with:
#
#   link      : "identity" (linear), "logit" (logistic), "loglinear" (Poisson)
#   f         : main FACTOR effects            (factor indices)
#   x         : main CONTINUOUS linear effects (continuous indices)
#   fx        : factor x continuous interactions (two-digit codes or index pairs)
#   ff        : factor x factor interactions   (two-digit codes or index pairs)
#   xx        : continuous quadratic / interaction terms
#   intercept : include an intercept (default TRUE)
#   coding    : factor contrast coding, "zero-sum" (default) or "baseline"
#
# Factors and continuous covariates are declared in design_box exactly as
# before: a single integer L >= 2 is a factor with levels 1..L; c(lo, hi) is a
# continuous range.  The default criterion is D-optimality for all parameters.
###############################################################################

library(owea)

# ---------------------------------------------------------------------------
# (1) Linear model, 3-level factor (zero-sum) + one continuous covariate.
#     Equivalent to the hand-coded mix_vec in demo_factors.R -- but with no
#     model function to write.
# ---------------------------------------------------------------------------
res_lin <- optimal_design(design_box    = list(c(3), c(0, 1)),
                          step_sequence = c(0.1, 0.05),
                          link = "identity", f = 1, x = 1, p = 0)
print_result(res_lin, "Linear: intercept + 3-level factor + continuous")

# ---------------------------------------------------------------------------
# (2) Logistic model with a factor x continuous interaction.
#     design_box: one 2-level factor + one continuous covariate on [-1, 1].
#     Model: intercept + factor + continuous + factor:continuous.
# ---------------------------------------------------------------------------
theta <- c(0.5, -0.8, 1.2, 0.4)
res_log <- optimal_design(design_box    = list(dose = c(2), conc = c(-1, 1)),
                          step_sequence = c(0.2, 0.1, 0.05),
                          link = "logit", f = 1, x = 1, fx = c(11),
                          theta = theta, p = 1)              # A-optimality
print_result(res_log, "Logistic: factor + continuous + interaction (A-opt)")

# model_summary() labels each coefficient by its term (naming the covariates via
# a NAMED design_box, and showing the theta used).
model_summary(res_log)

# theta may be omitted for a GLM: it is then drawn from N(0,1) (with a warning)
# and returned in the result.
res_auto <- optimal_design(design_box    = list(c(2), c(-1, 1)),
                           step_sequence = c(0.2, 0.1),
                           link = "logit", f = 1, x = 1, fx = c(11), p = 0)
cat("\nauto-generated theta:", paste(round(res_auto$theta, 3), collapse = ", "),
    "\n")

# ---------------------------------------------------------------------------
# (3) Log-linear (Poisson) model, three factors + three continuous with
#     quadratic terms.  design_box = 3 factors (2,2,3 levels) + 3 continuous on
#     [0, 6].  Per-covariate step_sequence (a list) gives each continuous
#     covariate its own grid; theta is given explicitly for reproducibility.
# ---------------------------------------------------------------------------
theta_p <- c(0.2, -0.15, 0.1, 0.05, -0.05, 0.05,
             0.06, -0.03, 0.04, 0.02, -0.02, 0.03)   # length k = 12
res_pois <- optimal_design(
  design_box    = list(c(2), c(2), c(3), c(0, 6), c(0, 6), c(0, 6)),
  step_sequence = list(c(1, 1, 0.5), c(0.5, 0.5, 0.25), c(0.25, 0.25, 0.1)),
  link = "loglinear", f = c(1, 3), x = c(1, 2, 3),
  fx = c(11, 32), xx = c(11, 23), theta = theta_p, p = 0)
print_result(res_pois, "Log-linear: factors + continuous + interactions + quadratics")

# ---------------------------------------------------------------------------
# (4) Exact design from the same spec, and the information matrix of a design.
# ---------------------------------------------------------------------------
ex <- exact_design(n = 16, design_box = list(c(2), c(-1, 1)),
                   step_sequence = c(0.2, 0.1),
                   link = "logit", f = 1, x = 1, fx = c(11),
                   theta = theta, p = 0, seed = 1)
print(ex)

# design_information() also accepts the model spec directly (factor columns of
# the design are auto-detected, or pass factor_levels to be explicit).
M <- design_information(res_log$support, res_log$weights,
                        link = "logit", f = 1, x = 1, fx = c(11),
                        theta = theta, factor_levels = c(2, NA))
cat("\ndet(M)^(1/k) =", det(M)^(1 / nrow(M)), "\n")

# ---------------------------------------------------------------------------
# (5) Factor x factor interaction (ff).  Two factors (2 and 3 levels) with all
#     main effects and their interaction; identity link, D-optimality.
# ---------------------------------------------------------------------------
res_ff <- optimal_design(design_box    = list(A = c(2), B = c(3)),
                         step_sequence = numeric(0),      # all-factor design
                         link = "identity", f = c(1, 2), ff = c(12), p = 0)
print_result(res_ff, "Two factors with an A x B interaction")
model_summary(res_ff)

# ---------------------------------------------------------------------------
# (6) Multi-category responses.  Both use the info-matrix path (rank J-1 per
#     point) and need ncat = J and a theta (local designs).
# ---------------------------------------------------------------------------
# Baseline-category multinomial (J = 3): theta = (beta_1, beta_2) stacked.
res_mn <- optimal_design(design_box = list(dose = c(-1, 1)),
                         step_sequence = c(0.2, 0.1),
                         link = "multinomial", ncat = 3, x = 1,
                         theta = c(0.5, -0.8, -0.3, 1.1), p = 0)
print_result(res_mn, "Baseline-category multinomial logit (J = 3)")
model_summary(res_mn)

# Proportional-odds ordinal (J = 4): theta = (alpha_1 < alpha_2 < alpha_3, beta);
# the J-1 thresholds replace the intercept and must be increasing.
res_cm <- optimal_design(design_box = list(dose = c(-1, 1)),
                         step_sequence = c(0.2, 0.1),
                         link = "cumulative", ncat = 4, x = 1,
                         theta = c(-1.0, 0.0, 1.0, 0.8), p = 0)
print_result(res_cm, "Proportional-odds ordinal logit (J = 4)")
model_summary(res_cm)
