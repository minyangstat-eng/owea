###############################################################################
# demo_factors.R
#
# FACTOR (categorical) covariates in optimal_design() / exact_design().
#
# A covariate in design_box is EITHER continuous -- a c(lo, hi) range stepped by
# step_sequence -- OR a factor: a single integer L >= 2 giving its number of
# levels (the levels are the integers 1..L).  A factor covariate is passed to
# info_vector / info_matrix as its RAW INTEGER LEVEL; any contrast coding is
# done INSIDE your own model function.  Reported support points show factor
# covariates as integer levels.
###############################################################################

library(owea)

# ---------------------------------------------------------------------------
# (1) The 2^k logistic screening model of Yang, Chen & Yang (2026), Example 3,
#     A-optimality.  Each of the k covariates is a TWO-LEVEL factor.
#
#     The model maps the factor LEVELS {1, 2} to the +/-1 coding it needs,
#     INSIDE the info_vector function.  Expressed as k two-level factors via
#     design_box, this reproduces the explicit 2^k +/-1 candidate set.
# ---------------------------------------------------------------------------

# information vector for the logistic model; x carries factor LEVELS in {1, 2}.
info_vec_ex3 <- function(x, theta) {
  xc  <- ifelse(x == 1, 1, -1)                  # level 1 -> +1, level 2 -> -1
  q   <- c(1, xc)                               # length k+1
  eta <- min(max(sum(q * theta), -500), 500)    # guard exp() overflow
  (exp(eta / 2) / (1 + exp(eta))) * q
}

k         <- 7
set.seed(2026)
theta_ex3 <- runif(k + 1, -3, 3)                # beta ~ U(-3, 3)^(k+1)

# --- baseline: explicit 2^k candidate set of +/-1 rows (no factor coding) --
info_vec_pm1 <- function(x, theta) {            # x already +/-1
  q   <- c(1, x)
  eta <- min(max(sum(q * theta), -500), 500)
  (exp(eta / 2) / (1 + exp(eta))) * q
}
X_ex3   <- as.matrix(expand.grid(rep(list(c(-1, 1)), k)))  # 2^k rows, k cols
res_set <- optimal_design(info_vector  = info_vec_pm1, theta = theta_ex3,
                          candidate_set = X_ex3, p = 1)     # A-optimality
print_result(res_set, "Example 3 -- explicit +/-1 candidate set")

# --- new way: k two-level FACTORS via design_box ---------------------------
#   design_box entry c(2) == a factor with 2 levels {1, 2}; the model encodes
#   the levels to +/-1 itself.  step_sequence = numeric(0): no continuous dims.
res_fac <- optimal_design(info_vector   = info_vec_ex3, theta = theta_ex3,
                          design_box    = rep(list(c(2)), k),
                          step_sequence = numeric(0), p = 1)
print_result(res_fac, "Example 3 -- k two-level factors (levels printed as integers)")

cat(sprintf("\n   -> criterion difference (set vs factor) = %.2e\n",
            abs(res_set$criterion - res_fac$criterion)))

# Model format #
# Zeor-sum#
res_fac <- optimal_design(f=c(1,2,3,4,5,6,7), link = "logit",  theta = theta_ex3,
                          design_box    = rep(list(c(2)), k),
                          step_sequence = numeric(0), p = 0)
print_result(res_fac, "Example 3 -- k two-level factors (levels printed as integers)")

#Baseline$
theta_ex4 = 2*theta_ex3
theta_ex4[1] = theta_ex3[1]-sum(theta_ex3[2:8])
res_fac0 <- optimal_design(f=c(1,2,3,4,5,6,7), link = "logit", coding="baseline", theta = theta_ex4,
                          design_box    = rep(list(c(2)), k),
                          step_sequence = numeric(0), p = 0)
print_result(res_fac0, "Example 3 -- k two-level factors (levels printed as integers)")
M0=design_information(f=c(1,2,3,4,5,6,7),link="logit",res_fac0$support,res_fac0$weights,coding="baseline", theta = theta_ex4)
M1=design_information(f=c(1,2,3,4,5,6,7),link="logit",res_fac$support,res_fac$weights,coding="baseline", theta = theta_ex4)

log(det(M0))-log(det(M1))

res_fac$support-res_fac0$support
res_fac$weights-res_fac0$weights

cat(sprintf("\n   -> criterion difference (set vs factor) = %.2e\n",
            abs(res_set$criterion - res_fac$criterion)))

# ---------------------------------------------------------------------------
# (2) Mixed factor + continuous covariates.
#     covariate 1: a 3-level factor;  covariate 2: continuous on [0, 1].
#     The model applies zero-sum contrasts to the 3-level factor ITSELF.
# ---------------------------------------------------------------------------
mix_vec <- function(x) {                        # x[1] in {1,2,3}; x[2] in [0,1]
  contr <- switch(as.character(x[1]),
                  "1" = c(1, 0), "2" = c(0, 1), "3" = c(-1, -1))  # zero-sum
  c(1, contr, x[2])                             # intercept + 2 contrasts + x2
}
theta_mix <- c(0.5, 1.0, -0.5, 2.0)             # (unused by D-optimality here)

res_mix <- optimal_design(info_vector   = mix_vec,
                          design_box    = list(c(3), c(0, 1)),
                          step_sequence = c(0.1, 0.05),
                          p             = 0)     # D-optimality
print_result(res_mix, "Mixed: 3-level factor + continuous [0,1]")

# exact design (n = 12) for the same mixed model
res_mix_exact <- exact_design(n = 12, info_vector = mix_vec,
                              design_box    = list(c(3), c(0, 1)),
                              step_sequence = c(0.1, 0.05),
                              p = 0, seed = 1)
print(res_mix_exact)

#Model version#
#Zero-sum#

res_mix <- optimal_design(f=c(1), x=c(1), link="identity",
                          design_box    = list(c(3), c(0, 1)),
                          step_sequence = c(0.1, 0.05),
                          p             = 0)     # D-optimality
print_result(res_mix, "Mixed: 3-level factor + continuous [0,1]")

#Baseline#

res_mix0 <- optimal_design(f=c(1), x=c(1), link="identity", coding = "baseline",
                          design_box    = list(c(3), c(0, 1)),
                          step_sequence = c(0.1, 0.05),
                          p             = 0)     # D-optimality
print_result(res_mix, "Mixed: 3-level factor + continuous [0,1]")

M0=design_information(f=c(1), x=c(1), link="identity",res_mix0$support,res_mix0$weights,coding="baseline")
M1=design_information(f=c(1), x=c(1), link="identity",res_mix$support,res_mix$weights,coding="baseline")

# ESD experiment (Lukemire et al.) -- logistic model with an ESD*Pulse interaction.
# x = c(A, B, ESD, Pulse, vol):
#   A, B, ESD, Pulse are FACTOR LEVELS in {1, 2}  (level 1 -> -1, level 2 -> +1)
#   vol is continuous in [25, 45]
# theta has length 7.
info_vec_esd <- function(x, theta) {
  code  <- c(-1, 1)                 # factor level 1 -> -1, level 2 -> +1
  A     <- code[x[1]]
  B     <- code[x[2]]
  ESD   <- code[x[3]]
  Pulse <- code[x[4]]
  vol   <- x[5]
  
  f <- c(1, A, B, ESD, Pulse, vol, ESD * Pulse)   # length 7 (out in the Julia code)
  xbeta <- sum(f * theta)
  xbeta <- min(max(xbeta, -500), 500)             # guard exp() overflow
  (exp(xbeta / 2) / (1 + exp(xbeta))) * f
}

theta0 <- c(-7.5, 1.50,-0.2,-0.15, 0.25, 0.35, 0.4)

res_fac <- optimal_design(info_vector   = info_vec_esd, theta = theta0,
                          design_box    = list(c(2),c(2),c(2),c(2),c(25,45)),
                          step_sequence = c(0.1, 0.01, 0.001), p = 0)
print_result(res_fac)

M1   <- infor_matrix(res_fac$support, res_fac$weights, info_vec_esd, theta = theta0)
det(M1)^(1/7)
sum(diag(solve(M1)))/7

#Model version#
#Zero-sum#
res_fac <- optimal_design(f=c(1,2,3,4),x=c(1),ff=c(34), link="logit", theta = theta0,
                          design_box    = list(c(2),c(2),c(2),c(2),c(25,45)),
                          step_sequence = c(0.1, 0.01, 0.001), p = 0)
print_result(res_fac)

M   <- infor_matrix(f=c(1,2,3,4),x=c(1),ff=c(34), link="logit", res_fac$support, res_fac$weights, theta = theta0)
