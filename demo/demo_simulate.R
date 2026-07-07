###############################################################################
# demo_simulate.R
#
# simulate_design(): given a design (support points + integer run counts), a
# formula-style model (link + f/x/fx/ff/xx, as in optimal_design) and the TRUE
# parameters theta, simulate the responses and fit the model back (MLE). With
# nsim > 1 it repeats the simulate-and-fit and compares the empirical covariance
# of the estimates with the design-based large-sample covariance M^{-1}/N.
###############################################################################

library(owea)

# ===========================================================================
# (A) The exact_design workflow: design -> exact allocation -> simulate + fit
# ===========================================================================
# Find the (approximate then integer) locally D-optimal design for a logistic
# model logit(p) = b0 + b1*dose on dose in [-1, 1], at theta = (0.3, 1.5).
ed <- exact_design(n = 40, design_box = list(dose = c(-1, 1)),
                   step_sequence = c(0.5, 0.25, 0.1),
                   link = "logit", x = 1, theta = c(0.3, 1.5), seed = 1)
print(ed)                       # integer counts n_i summing to 40

# Feed the exact_design OBJECT straight into simulate_design(): its support and
# counts are taken from the object.  You still supply the model (link/x/theta).
sim <- simulate_design(ed, link = "logit", x = 1,
                       design_box = list(dose = c(-1, 1)),
                       theta = c(0.3, 1.5), seed = 2)

cat("\n--- what simulate_design() returns (single run) ---\n")
cat("names(sim):", paste(names(sim), collapse = ", "), "\n")
cat("N (total runs):", sim$N, "\n")
cat("first rows of the simulated data:\n"); print(head(sim$data, 4))
cat("theta_hat:\n"); print(round(sim$theta_hat, 4))
cat("se:\n");        print(round(sim$se, 4))
cat("logLik =", round(sim$loglik, 3), "  converged =", sim$converged, "\n")
cat("vcov (2x2):\n"); print(round(sim$vcov, 5))

# ===========================================================================
# (B) One example per LINK (single simulated dataset; theta vs estimate)
#     A design here is just any support matrix + counts -- it need not be
#     optimal.  design_box tells simulate_design which columns are factors.
# ===========================================================================
sup <- matrix(c(-1, -0.5, 0, 0.5, 1), ncol = 1)   # 5 support points, 1 covariate
cnt <- rep(200, 5)                                # 200 runs each
box <- list(dose = c(-1, 1))
show <- function(r, label) {
  cat(sprintf("\n--- %s ---\n", label))
  print(round(rbind(true = r$theta, est = r$theta_hat, se = r$se), 3))
  cat("logLik =", round(r$loglik, 2), "  converged =", r$converged, "\n")
}

# 1. identity  (linear / Gaussian):  y ~ N(b0 + b1*dose, sigma^2)
r_id <- simulate_design(sup, cnt, link = "identity", x = 1, design_box = box,
                        theta = c(2, -3), sigma = 0.5, seed = 10)
show(r_id, "identity (Gaussian), sigma = 0.5")
cat("sigma_hat =", round(r_id$sigma_hat, 3), "\n")

# 2. logit  (logistic):  y ~ Bernoulli(plogis(b0 + b1*dose))
r_lg <- simulate_design(sup, cnt, link = "logit", x = 1, design_box = box,
                        theta = c(0.3, 1.5), seed = 11)
show(r_lg, "logit (logistic)")
cat("response 0/1 counts:", paste(table(r_lg$data$y), collapse = " / "), "\n")

# 3. loglinear  (Poisson log link):  y ~ Poisson(exp(b0 + b1*dose))
r_ll <- simulate_design(sup, cnt, link = "loglinear", x = 1, design_box = box,
                        theta = c(1.0, 0.7), seed = 12)
show(r_ll, "loglinear (Poisson)")
cat("mean count =", round(mean(r_ll$data$y), 2), "\n")

# 4. multinomial  (baseline-category, J = 3):
#    theta = (beta_1, beta_2) stacked; baseline = category 3.
r_mn <- simulate_design(sup, cnt, link = "multinomial", ncat = 3, x = 1,
                        design_box = box, theta = c(0.3, 1.0, -0.4, 0.8),
                        seed = 13)
show(r_mn, "multinomial (J = 3)")
cat("category 1/2/3 counts:", paste(table(r_mn$data$y), collapse = " / "), "\n")

# 5. cumulative  (proportional-odds ordinal, J = 3):
#    theta = (alpha_1 < alpha_2, beta); thresholds replace the intercept.
r_cm <- simulate_design(sup, cnt, link = "cumulative", ncat = 3, x = 1,
                        design_box = box, theta = c(-0.8, 0.6, 1.2), seed = 14)
show(r_cm, "cumulative / ordinal (J = 3)")
cat("ordered 1/2/3 counts:", paste(table(r_cm$data$y), collapse = " / "), "\n")

# ===========================================================================
# (C) Monte Carlo: check the design delivers its predicted precision.
#     nsim > 1 returns estimate summaries instead of a single fit.
# ===========================================================================
mc <- simulate_design(sup, cnt, link = "logit", x = 1, design_box = box,
                      theta = c(0.3, 1.5), nsim = 500, seed = 20)
cat("\n--- Monte Carlo (nsim = 500), what it returns ---\n")
cat("names(mc):", paste(names(mc), collapse = ", "), "\n")
cat("n_converged:", mc$n_converged, "of", mc$nsim, "\n")
print(round(rbind(true          = mc$theta,
                  mean_estimate = mc$theta_hat_mean,
                  bias          = mc$bias,
                  se_empirical  = mc$se_empirical,
                  se_design     = mc$se_design), 4))
cat("empirical vs design covariance:\n")
print(round(mc$cov_empirical, 5)); print(round(mc$cov_design, 5))
