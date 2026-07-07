# Shared models used across tests.

theta_biexp <- c(1.0, 1.0, 1.0, 2.0)

biexp_vec <- function(x, th) {
  xx <- x[1]
  c(            exp(-th[2] * xx),
    -th[1] * xx * exp(-th[2] * xx),
                exp(-th[4] * xx),
    -th[3] * xx * exp(-th[4] * xx))
}
biexp_info <- function(x) tcrossprod(biexp_vec(x, theta_biexp))

beta_logit <- c(1.0, -0.5, 0.5, 1.0)
logistic_info <- function(x) {
  q <- c(1.0, x[1], x[2], x[3])
  e <- exp(sum(q * beta_logit))
  (e / (1 + e)^2) * tcrossprod(q)
}
logistic_vec <- function(x, th) {
  q <- c(1.0, x[1], x[2], x[3])
  e <- exp(sum(q * th))
  sqrt(e / (1 + e)^2) * q          # f f' = (e/(1+e)^2) q q'
}
