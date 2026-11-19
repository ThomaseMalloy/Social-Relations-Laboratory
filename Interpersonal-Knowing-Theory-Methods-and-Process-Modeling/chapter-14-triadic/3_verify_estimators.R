# ==============================================================================
# Chapter 14, Step 3 (optional): Verify the estimators by Monte Carlo
#
# 1. Unbiasedness. Simulates many round robins with known SRM variances and
#    covariances and checks that the average estimates equal the true values,
#    for n = 4 (KP + 3 nominees, the minimal triadic round robin) and n = 6.
# 2. Type I error. Simulates data sets in which every null hypothesis of
#    Tables 14.6 and 14.7 is true and checks that each test rejects about 5%
#    of the time at alpha = .05.
#
# Writes output/verification.txt. Base R only. Takes about a minute.
# ==============================================================================

if (!dir.exists("data")) stop("Set the working directory to the chapter-14-triadic folder first.")
source("srm_triadic_functions.R")
set.seed(2014)
txt <- character(0); say <- function(...) { s <- sprintf(...); txt <<- c(txt, s); cat(s, "\n") }

# One round robin of size n with known components for X, and Y = .5 X + noise
sim_rr <- function(n, Sab, Sg) {
  ab <- matrix(rnorm(2 * n), n) %*% chol(Sab)
  G  <- matrix(NA, n, n)
  for (i in 1:(n - 1)) for (j in (i + 1):n) { z <- rnorm(2) %*% chol(Sg); G[i, j] <- z[1]; G[j, i] <- z[2] }
  X <- 5 + outer(ab[, 1], rep(1, n)) + outer(rep(1, n), ab[, 2]) + G
  Y <- .5 * X + .5 * (outer(rnorm(n), rep(1, n)) + outer(rep(1, n), rnorm(n)) + matrix(rnorm(n * n), n))
  diag(X) <- NA; diag(Y) <- NA
  list(X = X, Y = Y)
}

# ---- 1. Unbiasedness --------------------------------------------------------
Sab <- matrix(c(.30, .10, .10, .20), 2)   # var perceiver, cov, var target
Sg  <- matrix(c(.40, .15, .15, .40), 2)   # var relationship, dyadic cov
truth <- c(perceiver = .30, target = .20, relationship = .40, perceiver_target = .10, dyadic = .15,
           cov_aa = .15, cov_bb = .10, cov_gg_intra = .20, cov_gg_inter = .075)
say("1. Unbiasedness of the round robin estimators (Warner, Kenny, & Stoto, 1979)")
for (n in c(4, 6)) {
  R <- if (n == 4) 30000 else 10000
  est <- matrix(NA, R, 11)
  for (r in 1:R) {
    s <- sim_rr(n, Sab, Sg); f <- srm_one(s$X)
    est[r, ] <- c(f$var, srm_cov(s$X, s$Y), f$mean, grand_mean_sampling_var(f$var, n))
  }
  m  <- colMeans(est); se <- apply(est, 2, sd) / sqrt(R)
  say("\n  n = %d (%d round robins)        true   mean estimate   (MC SE)", n, R)
  lab <- c("perceiver variance", "target variance", "relationship variance", "perceiver-target cov",
           "dyadic (reciprocity) cov", "Cov(alpha_x, alpha_y)", "Cov(beta_x, beta_y)",
           "Cov(gamma_x, gamma_y) intra", "Cov(gamma_ij_x, gamma_ji_y)")
  for (k in 1:9) say("  %-30s %6.3f   %6.3f   (%.3f)", lab[k], truth[k], m[k], se[k])
  tv <- ((n - 1) * (.30 + .20 + 2 * .10) + .40 + .15) / (n * (n - 1))
  say("  %-30s %6.3f   %6.3f   (empirical variance of grand means %.3f)", "grand mean sampling variance", tv, m[11], var(est[, 10]))
  ok <- abs(m[1:9] - truth) < 4 * se[1:9] + .003
  if (!all(ok)) stop("Estimator check failed for n = ", n, ": ", paste(lab[!ok], collapse = ", "))
}
say("\n  All estimates are within Monte Carlo error of the true values.")

# ---- 2. Type I error under the null hypotheses ------------------------------
reps <- 1000; K <- 80; n <- 4
say("\n2. Type I error rates at alpha = .05 (%d simulated data sets of %d KPs, n = %d, two groups, all nulls true)", reps, K, n)
rej <- matrix(NA, reps, 8, dimnames = list(NULL, c("alpha_KP = mean alpha_N", "beta_KP = mean beta_N",
       "mu_g = mu_h", "sigma2_alpha_g = sigma2_alpha_h", "sigma2_gamma_g = sigma2_gamma_h",
       "r(beta_KP_g, beta_KP_h) = 0", "sigma2_total_g = sigma2_total_h",
       "sigma2_mu = 0 (t test of d_k)")))
for (r in 1:reps) {
  st <- array(NA, c(K, 2, 6))   # KP x group x (alpha_KP, mean alpha_N, beta_KP, mean beta_N, mean, var_a)
  g  <- matrix(NA, K, 2)        # relationship variance
  tv <- matrix(NA, K, 2)        # total variance (perceiver + target + relationship)
  vm <- matrix(NA, K, 2)        # sampling variance of the grand mean
  for (k in 1:K) for (h in 1:2) {
    f <- srm_one(sim_rr(n, Sab, Sg)$X)
    st[k, h, ] <- c(f$a[1], mean(f$a[-1]), f$b[1], mean(f$b[-1]), f$mean, f$var[["perceiver"]])
    g[k, h]    <- f$var[["relationship"]]
    tv[k, h]   <- f$var[["perceiver"]] + f$var[["target"]] + f$var[["relationship"]]
    vm[k, h]   <- grand_mean_sampling_var(f$var, n)
  }
  rej[r, ] <- c(t.test(st[, 1, 1], st[, 1, 2], paired = TRUE)$p.value,
                t.test(st[, 1, 3], st[, 1, 4], paired = TRUE)$p.value,
                t.test(st[, 1, 5], st[, 2, 5], paired = TRUE)$p.value,
                t.test(st[, 1, 6], st[, 2, 6], paired = TRUE)$p.value,
                t.test(g[, 1], g[, 2], paired = TRUE)$p.value,
                cor.test(st[, 1, 3], st[, 2, 3])$p.value,
                t.test(tv[, 1], tv[, 2], paired = TRUE)$p.value,
                t.test(K / (K - 1) * (st[, 1, 5] - mean(st[, 1, 5]))^2 - vm[, 1],
                       alternative = "greater")$p.value) < .05
}
rate <- colMeans(rej); mcse <- sqrt(.05 * .95 / reps)
for (k in seq_along(rate)) say("  %-34s rejection rate %.3f", colnames(rej)[k], rate[k])
say("  (expected .050; Monte Carlo 95%% range about %.3f to %.3f)", .05 - 2 * mcse, .05 + 2 * mcse)
writeLines(txt, "output/verification.txt")
cat("\nWrote output/verification.txt\n")
