# Monte Carlo check of the ABD estimators: bias and Type I error.
source("abd_functions.R")
simulate_abd <- function(K = 50, nA = 3, nB = 3, vA = c(.20, .10), vP = c(.15, .25), rAP = c(.30, .40),
                         vG = .30, rGG = .25, err = .08) {
  # vA, vP: actor and partner variances for direction A->B (first) and B->A (second)
  # rAP: actor-partner correlation for A members (first) and B members (second)
  # vG, rGG: relationship variance and dyadic reciprocity correlation; err: error per indicator
  bv <- function(n, v1, v2, r) { z <- matrix(rnorm(2 * n), n); cbind(z[, 1] * sqrt(v1), (r * z[, 1] + sqrt(1 - r^2) * z[, 2]) * sqrt(v2)) }
  out <- vector("list", K)
  for (k in 1:K) {
    EA <- bv(nA, vA[1], vP[2], rAP[1])   # A: actor (A->B), partner (B->A)
    EB <- bv(nB, vA[2], vP[1], rAP[2])   # B: actor (B->A), partner (A->B)
    g <- expand.grid(i = 1:nA, j = 1:nB); G <- bv(nrow(g), vG, vG, rGG)
    mk <- function(v) rnorm(length(v), v, sqrt(err))
    ab <- 3 + EA[g$i, 1] + EB[g$j, 2] + G[, 1]; ba <- 3 + EB[g$j, 1] + EA[g$i, 2] + G[, 2]
    out[[k]] <- rbind(
      data.frame(block = k, actor = paste0("A", g$i), actor_group = "A", partner = paste0("B", g$j), w1 = mk(ab), w2 = mk(ab)),
      data.frame(block = k, actor = paste0("B", g$j), actor_group = "B", partner = paste0("A", g$i), w1 = mk(ba), w2 = mk(ba)))
  }
  do.call(rbind, out)
}
truth <- c(A_B.c.actor = .20, A_B.c.partner = .15, A_B.c.relationship = .30,
           B_A.c.actor = .10, B_A.c.partner = .25, B_A.c.relationship = .30,
           recip.c.A = .30 * sqrt(.20 * .25), recip.c.B = .40 * sqrt(.10 * .15), recip.c.relationship = .25 * .30)
set.seed(2026)
for (sz in list(c(3, 3), c(4, 6))) {
  R <- t(replicate(300, { s <- abd_all(simulate_abd(nA = sz[1], nB = sz[2]), list(c = c("w1", "w2")), c("A", "B"))$summary
                          setNames(s$estimate, s$term)[names(truth)] }))
  cat(sprintf("\nBlocks of %d + %d (300 replications of 50 blocks)\n", sz[1], sz[2]))
  print(round(cbind(true = truth, mean_estimate = colMeans(R), mc_se = apply(R, 2, sd) / sqrt(nrow(R))), 4))
}
p <- replicate(1000, { s <- abd_all(simulate_abd(vA = c(0, .1), rAP = c(0, 0), rGG = 0), list(c = c("w1", "w2")), c("A", "B"))$summary
                       setNames(s$p, s$term)[c("A_B.c.actor", "recip.c.A", "recip.c.relationship")] })
cat("\nType I error rates (alpha = .05) when the true value is zero:\n"); print(rowMeans(p < .05))
