# Monte Carlo check: pooled half-block estimates are unbiased, and the tests hold their
# Type I rate when the true value is zero.
source("half_block_functions.R")
e <- new.env(); sys.source("0_make_example_data.R", envir = e); sim <- e$simulate_half_block
constructs <- list(warmth = c("warm", "kind"), competence = c("capable", "effective"))
truth <- c(warmth.actor = .10, warmth.partner = .60, warmth.relationship = .15, warmth.error = .20,
           competence.actor = .05, competence.partner = .40, competence.relationship = .10, competence.error = .20,
           warmth_competence.actor_cov = .50 * sqrt(.10 * .05),
           warmth_competence.partner_cov = -.40 * sqrt(.60 * .40),
           warmth_competence.relationship_cov = .20 * sqrt(.15 * .10))
set.seed(2026)
for (sz in list(c(4, 4), c(6, 6), c(3, 8))) {
  R <- t(replicate(500, { s <- hb_all(sim(K = 50, na = sz[1], np = sz[2]), constructs)$summary
                          setNames(s$estimate, s$term)[names(truth)] }))
  cat(sprintf("\nBlocks %d x %d (500 replications of 50 blocks)\n", sz[1], sz[2]))
  print(round(cbind(true = truth, mean_estimate = colMeans(R), mc_se = apply(R, 2, sd) / sqrt(nrow(R))), 4))
}
# Type I error: all actor variance and all covariances zero
p <- replicate(1000, { s <- hb_all(sim(K = 50, s2 = list(actor = c(0, 0), partner = c(.6, .4), relationship = c(.15, .1)),
                                        r = c(actor = 0, partner = 0, relationship = 0)), constructs)$summary
                       setNames(s$p, s$term)[c("warmth.actor", "warmth_competence.partner_cov")] })
cat("\nType I error rates (alpha = .05):\n"); print(rowMeans(p < .05))
# Own-twin test Type I error
o <- replicate(1000, { x <- sim(K = 50, own_effect = 0); x$partner <- sub("^P", "T", x$partner)
                       hb_own_compare(x, "warm")$p_one_tailed })
cat("Own-twin test Type I error:", mean(o < .05), "\n")
