# ==============================================================================
# Chapter 14, Step 2: Tests of the triadic hypotheses (Tables 14.6, 14.7, 14.9)
#
# Reads the second-stage data sets written by Step 1 and tests:
#   Table 14.6 (within each group): variance of KP grand means; perceiver,
#     target, and relationship variances; KP vs average nominee perceiver and
#     target effects (including the directional egocentric bias test)
#   Table 14.7 (between groups): grand means; total variance (perceived
#     heterogeneity of the group); perceiver, target, and relationship variances; KP perceiver and target effects; consistency
#     (correlation) of KP effects across groups
#   Table 14.9 (bivariate): perceiver-perceiver, target-target, and
#     intrapersonal relationship-relationship covariances of two constructs
#
# The unit of analysis is the key person. Each KP contributes one estimate of
# each parameter per group; tests are t tests across KPs (the idiographic
# approach used throughout the book and by TripleR for multiple groups).
#
# Variance of the grand mean (sigma^2_mu = 0): each KP has only one grand mean
# per group, so a t test of the grand means tests their mean, not their
# variance. The test asks whether KP grand means vary more than their sampling
# error allows. Three versions:
#   (a) a multilevel model of the raw scores, KP random intercept with
#       perceiver and target nested in KP, likelihood ratio test (p halved
#       because the variance is tested at its boundary); needs lme4
#   (b) a moment estimate: observed variance of KP grand means minus the
#       average sampling variance implied by each KP's SRM components
#   (c) a one-sample t test (df = K - 1) of d_k = K/(K-1)(M_k - mean M)^2 minus
#       KP k's sampling variance; the mean of d_k is the estimate in (b)
#
# Output (output/): triadic_results.txt, triadic_results.csv
# Run from the chapter-14-triadic folder after Step 1. Base R (+ lme4 for (a)).
# ==============================================================================

if (!dir.exists("data")) stop("Set the working directory to the chapter-14-triadic folder first.")

# ---- Settings ---------------------------------------------------------------
raw_file   <- "data/example_triadic_long.csv"   # long data used in Step 1 (for test a)
constructs <- c("listening", "satisfaction")
groups     <- c("friends", "coworkers")         # g and h, in that order
bivariate  <- c("listening", "satisfaction")    # X and Y for Table 14.9
# -----------------------------------------------------------------------------

kp <- read.csv("data/triadic_kp_level.csv", stringsAsFactors = FALSE)
w  <- read.csv("data/triadic_kp_wide.csv",  stringsAsFactors = FALSE)
out <- list(); txt <- character(0)
say <- function(...) { s <- sprintf(...); txt <<- c(txt, s); cat(s, "\n") }
add <- function(table, group, construct, parameter, hypothesis, est, test, stat, df, p, note = "") {
  out[[length(out) + 1]] <<- data.frame(table, group, construct, parameter, hypothesis,
    estimate = est, test, statistic = stat, df = df, p = p, note, stringsAsFactors = FALSE)
}
t1 <- function(x, alt = "two.sided") { tt <- t.test(x, alternative = alt); c(mean(x), tt$statistic, tt$parameter, tt$p.value) }
tp <- function(x, y, alt = "two.sided") { tt <- t.test(x, y, paired = TRUE, alternative = alt); c(mean(x - y), tt$statistic, tt$parameter, tt$p.value) }
fmt_p <- function(p) ifelse(p < .001, "< .001", sprintf("= %.3f", p))
line <- function(label, r, extra = "") say("  %-46s est = %7.3f  t(%d) = %6.2f  p %s %s", label, r[1], as.integer(r[3]), r[2], fmt_p(r[4]), extra)

have_lme4 <- requireNamespace("lme4", quietly = TRUE)
raw <- if (file.exists(raw_file)) read.csv(raw_file, stringsAsFactors = FALSE) else NULL

# ============================ Table 14.6: within =============================
say("TABLE 14.6  Within-group triadic hypotheses (unit = key person; N = %d KPs)", length(unique(kp$kp.id)))
for (g in groups) for (v in constructs) {
  k <- kp[kp$group == g, ]
  say("\n[%s | %s]", g, v)
  # Grand mean variance
  m  <- k[[paste0(v, ".mean")]]
  vm <- var(m) - mean(k[[paste0(v, ".var.mean")]])
  say("  Grand mean: mean of KP grand means = %.3f; SD = %.3f", mean(m), sd(m))
  say("  sigma2_mu (moment estimate: var of KP means - mean sampling var) = %.3f", vm)
  add("14.6", g, v, "grand mean", "mean of KP grand means", mean(m), "descriptive", NA, NA, NA)
  add("14.6", g, v, "sigma2_mu", "sigma2_mu = 0", vm, "moment estimate", NA, NA, NA)
  # One-sample t test: d_k = K/(K-1) (M_k - mean M)^2 - sampling var of M_k; mean d_k estimates sigma2_mu
  K  <- length(m)
  dk <- K / (K - 1) * (m - mean(m))^2 - k[[paste0(v, ".var.mean")]]
  r  <- t1(dk, "greater")
  line("sigma2_mu, t test of d_k (one-sided, > 0)", r)
  add("14.6", g, v, "sigma2_mu", "sigma2_mu = 0", r[1], "one-sample t of d_k, one-sided", r[2], r[3], r[4],
      "d_k = K/(K-1)(M_k - mean M)^2 - sampling variance of M_k")
  if (have_lme4 && !is.null(raw)) {
    r <- raw[raw$group == g, ]; r$y <- r[[v]]
    f1 <- lme4::lmer(y ~ 1 + (1 | kp.id) + (1 | kp.id:perceiver) + (1 | kp.id:target), data = r, REML = FALSE)
    f0 <- lme4::lmer(y ~ 1 +                 (1 | kp.id:perceiver) + (1 | kp.id:target), data = r, REML = FALSE)
    lr <- as.numeric(2 * (logLik(f1) - logLik(f0))); p <- 0.5 * pchisq(lr, 1, lower.tail = FALSE)
    vc <- as.data.frame(lme4::VarCorr(f1)); vk <- vc$vcov[vc$grp == "kp.id"]
    say("  sigma2_mu (multilevel, KP random intercept) = %.3f  LR chi2(1) = %.2f  p %s (halved)", vk, lr, fmt_p(p))
    add("14.6", g, v, "sigma2_mu", "sigma2_mu = 0", vk, "multilevel LR test (p halved)", lr, 1, p)
  }
  # Perceiver, target, relationship variances
  for (s in list(c("var.a", "sigma2_alpha", "Perceiver variance"),
                 c("var.b", "sigma2_beta",  "Target variance"),
                 c("var.g", "sigma2_gamma", "Relationship variance"))) {
    r <- t1(k[[paste0(v, ".", s[1])]], "greater")
    line(paste0(s[3], " (one-sided, > 0)"), r)
    add("14.6", g, v, s[2], paste(s[2], "= 0"), r[1], "one-sample t, one-sided", r[2], r[3], r[4])
  }
  # Relative variances (percentage of the summed components)
  tot <- sum(sapply(c("var.a", "var.b", "var.g"), function(s) mean(k[[paste0(v, ".", s)]])))
  say("  Relative variance: perceiver %.0f%%, target %.0f%%, relationship (incl. error) %.0f%%",
      100 * mean(k[[paste0(v, ".var.a")]]) / tot, 100 * mean(k[[paste0(v, ".var.b")]]) / tot,
      100 * mean(k[[paste0(v, ".var.g")]]) / tot)
  # KP vs average nominee
  r <- tp(k[[paste0(v, ".a.kp")]], k[[paste0(v, ".a.nom")]])
  line("alpha_KP = mean nominee alpha (paired)", r)
  add("14.6", g, v, "alpha_KP - mean alpha_N", "alpha_KP = mean alpha_N", r[1], "paired t", r[2], r[3], r[4])
  r <- tp(k[[paste0(v, ".b.kp")]], k[[paste0(v, ".b.nom")]])
  line("beta_KP = mean nominee beta (paired)", r)
  add("14.6", g, v, "beta_KP - mean beta_N", "beta_KP = mean beta_N", r[1], "paired t", r[2], r[3], r[4])
  r <- tp(k[[paste0(v, ".b.kp")]], k[[paste0(v, ".b.nom")]], "greater")
  line("Egocentric bias beta_KP > mean beta_N", r, "(one-sided)")
  add("14.6", g, v, "beta_KP - mean beta_N", "beta_KP > mean beta_N (egocentric bias)", r[1], "paired t, one-sided", r[2], r[3], r[4])
}

# ============================ Table 14.7: between ============================
gg <- groups[1]; hh <- groups[2]
say("\nTABLE 14.7  Between-group triadic hypotheses (%s = g, %s = h; paired over KPs)", gg, hh)
col <- function(v, s, grp) w[[paste0(v, ".", s, "_", grp)]]
for (v in constructs) {
  say("\n[%s]", v)
  for (s in list(c("mean", "mu_g = mu_h", "Grand means"),
                 c("var.total", "sigma2_total_g = sigma2_total_h", "Total variances (heterogeneity)"),
                 c("var.a", "sigma2_alpha_g = sigma2_alpha_h", "Perceiver variances"),
                 c("var.b", "sigma2_beta_g = sigma2_beta_h", "Target variances"),
                 c("var.g", "sigma2_gamma_g = sigma2_gamma_h", "Relationship variances"),
                 c("a.kp", "alpha_KP_g = alpha_KP_h", "KP perceiver effects"),
                 c("b.kp", "beta_KP_g = beta_KP_h", "KP target effects"))) {
    r <- tp(col(v, s[1], gg), col(v, s[1], hh))
    line(paste0(s[3], " (g - h)"), r)
    add("14.7", paste(gg, "vs", hh), v, s[1], s[2], r[1], "paired t", r[2], r[3], r[4])
  }
  for (s in list(c("a.kp", "KP perceiver effect"), c("b.kp", "KP target effect"))) {
    ct <- cor.test(col(v, s[1], gg), col(v, s[1], hh))
    say("  Consistency of %-24s r = %6.3f  t(%d) = %6.2f  p %s", paste0(s[2], ":"), ct$estimate,
        as.integer(ct$parameter), ct$statistic, fmt_p(ct$p.value))
    add("14.7", paste(gg, "with", hh), v, paste("r", s[1]), "consistency across groups (r = 0)",
        unname(ct$estimate), "Pearson r", unname(ct$statistic), unname(ct$parameter), ct$p.value)
  }
}

# ============================ Table 14.9: bivariate ==========================
X <- bivariate[1]; Y <- bivariate[2]; tag <- paste0(X, ".", Y)
say("\nTABLE 14.9  Bivariate associations: X = %s, Y = %s (one-sided, > 0)", X, Y)
for (g in groups) {
  k <- kp[kp$group == g, ]
  say("\n[%s]", g)
  for (s in list(c("cov.aa", "var.a", "Cov(alpha_x, alpha_y)"),
                 c("cov.bb", "var.b", "Cov(beta_x, beta_y)"),
                 c("cov.gg.intra", "var.g", "Cov(gamma_ij_x, gamma_ij_y)"))) {
    r  <- t1(k[[paste0(tag, ".", s[1])]], "greater")
    vx <- mean(k[[paste0(X, ".", s[2])]]); vy <- mean(k[[paste0(Y, ".", s[2])]])
    rr <- if (vx > 0 && vy > 0) r[1] / sqrt(vx * vy) else NA
    line(s[3], r, sprintf("(r = %.2f)", rr))
    add("14.9", g, tag, s[3], paste(s[3], "> 0"), r[1], "one-sample t, one-sided", r[2], r[3], r[4],
        sprintf("standardized r = %.3f (mean cov / sqrt(mean var x * mean var y))", rr))
  }
}

res <- do.call(rbind, out)
write.csv(res, "output/triadic_results.csv", row.names = FALSE)
writeLines(txt, "output/triadic_results.txt")
cat("\nWrote output/triadic_results.txt and output/triadic_results.csv\n")
