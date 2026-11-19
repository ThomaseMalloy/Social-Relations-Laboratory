# ==============================================================================
# Chapter 14, Step 0 (optional): Example triadic idiographic key person data
#
# Simulates 150 key people (KPs). Each KP nominates 3 friends (group g) and 3
# coworkers (group h), all with romantic partners. In each group the KP makes
# the 12 triadic assessments of Table 14.2 on two constructs:
#   listening    (X): how well the target listens to the target's partner
#   satisfaction (Y): the target's satisfaction with the target's partner
# Each assessment is "the KP's view of how <perceiver> assesses <target>'s
# relationship with <target>'s partner"; the KP row is the KP's own view.
#
# The data are generated from known values (listed at the bottom and written
# to output/true_values.txt) so Step 2 can be checked against the truth.
#
# Output (data/):
#   example_triadic_long.csv  one row per KP x group x perceiver x target
#   example_triadic_wide.csv  one row per KP x group, cells named
#                             <construct>_<perceiver>_<target> (e.g. satisfaction_KP_N1)
#
# Run from the chapter-14-triadic folder. Base R only.
# ==============================================================================

if (!dir.exists("data")) stop("Set the working directory to the chapter-14-triadic folder first.")
set.seed(14)

K   <- 150                       # key people
ids <- c("KP", "N1", "N2", "N3")  # focal people in each round robin
n   <- length(ids)

# ---- True values ------------------------------------------------------------
mu      <- c(friends = 6.5, coworkers = 5.5)   # satisfaction grand means (in-group > out-group)
mu_x    <- c(friends = 6.0, coworkers = 5.5)   # listening grand means
sd_kp_mean <- 0.6                               # SD of KP grand means (sigma^2_mu = .36)
var_a   <- c(friends = .50, coworkers = .25)   # perceiver variance (friends > coworkers)
var_b   <- c(friends = .40, coworkers = .40)   # target variance (equal)
var_g   <- c(friends = .60, coworkers = .25)   # relationship variance (friends > coworkers)
var_e   <- .10                                  # error variance
r_gg    <- .30                                  # dyadic (interpersonal) correlation of gamma
ego     <- 0.80                                 # egocentric bias: beta_KP - mean nominee beta
ego_sd  <- 0.60                                 # KP differences in egocentric bias
r_ego   <- .80                                  # consistency of KP bias across groups
r_akp   <- .70                                  # consistency of KP perceiver effect across groups
# Listening (X) and satisfaction (Y): correlations of the SRM effects
r_aa <- .40; r_bb <- .50; r_gxy <- .30

# Correlations of the KP's perceiver effects across groups (r_akp) and constructs (r_aa)
R_kp <- matrix(c(1,             r_akp,         r_aa,          r_aa * r_akp,
                 r_akp,         1,             r_aa * r_akp,  r_aa,
                 r_aa,          r_aa * r_akp,  1,             r_akp,
                 r_aa * r_akp,  r_aa,          r_akp,         1), 4, 4)
L_kp <- chol(R_kp)

rows <- list()
for (k in seq_len(K)) {
  m_k <- rnorm(1, 0, sd_kp_mean)                              # KP's general view
  z1 <- rnorm(1); z2 <- rnorm(1)
  ego_k <- ego + ego_sd * c(z1, r_ego * z1 + sqrt(1 - r_ego^2) * z2)
  # KP perceiver effects (z units): listening g, listening h, satisfaction g, satisfaction h
  akp_z <- as.vector(rnorm(4) %*% L_kp)
  for (h in 1:2) {
    grp <- names(mu)[h]
    sa <- sqrt(var_a[[grp]]); sb <- sqrt(var_b[[grp]]); sg <- sqrt(var_g[[grp]])
    # perceiver and target effects for X and Y (correlated r_aa, r_bb)
    ax <- rnorm(n); ay <- r_aa * ax + sqrt(1 - r_aa^2) * rnorm(n)
    ax[1] <- akp_z[h]; ay[1] <- akp_z[2 + h]
    bx <- rnorm(n); by <- r_bb * bx + sqrt(1 - r_bb^2) * rnorm(n)
    # relationship effects: g_ij with g_ji correlated r_gg; X with Y correlated r_gxy
    gx <- matrix(NA, n, n); gy <- matrix(NA, n, n)
    for (i in 1:(n - 1)) for (j in (i + 1):n) {
      u <- rnorm(1); v <- r_gg * u + sqrt(1 - r_gg^2) * rnorm(1)
      gx[i, j] <- u; gx[j, i] <- v
      gy[i, j] <- r_gxy * u + sqrt(1 - r_gxy^2) * rnorm(1)
      gy[j, i] <- r_gxy * v + sqrt(1 - r_gxy^2) * rnorm(1)
    }
    for (i in 1:n) for (j in 1:n) if (i != j) {
      bias_y <- if (j == 1) ego_k[h] else 0            # KP as target: egocentric bias
      bias_x <- if (j == 1) ego_k[h] / 2 else 0
      rows[[length(rows) + 1]] <- data.frame(
        kp.id = sprintf("KP%03d", k), group = grp, perceiver = ids[i], target = ids[j],
        listening    = mu_x[[grp]] + m_k + sa * ax[i] + sb * bx[j] + bias_x + sg * gx[i, j] + sqrt(var_e) * rnorm(1),
        satisfaction = mu[[grp]]   + m_k + sa * ay[i] + sb * by[j] + bias_y + sg * gy[i, j] + sqrt(var_e) * rnorm(1))
    }
  }
}
long <- do.call(rbind, rows)
long$listening <- round(long$listening, 3); long$satisfaction <- round(long$satisfaction, 3)
write.csv(long, "data/example_triadic_long.csv", row.names = FALSE)

# The same data in wide form (one row per KP x group), as survey software exports it
cells <- paste(long$perceiver, long$target, sep = "_")
wide <- unique(long[, c("kp.id", "group")])
for (v in c("listening", "satisfaction")) for (cc in unique(cells)) {
  sel <- cells == cc
  wide[[paste(v, cc, sep = "_")]] <- long[[v]][sel][match(paste(wide$kp.id, wide$group),
                                                         paste(long$kp.id, long$group)[sel])]
}
rownames(wide) <- NULL
write.csv(wide, "data/example_triadic_wide.csv", row.names = FALSE)

truth <- c(
  "True values used to simulate data/example_triadic_long.csv",
  sprintf("Key people: %d; focal people per round robin: %d; groups: friends (g), coworkers (h)", K, n),
  "",
  "Population values, and what the estimates in Step 2 should be close to:",
  sprintf("  Grand means, satisfaction: friends %.2f, coworkers %.2f; listening: friends %.2f, coworkers %.2f",
          mu[1], mu[2], mu_x[1], mu_x[2]),
  sprintf("    (plus a small upward shift from the KP egocentric bias: bias / n = +%.2f satisfaction, +%.2f listening)",
          ego / n, ego / 2 / n),
  sprintf("  Variance of KP grand means (sigma2_mu): %.2f (+ about %.3f from KP differences in bias)",
          sd_kp_mean^2, (ego_sd / n)^2),
  sprintf("  Perceiver variance: friends %.2f, coworkers %.2f", var_a[1], var_a[2]),
  sprintf("  Target variance: friends %.2f, coworkers %.2f; estimates also contain the KP bias, E(bias^2)/n:",
          var_b[1], var_b[2]),
  sprintf("    expected estimate %.2f for satisfaction and %.2f for listening in both groups",
          var_b[1] + (ego^2 + ego_sd^2) / n, var_b[1] + ((ego / 2)^2 + (ego_sd / 2)^2) / n),
  sprintf("  Relationship variance: friends %.2f, coworkers %.2f; with one indicator the estimate also",
          var_g[1], var_g[2]),
  sprintf("    contains error %.2f: expected %.2f and %.2f", var_e, var_g[1] + var_e, var_g[2] + var_e),
  sprintf("  Egocentric bias, beta_KP - mean nominee beta: %.2f satisfaction, %.2f listening (both groups)",
          ego, ego / 2),
  sprintf("  KP bias correlation across groups: %.2f; KP perceiver effect correlation across groups: %.2f",
          r_ego, r_akp),
  "    (observed correlations of the estimates are lower: they also contain sampling error)",
  sprintf("  Bivariate correlations of X and Y effects: perceiver %.2f, target %.2f, relationship %.2f",
          r_aa, r_bb, r_gxy),
  sprintf("    covariances: perceiver %.3f / %.3f, relationship %.3f / %.3f (friends / coworkers)",
          r_aa * var_a[1], r_aa * var_a[2], r_gxy * var_g[1], r_gxy * var_g[2]))
writeLines(truth, "output/true_values.txt")
cat("Wrote data/example_triadic_long.csv (", nrow(long), "rows ) and data/example_triadic_wide.csv (",
    nrow(wide), "rows )\n")
cat(truth, sep = "\n")
