# srm_functions.R
# Round-robin SRM estimation for ONE idiographic round robin (one key person in one
# family), following the TripleR conventions (Schönbrodt, Back, & Schmukle, 2012;
# Kenny, 1994; Warner, Kenny, & Stoto, 1979):
#   * actor, partner, relationship effects from the round-robin ANOVA formulas
#   * variance and covariance components with the standard bias corrections
#   * latent constructs with two indicators: stable components are cross-indicator
#     covariances; error variance = single-indicator relationship variance minus the
#     latent (cross-indicator) relationship variance
#   * standardized variances = component / (actor + partner + relationship + error);
#     covariances standardized as correlations (NA if a variance is <= 0; values outside
#     [-1, 1] are reported as they are, so out-of-bounds estimates stay visible)
# Used when TripleR is not installed, and to check the TripleR results.

srm_effects <- function(X) {
  # X: n x n matrix (rows = actors, columns = partners, NA diagonal)
  n <- nrow(X); diag(X) <- NA
  rm <- rowMeans(X, na.rm = TRUE); cm <- colMeans(X, na.rm = TRUE); M <- mean(X, na.rm = TRUE)
  k1 <- (n - 1)^2 / (n * (n - 2)); k2 <- (n - 1) / (n * (n - 2)); k3 <- (n - 1) / (n - 2)
  a <- k1 * rm + k2 * cm - k3 * M
  b <- k1 * cm + k2 * rm - k3 * M
  g <- X - outer(a, rep(1, n)) - outer(rep(1, n), b) - M
  diag(g) <- NA
  list(mean = M, a = a, b = b, g = g, n = n)
}

# bivariate components between two round robins measured on the same people
srm_bivariate <- function(ex, ey) {
  n <- ex$n; d <- (n - 1) * (n - 2) / 2
  ij <- which(upper.tri(ex$g), arr.ind = TRUE)
  gx1 <- ex$g[ij]; gx2 <- ex$g[ij[, 2:1]]; gy1 <- ey$g[ij]; gy2 <- ey$g[ij[, 2:1]]
  MSs <- sum((gx1 + gx2) * (gy1 + gy2) / 2) / (d - 1)
  MSd <- sum((gx1 - gx2) * (gy1 - gy2) / 2) / d
  intra <- (MSs + MSd) / 2            # relationship covariance, same direction (ij, ij)
  inter <- (MSs - MSd) / 2            # relationship covariance, reversed (ij, ji)
  c1 <- (n - 1) / (n * (n - 2)); c2 <- 1 / (n * (n - 2))
  AA <- sum(ex$a * ey$a) / (n - 1) - c1 * intra - c2 * inter
  BB <- sum(ex$b * ey$b) / (n - 1) - c1 * intra - c2 * inter
  AB <- sum(ex$a * ey$b) / (n - 1) - c1 * inter - c2 * intra   # actor x, partner y
  BA <- sum(ex$b * ey$a) / (n - 1) - c1 * inter - c2 * intra   # partner x, actor y
  c(AA = AA, BB = BB, AB = AB, BA = BA, RR_intra = intra, RR_inter = inter)
}

safe_r <- function(cov, v1, v2) ifelse(v1 > 0 & v2 > 0, cov / sqrt(v1 * v2), NA_real_)

# latent construct from two indicator matrices X1, X2
srm_latent <- function(X1, X2) {
  e1 <- srm_effects(X1); e2 <- srm_effects(X2)
  b12 <- srm_bivariate(e1, e2); b11 <- srm_bivariate(e1, e1); b22 <- srm_bivariate(e2, e2)
  vA <- b12[["AA"]]; vP <- b12[["BB"]]; vR <- b12[["RR_intra"]]
  vE <- (b11[["RR_intra"]] + b22[["RR_intra"]]) / 2 - vR
  cAP <- (b12[["AB"]] + b12[["BA"]]) / 2
  cRR <- b12[["RR_inter"]]
  tot <- vA + vP + vR + vE
  em <- srm_effects((X1 + X2) / 2)         # effects of the construct (mean of indicators)
  list(mean = em$mean, a = em$a, b = em$b, g = em$g,
       comp = c(varA = vA, varP = vP, varR = vR, varE = vE, covAP = cAP, covR = cRR,
                stdA = vA / tot, stdP = vP / tot, stdR = vR / tot, stdE = vE / tot,
                genrec = safe_r(cAP, vA, vP), dyrec = safe_r(cRR, vR, vR)),
       e1 = e1, e2 = e2)
}

# latent bivariate (construct x construct), e.g. positivity with negativity
srm_latent_bivariate <- function(LX, LY) {
  # average the four cross-indicator combinations (x1,y1), (x1,y2), (x2,y1), (x2,y2)
  B <- rowMeans(sapply(list(list(LX$e1, LY$e1), list(LX$e1, LY$e2),
                            list(LX$e2, LY$e1), list(LX$e2, LY$e2)),
                       function(p) srm_bivariate(p[[1]], p[[2]])))
  cx <- LX$comp; cy <- LY$comp
  c(covAA = B[["AA"]], covPP = B[["BB"]], covAP = B[["AB"]], covPA = B[["BA"]],
    covRintra = B[["RR_intra"]], covRinter = B[["RR_inter"]],
    rAA = safe_r(B[["AA"]], cx[["varA"]], cy[["varA"]]),
    rPP = safe_r(B[["BB"]], cx[["varP"]], cy[["varP"]]),
    rAP = safe_r(B[["AB"]], cx[["varA"]], cy[["varP"]]),
    rPA = safe_r(B[["BA"]], cx[["varP"]], cy[["varA"]]),
    rRintra = safe_r(B[["RR_intra"]], cx[["varR"]], cy[["varR"]]),
    rRinter = safe_r(B[["RR_inter"]], cx[["varR"]], cy[["varR"]]))
}

# first-stage SRM for every family x key person in raw long data -> long parameter file
# (needs param_functions.R and sim_functions.R for rr_matrix / srm_to_long)
srm_mirr_builtin <- function(raw) {
  out <- list()
  for (f in sort(unique(raw$family))) for (k in ROLES) {
    dk <- raw[raw$family == f & raw$source == k, ]
    if (nrow(dk) == 0) next
    est <- list()
    for (c in CONSTRUCTS) est[[c]] <- srm_latent(rr_matrix(dk, paste0(c, 1)), rr_matrix(dk, paste0(c, 2)))
    est$biv <- srm_latent_bivariate(est$POS, est$NEG)
    out[[length(out) + 1]] <- srm_to_long(f, k, est)
  }
  do.call(rbind, out)
}
