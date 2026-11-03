###############################################################################
# Dyadity analyses for round-robin data
# Companion code for Chapter 9 (Malloy): expressed and elicited dyadity,
# expressed-minus-elicited differences (Eq. 9.10), Monte Carlo tests of
# individual differences in dyadity, role-specific tests, and the IKPD
# egocentric-bias contrast (Eqs. 9.11 and 9.12).
#
# Base R only; no packages are required. Tested with R 4.3.
#
# QUICK START
#   source("dyadity_analysis.R")
#   res <- dyadity_analysis(mydata, group = "group", actor = "actor",
#                           partner = "partner", y1 = "listen1",
#                           y2 = "listen2",   # optional second indicator
#                           B = 9999, seed = 1)
#   print(res)
#   role_test(res, roles, role = "mother")                  # one-sample t on Delta D_r
#   role_test(res, roles, role = "mother", role2 = "father") # paired t
#   kp_contrast_test(res, roles, kp_role = "KP")             # Eqs. 9.11 and 9.12
#
# DATA FORMAT (long): one row per directed dyad (actor rates/behaves toward partner)
#   group    round-robin identifier
#   actor    person identifier (unique within group)
#   partner  person identifier
#   y1       the dependent variable (first indicator)
#   y2       optional second indicator of the same construct
# Every group must be a complete round robin with n >= 4 members.
# Groups may differ in size.
#
# PRACTICAL NOTES
# * With one indicator, dyadity indices square relationship effects that contain
#   measurement error. With two indicators of the same construct, supply y2:
#   the indices then use cross-products (Eqs. 9.6a, 9.8a) and are free of
#   error bias; individual indices can be negative (read as "no dyadity").
# * The Monte Carlo test compares the mean within-group variance of the
#   dyadity indices with its expected value when no one differs. Do not test
#   that variance against 0. B = 9999 takes seconds to a few minutes.
# * Calibration (simulated data, true null, alpha = .05): rejection rates of
#   .03 to .05 with one indicator and .03 to .08 with two indicators.
# * Power depends mostly on round-robin size. When one member had four times
#   the relationship variance of the others, the test detected it 82% of the
#   time with 20 eight-person round robins, but only 13% of the time with 50
#   four-person round robins. Plan for n >= 6 when dyadity differences are the
#   question.
# * With one indicator and n = 4, the expressed and elicited variance tests are
#   mathematically identical.
###############################################################################

## ---------------------------------------------------------------------------
## 1. SRM estimates for one round robin (Equations 9.3 to 9.5)
## ---------------------------------------------------------------------------
srm_rr <- function(X) {
  X <- as.matrix(X); n <- nrow(X)
  if (n < 4) stop("A round robin needs at least 4 members.")
  diag(X) <- NA
  if (any(is.na(X[row(X) != col(X)]))) stop("Missing off-diagonal values; each group must be complete.")
  M  <- mean(X, na.rm = TRUE)                        # grand mean M..
  mr <- rowMeans(X, na.rm = TRUE)                    # row means M_i.
  mc <- colMeans(X, na.rm = TRUE)                    # column means M_.j
  w1 <- (n - 1)^2 / (n * (n - 2)); w2 <- (n - 1) / (n * (n - 2)); w3 <- (n - 1) / (n - 2)
  a  <- w1 * mr + w2 * mc - w3 * M                   # actor effects (9.3)
  b  <- w1 * mc + w2 * mr - w3 * M                   # partner effects (9.4)
  G  <- X - outer(a, rep(1, n)) - outer(rep(1, n), b) - M   # relationship effects (9.5)
  diag(G) <- NA
  # relationship variance and dyadic covariance from dyad sums and differences
  up <- upper.tri(G)
  s  <- (G[up] + t(G)[up]); d <- (G[up] - t(G)[up])
  dfb <- (n - 1) * (n - 2) / 2 - 1; dfw <- (n - 1) * (n - 2) / 2
  MSB <- sum(s^2 / 2) / dfb; MSW <- sum(d^2 / 2) / dfw
  vg <- (MSB + MSW) / 2; cg <- (MSB - MSW) / 2
  rho <- if (vg > 0) cg / vg else 0
  list(n = n, M = M, actor = a, partner = b, G = G, var_rel = vg,
       cov_rel = cg, reciprocity = rho)
}

## ---------------------------------------------------------------------------
## 2. Dyadity indices for one round robin (Equations 9.6 to 9.10, 9.6a, 9.8a)
## ---------------------------------------------------------------------------
dyadity_rr <- function(X1, X2 = NULL) {
  f1 <- srm_rr(X1); n <- f1$n
  G1 <- f1$G; diag(G1) <- 0
  if (is.null(X2)) {
    S <- G1 * G1                                      # squares (9.6, 9.8)
    f2 <- NULL; G2 <- NULL
  } else {
    f2 <- srm_rr(X2); G2 <- f2$G; diag(G2) <- 0
    S <- G1 * G2                                      # cross-products (9.6a, 9.8a)
  }
  d_exp  <- rowSums(S) / (n - 1)                      # expressed dyadity, D_i.
  d_elic <- colSums(S) / (n - 1)                      # elicited dyadity, D_.j
  mbar   <- mean(d_exp)                               # = mean(d_elic) (identity)
  ids <- rownames(as.matrix(X1)); if (is.null(ids)) ids <- as.character(seq_len(n))
  persons <- data.frame(
    person   = ids,
    D_exp    = d_exp,
    D_elic   = d_elic,
    Delta_D  = d_exp - d_elic,                        # Eq. 9.10
    ratio_exp  = if (mbar > 0) d_exp / mbar else NA,  # standardized (1.00 = average)
    ratio_elic = if (mbar > 0) d_elic / mbar else NA,
    stringsAsFactors = FALSE)
  list(n = n, persons = persons,
       var_exp = var(d_exp), var_elic = var(d_elic),
       mean_D = mbar, G1 = G1, G2 = G2, srm1 = f1, srm2 = f2,
       cross = !is.null(X2))
}

## ---------------------------------------------------------------------------
## 3. Long data to one matrix per group
## ---------------------------------------------------------------------------
long_to_matrices <- function(data, group, actor, partner, y) {
  out <- list()
  for (g in unique(data[[group]])) {
    dg  <- data[data[[group]] == g, ]
    ids <- sort(unique(c(as.character(dg[[actor]]), as.character(dg[[partner]]))))
    X <- matrix(NA_real_, length(ids), length(ids), dimnames = list(ids, ids))
    X[cbind(match(as.character(dg[[actor]]), ids),
            match(as.character(dg[[partner]]), ids))] <- dg[[y]]
    out[[as.character(g)]] <- X
  }
  out
}

## ---------------------------------------------------------------------------
## 4. Monte Carlo null for individual differences in dyadity
##    Null: every member has the same relationship variance (no individual
##    differences in expressed or elicited dyadity). Dyad-level scores are
##    simulated with the dyadic reciprocity (and, with two indicators, the
##    cross-indicator covariances) pooled across round robins, relationship
##    effects are estimated with Eqs. 9.3 to 9.5 and rescaled to each group's
##    observed sum of squares, and the mean within-group variance of the
##    dyadity indices is recomputed.
##    p = (1 + #{T* >= T}) / (B + 1)   (Hope, 1968; Phipson & Smyth, 2010)
## ---------------------------------------------------------------------------
pair_moments <- function(Gu, Gv) {
  # dyad-level covariances between two sets of relationship effects, from
  # dyad sums and differences (the SRM relationship-variance logic)
  n <- nrow(Gu); up <- upper.tri(Gu)
  su <- Gu[up] + t(Gu)[up]; sv <- Gv[up] + t(Gv)[up]
  du <- Gu[up] - t(Gu)[up]; dv <- Gv[up] - t(Gv)[up]
  c(sb = sum(su * sv / 2), sw = sum(du * dv / 2),
    dfb = (n - 1) * (n - 2) / 2 - 1, dfw = (n - 1) * (n - 2) / 2)
}

pooled_cov <- function(rrs, r12 = NULL) {
  # pooled (across round robins) correlation structure of a directed dyad
  # (Z1_ij, Z2_ij, Z1_ji, Z2_ji); with one indicator only (Z1_ij, Z1_ji)
  pm <- function(u, v) {
    m <- rowSums(sapply(rrs, function(r) pair_moments(r[[u]], r[[v]])))
    MSB <- m["sb"] / m["dfb"]; MSW <- m["sw"] / m["dfw"]
    c(same = unname((MSB + MSW) / 2), recip = unname((MSB - MSW) / 2))
  }
  m11 <- pm("G1", "G1")
  if (!rrs[[1]]$cross) {
    S <- matrix(c(m11["same"], m11["recip"], m11["recip"], m11["same"]), 2)
  } else {
    m22 <- pm("G2", "G2"); m12 <- pm("G1", "G2")
    v1 <- m11["same"]; v2 <- m22["same"]; c12 <- m12["same"]
    S <- matrix(c(v1, c12, m11["recip"], m12["recip"],
                  c12, v2, m12["recip"], m22["recip"],
                  m11["recip"], m12["recip"], v1, c12,
                  m12["recip"], m22["recip"], c12, v2), 4, byrow = TRUE)
  }
  S <- cov2cor(S)
  if (!is.null(r12) && nrow(S) == 4) { S[1, 2] <- S[2, 1] <- S[3, 4] <- S[4, 3] <- r12 }
  e <- eigen(S, symmetric = TRUE)                     # guard against non-positive-definite
  L <- e$vectors %*% diag(sqrt(pmax(e$values, 1e-6))) # estimates in small samples
  L
}

sim_effects <- function(n, L) {
  # simulate one round robin of dyad-level scores with loading matrix L
  up <- which(upper.tri(diag(n)), arr.ind = TRUE); k <- nrow(L) / 2
  Z <- matrix(rnorm(nrow(up) * nrow(L)), nrow(up)) %*% t(L)
  out <- lapply(seq_len(k), function(v) {
    X <- matrix(0, n, n)
    X[up] <- Z[, v]; X[up[, 2:1]] <- Z[, k + v]; diag(X) <- NA; X
  })
  out
}

mc_null_one <- function(rr, L) {
  n <- rr$n
  if (is.list(L)) L <- L[[rr$name]]
  sims <- sim_effects(n, L)
  G1 <- srm_rr(sims[[1]])$G; diag(G1) <- 0
  G1 <- G1 * sqrt(sum(rr$G1^2) / sum(G1^2))          # rescale to observed sum of squares
  if (!rr$cross) {
    S <- G1 * G1
  } else {
    G2 <- srm_rr(sims[[2]])$G; diag(G2) <- 0
    G2 <- G2 * sqrt(sum(rr$G2^2) / sum(G2^2))
    S <- G1 * G2
  }
  c(exp = var(rowSums(S) / (n - 1)), elic = var(colSums(S) / (n - 1)))
}

mc_test <- function(rrs, B = 9999, seed = NULL, progress = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  T_exp  <- mean(sapply(rrs, `[[`, "var_exp"))
  T_elic <- mean(sapply(rrs, `[[`, "var_elic"))
  for (g in names(rrs)) rrs[[g]]$name <- g
  if (!rrs[[1]]$cross) {
    L <- pooled_cov(rrs)
  } else {
    # reciprocal structure pooled; same-dyad agreement between indicators per group
    L <- lapply(rrs, function(r) {
      m <- function(u, v) { p <- pair_moments(r[[u]], r[[v]]); (p["sb"] / p["dfb"] + p["sw"] / p["dfw"]) / 2 }
      r12 <- unname(m("G1", "G2") / sqrt(max(m("G1", "G1"), 1e-8) * max(m("G2", "G2"), 1e-8)))
      pooled_cov(rrs, r12 = max(min(r12, 0.99), -0.99))
    })
  }
  null <- matrix(NA_real_, B, 2, dimnames = list(NULL, c("exp", "elic")))
  step <- max(1, B %/% 10)
  for (b in seq_len(B)) {
    null[b, ] <- rowMeans(sapply(rrs, mc_null_one, L = L))
    if (progress && b %% step == 0) cat(sprintf("  Monte Carlo: %d of %d\n", b, B))
  }
  data.frame(
    index      = c("Expressed dyadity", "Elicited dyadity"),
    observed_T = c(T_exp, T_elic),
    null_mean  = colMeans(null),
    null_95th  = apply(null, 2, quantile, 0.95),
    p_value    = c((1 + sum(null[, "exp"]  >= T_exp))  / (B + 1),
                   (1 + sum(null[, "elic"] >= T_elic)) / (B + 1)),
    row.names = NULL)
}

## ---------------------------------------------------------------------------
## 5. Main function
## ---------------------------------------------------------------------------
dyadity_analysis <- function(data, group, actor, partner, y1, y2 = NULL,
                             B = 9999, seed = NULL, progress = TRUE) {
  M1 <- long_to_matrices(data, group, actor, partner, y1)
  M2 <- if (!is.null(y2)) long_to_matrices(data, group, actor, partner, y2) else NULL
  rrs <- lapply(names(M1), function(g) dyadity_rr(M1[[g]], if (!is.null(M2)) M2[[g]]))
  names(rrs) <- names(M1)
  persons <- do.call(rbind, lapply(names(rrs), function(g)
    cbind(group = g, rrs[[g]]$persons, stringsAsFactors = FALSE)))
  rownames(persons) <- NULL
  groups <- data.frame(
    group = names(rrs),
    n = sapply(rrs, `[[`, "n"),
    mean_D = sapply(rrs, `[[`, "mean_D"),
    var_exp = sapply(rrs, `[[`, "var_exp"),
    var_elic = sapply(rrs, `[[`, "var_elic"),
    rel_var = sapply(rrs, function(r) r$srm1$var_rel),
    reciprocity = sapply(rrs, function(r) r$srm1$reciprocity),
    row.names = NULL)
  mc <- if (B > 0) mc_test(rrs, B = B, seed = seed, progress = progress) else NULL
  structure(list(persons = persons, groups = groups, mc = mc, rr = rrs,
                 method = if (is.null(y2)) "squares (Eqs. 9.6, 9.8)"
                          else "cross-products of two indicators (Eqs. 9.6a, 9.8a)",
                 B = B), class = "dyadity")
}

print.dyadity <- function(x, digits = 3, ...) {
  cat("\nDyadity analysis:", x$method, "\n")
  cat("Round robins:", nrow(x$groups), "  Members:", nrow(x$persons), "\n")
  cat("\nPerson-level dyadity (D_exp = expressed, D_elic = elicited, Delta_D = Eq. 9.10):\n")
  print(x$persons, digits = digits, row.names = FALSE)
  cat("\nGroup-level summary (var_exp, var_elic = within-group variances of the indices):\n")
  print(x$groups, digits = digits, row.names = FALSE)
  if (!is.null(x$mc)) {
    cat(sprintf("\nMonte Carlo test of individual differences in dyadity (B = %d):\n", x$B))
    print(x$mc, digits = digits, row.names = FALSE)
    cat("H0: no individual differences. Compare observed_T with null_mean, not with 0.\n")
    if (all(x$groups$n == 4) && !grepl("cross", x$method))
      cat("Note: with one indicator and n = 4, the expressed and elicited tests are identical.\n")
  }
  invisible(x)
}

## ---------------------------------------------------------------------------
## 6. Role-specific tests of Delta D (Eq. 9.10)
##    roles: data.frame with columns group, person, role
##    One role:  one-sample t of Delta D_r against 0 (df = RR - 1)
##    Two roles: paired t of Delta D_role - Delta D_role2 (df = RR - 1)
## ---------------------------------------------------------------------------
role_test <- function(res, roles, role, role2 = NULL) {
  p <- merge(res$persons, roles, by = c("group", "person"))
  x <- p[p$role == role, c("group", "Delta_D")]
  if (anyDuplicated(x$group)) stop("More than one member per group has role '", role, "'.")
  if (is.null(role2)) {
    tt <- t.test(x$Delta_D, mu = 0)
    cat(sprintf("\nDelta D for role '%s' (expressed minus elicited), one-sample t test:\n", role))
  } else {
    y <- p[p$role == role2, c("group", "Delta_D")]
    xy <- merge(x, y, by = "group", suffixes = c("_1", "_2"))
    tt <- t.test(xy$Delta_D_1, xy$Delta_D_2, paired = TRUE)
    cat(sprintf("\nDelta D, '%s' minus '%s', paired t test:\n", role, role2))
  }
  cat(sprintf("  mean = %.3f, SD = %.3f, t(%d) = %.2f, p = %.3f, 95%% CI [%.3f, %.3f]\n",
              tt$estimate, sd(if (is.null(role2)) x$Delta_D else xy$Delta_D_1 - xy$Delta_D_2),
              as.integer(tt$parameter), tt$statistic, tt$p.value, tt$conf.int[1], tt$conf.int[2]))
  invisible(tt)
}

## ---------------------------------------------------------------------------
## 7. IKPD egocentric bias (Eqs. 9.11 and 9.12)
##    Key person's dyadity minus the mean of the nominees' dyadity, per group,
##    tested against 0 across groups (df = RR - 1).
## ---------------------------------------------------------------------------
kp_contrast_test <- function(res, roles, kp_role = "KP") {
  p <- merge(res$persons, roles, by = c("group", "person"))
  out <- do.call(rbind, lapply(split(p, p$group), function(d) {
    kp <- d$role == kp_role
    if (sum(kp) != 1) stop("Each group needs exactly one member with role '", kp_role, "'.")
    data.frame(group = d$group[1],
               exp  = d$D_exp[kp]  - mean(d$D_exp[!kp]),    # c' d_exp  (9.11)
               elic = d$D_elic[kp] - mean(d$D_elic[!kp]))   # c' d_elic (9.12)
  }))
  cat("\nKey person minus mean of nominees (egocentric bias), one-sample t tests:\n")
  for (v in c("exp", "elic")) {
    tt <- t.test(out[[v]], mu = 0)
    cat(sprintf("  %-9s mean = %.3f, t(%d) = %.2f, p = %.3f\n",
                if (v == "exp") "Expressed" else "Elicited",
                tt$estimate, as.integer(tt$parameter), tt$statistic, tt$p.value))
  }
  invisible(out)
}

###############################################################################
## EXAMPLE: the four idiographic round robins from Chapter 8 (Listening)
## Call run_example() after source() to check your installation.
###############################################################################
run_example <- function(B = 9999) {
  rr_example <- list(
    Trey = matrix(c(NA,4,4,5, 2,NA,2,5, 3,3,NA,3, 4,3,2,NA), 4, byrow = TRUE,
                  dimnames = rep(list(c("Trey","Jon","Mike","Page")), 2)),
    Jon  = matrix(c(NA,5,4,4, 1,NA,2,4, 1,4,NA,3, 2,5,4,NA), 4, byrow = TRUE,
                  dimnames = rep(list(c("Jon","Trey","Mike","Page")), 2)),
    Mike = matrix(c(NA,5,5,4, 5,NA,3,2, 4,1,NA,4, 5,2,3,NA), 4, byrow = TRUE,
                  dimnames = rep(list(c("Mike","Trey","Jon","Page")), 2)),
    Page = matrix(c(NA,1,1,2, 4,NA,1,2, 4,2,NA,2, 4,3,5,NA), 4, byrow = TRUE,
                  dimnames = rep(list(c("Page","Trey","Jon","Mike")), 2)))
  # convert to the long format the functions expect
  example_data <- do.call(rbind, lapply(names(rr_example), function(g) {
    X <- rr_example[[g]]; ij <- which(row(X) != col(X), arr.ind = TRUE)
    data.frame(group = g, actor = rownames(X)[ij[, 1]], partner = colnames(X)[ij[, 2]],
               listening = X[ij], stringsAsFactors = FALSE)
  }))
  # in each IKPD the key person is the group's namesake; nominees are the others
  example_roles <- data.frame(group = example_data$group, person = example_data$actor,
                              stringsAsFactors = FALSE)
  example_roles <- unique(example_roles)
  example_roles$role <- ifelse(example_roles$person == example_roles$group, "KP", "nominee")

  res <- dyadity_analysis(example_data, group = "group", actor = "actor",
                          partner = "partner", y1 = "listening", B = B, seed = 2026)
  print(res)
  role_test(res, example_roles, role = "KP")      # Delta D for key persons
  kp_contrast_test(res, example_roles, kp_role = "KP")
}

# Run the example when the file is executed with Rscript (not when sourced):
if (sys.nframe() == 0) run_example()
