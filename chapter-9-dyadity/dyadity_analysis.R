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
#                           partner = "partner",
#                           y = "listen1",                 # one indicator, or
#                           # y = c("listen1", "listen2"), # two or more indicators
#                           B = 9999, seed = 1)
#   print(res)
#   verify_matrix(res)                                      # matrix-algebra checks
#   role_test(res, roles, role = "mother")                  # one-sample t on Delta D_r
#   role_test(res, roles, role = "mother", role2 = "father") # paired t
#   kp_contrast_test(res, roles, kp_role = "KP")             # Eqs. 9.11 and 9.12
#
# DATA FORMAT (long): one row per directed dyad (actor rates/behaves toward partner)
#   group    round-robin identifier
#   actor    person identifier (unique within group)
#   partner  person identifier
#   y        the dependent variable; with K >= 2 indicators of the same
#            construct, one column per indicator
# Every group must be a complete round robin with n >= 4 members.
# Groups may differ in size.
#
# PRACTICAL NOTES
# * With one indicator, dyadity indices square relationship effects that contain
#   measurement error. With K >= 2 indicators of the same construct, supply
#   them all in y: the indices then use cross-products averaged over the
#   K(K - 1)/2 pairs of indicators (Eqs. 9.6a, 9.8a) and are free of error
#   bias; individual indices can be negative (read as "no dyadity").
# * The Monte Carlo test compares the mean within-group variance of the
#   dyadity indices with its expected value when no one differs. Do not test
#   that variance against 0. B = 9999 takes seconds to a few minutes.
# * Calibration (simulated data, true null, alpha = .05): rejection rates of
#   .03 to .05 with one indicator and .02 to .08 with two or three indicators.
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
##    Xs: list of K >= 1 matrices (indicators of the same construct) with the
##    same members. K = 1: squared effects (9.6, 9.8). K >= 2: products of the
##    indicators' relationship effects averaged over the K(K - 1)/2 pairs
##    (9.6a, 9.8a), which removes measurement error.
## ---------------------------------------------------------------------------
dyadity_rr <- function(Xs) {
  if (is.matrix(Xs) || is.data.frame(Xs)) Xs <- list(Xs)
  fits <- lapply(Xs, srm_rr); n <- fits[[1]]$n; K <- length(Xs)
  Gs <- lapply(fits, function(f) { G <- f$G; diag(G) <- 0; G })
  S <- product_matrix(Gs)
  d_exp  <- rowSums(S) / (n - 1)                      # expressed dyadity, D_i.
  d_elic <- colSums(S) / (n - 1)                      # elicited dyadity, D_.j
  mbar   <- mean(d_exp)                               # = mean(d_elic) (identity)
  ids <- rownames(as.matrix(Xs[[1]])); if (is.null(ids)) ids <- as.character(seq_len(n))
  persons <- data.frame(
    person   = ids,
    D_exp    = d_exp,
    D_elic   = d_elic,
    Delta_D  = d_exp - d_elic,                        # Eq. 9.10
    ratio_exp  = if (mbar > 0) d_exp / mbar else NA,  # standardized (1.00 = average)
    ratio_elic = if (mbar > 0) d_elic / mbar else NA,
    stringsAsFactors = FALSE)
  list(n = n, K = K, persons = persons, S = S,
       var_exp = var(d_exp), var_elic = var(d_elic),
       mean_D = mbar, G = Gs, srm = fits, X = Xs)
}

product_matrix <- function(Gs) {
  # S = G o G (one indicator) or the average of G_k o G_l over pairs k < l
  K <- length(Gs)
  if (K == 1) return(Gs[[1]] * Gs[[1]])
  pr <- combn(K, 2)
  Reduce(`+`, lapply(seq_len(ncol(pr)), function(p) Gs[[pr[1, p]]] * Gs[[pr[2, p]]])) / ncol(pr)
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
##    simulated with the dyadic reciprocity (and, with several indicators, the
##    cross-indicator covariances) pooled across round robins; the agreement
##    between indicators within a dyad is taken from each group. Relationship
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

dyad_cov <- function(rrs, k, l) {
  # same-dyad and reciprocal covariance of indicators k and l, pooled over groups
  m <- rowSums(sapply(rrs, function(r) pair_moments(r$G[[k]], r$G[[l]])))
  MSB <- m["sb"] / m["dfb"]; MSW <- m["sw"] / m["dfw"]
  c(same = unname((MSB + MSW) / 2), recip = unname((MSB - MSW) / 2))
}

group_same_cor <- function(r) {
  # within-group correlation between indicators for the same directed dyad
  K <- r$K; R <- diag(K)
  if (K == 1) return(R)
  v <- sapply(1:K, function(k) { p <- pair_moments(r$G[[k]], r$G[[k]])
                                 max((p["sb"] / p["dfb"] + p["sw"] / p["dfw"]) / 2, 1e-8) })
  for (k in 1:(K - 1)) for (l in (k + 1):K) {
    p <- pair_moments(r$G[[k]], r$G[[l]])
    R[k, l] <- R[l, k] <- max(min(((p["sb"] / p["dfb"] + p["sw"] / p["dfw"]) / 2) / sqrt(v[k] * v[l]), 0.99), -0.99)
  }
  R
}

null_loadings <- function(rrs) {
  # one loading matrix per group for the 2K-vector (Z_1ij..Z_Kij, Z_1ji..Z_Kji)
  K <- rrs[[1]]$K
  same <- recip <- matrix(0, K, K)
  for (k in 1:K) for (l in k:K) {
    m <- dyad_cov(rrs, k, l)
    same[k, l] <- same[l, k] <- m["same"]; recip[k, l] <- recip[l, k] <- m["recip"]
  }
  sd <- sqrt(pmax(diag(same), 1e-8))
  Rrec <- recip / outer(sd, sd)
  lapply(rrs, function(r) {
    Rsame <- if (K == 1) diag(1) else group_same_cor(r)
    Sig <- rbind(cbind(Rsame, Rrec), cbind(Rrec, Rsame))
    e <- eigen(Sig, symmetric = TRUE)                  # guard against non-positive-definite
    e$vectors %*% diag(sqrt(pmax(e$values, 1e-6)), nrow(Sig))   # small-sample estimates
  })
}

sim_effects <- function(n, L) {
  # simulate one round robin of dyad-level scores for K indicators
  up <- which(upper.tri(diag(n)), arr.ind = TRUE); K <- nrow(L) / 2
  Z <- matrix(rnorm(nrow(up) * nrow(L)), nrow(up)) %*% t(L)
  lapply(seq_len(K), function(v) {
    X <- matrix(0, n, n)
    X[up] <- Z[, v]; X[up[, 2:1]] <- Z[, K + v]; diag(X) <- NA; X
  })
}

mc_null_one <- function(rr, L) {
  n <- rr$n
  sims <- sim_effects(n, L)
  Gs <- lapply(seq_len(rr$K), function(k) {
    G <- srm_rr(sims[[k]])$G; diag(G) <- 0
    G * sqrt(sum(rr$G[[k]]^2) / sum(G^2))             # rescale to observed sum of squares
  })
  S <- product_matrix(Gs)
  c(exp = var(rowSums(S) / (n - 1)), elic = var(colSums(S) / (n - 1)))
}

mc_test <- function(rrs, B = 9999, seed = NULL, progress = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  T_exp  <- mean(sapply(rrs, `[[`, "var_exp"))
  T_elic <- mean(sapply(rrs, `[[`, "var_elic"))
  Ls <- null_loadings(rrs)
  null <- matrix(NA_real_, B, 2, dimnames = list(NULL, c("exp", "elic")))
  step <- max(1, B %/% 10)
  for (b in seq_len(B)) {
    null[b, ] <- rowMeans(mapply(mc_null_one, rrs, Ls))
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
##    y: name of the dependent variable, or a vector of names of K >= 2
##       indicators of the same construct (e.g., y = c("listen1", "listen2")).
##    y2: kept for convenience; y = "a", y2 = "b" is the same as y = c("a", "b").
## ---------------------------------------------------------------------------
dyadity_analysis <- function(data, group, actor, partner, y, y2 = NULL,
                             B = 9999, seed = NULL, progress = TRUE) {
  ys <- c(y, y2)
  Ms <- lapply(ys, function(v) long_to_matrices(data, group, actor, partner, v))
  gnames <- names(Ms[[1]])
  rrs <- lapply(gnames, function(g) dyadity_rr(lapply(Ms, `[[`, g)))
  names(rrs) <- gnames
  persons <- do.call(rbind, lapply(names(rrs), function(g)
    cbind(group = g, rrs[[g]]$persons, stringsAsFactors = FALSE)))
  rownames(persons) <- NULL
  groups <- data.frame(
    group = names(rrs),
    n = sapply(rrs, `[[`, "n"),
    mean_D = sapply(rrs, `[[`, "mean_D"),
    var_exp = sapply(rrs, `[[`, "var_exp"),
    var_elic = sapply(rrs, `[[`, "var_elic"),
    rel_var = sapply(rrs, function(r) r$srm[[1]]$var_rel),
    reciprocity = sapply(rrs, function(r) r$srm[[1]]$reciprocity),
    row.names = NULL)
  mc <- if (B > 0) mc_test(rrs, B = B, seed = seed, progress = progress) else NULL
  K <- length(ys)
  structure(list(persons = persons, groups = groups, mc = mc, rr = rrs, K = K,
                 method = if (K == 1) "one indicator, squared effects (Eqs. 9.6, 9.8)"
                          else sprintf("%d indicators, cross-products (Eqs. 9.6a, 9.8a)", K),
                 B = B), class = "dyadity")
}

print.dyadity <- function(x, digits = 3, ...) {
  cat("\nDyadity analysis:", x$method, "\n")
  cat("Round robins:", nrow(x$groups), "  Members:", nrow(x$persons), "\n")
  cat("\nPerson-level dyadity (D_exp = expressed, D_elic = elicited, Delta_D = Eq. 9.10):\n")
  print(x$persons, digits = digits, row.names = FALSE)
  cat("\nGroup-level summary (var_exp, var_elic = within-group variances of the indices;\n")
  cat("rel_var and reciprocity are for the first indicator):\n")
  print(x$groups, digits = digits, row.names = FALSE)
  if (!is.null(x$mc)) {
    cat(sprintf("\nMonte Carlo test of individual differences in dyadity (B = %d):\n", x$B))
    print(x$mc, digits = digits, row.names = FALSE)
    cat("H0: no individual differences. Compare observed_T with null_mean, not with 0.\n")
    if (all(x$groups$n == 4) && x$K == 1)
      cat("Note: with one indicator and n = 4, the expressed and elicited tests are identical.\n")
  }
  invisible(x)
}

## ---------------------------------------------------------------------------
## 5b. Matrix-algebra verification (matrix specification of Chapter 9)
##     Recomputes every quantity with the matrix formulas and checks the
##     identities, round robin by round robin.
## ---------------------------------------------------------------------------
verify_matrix <- function(res, tol = 1e-8, show_C = TRUE) {
  ok_all <- TRUE
  for (g in names(res$rr)) {
    r <- res$rr[[g]]; n <- r$n
    one <- rep(1, n); J <- matrix(1, n, n); I <- diag(n); A <- J - I
    Cm <- matrix(c(1, -1 / (n - 1), -1 / (n - 1), 1), 2)          # mixing matrix
    W  <- (n - 1) / (n * (n - 2)) * matrix(c(n - 1, 1, 1, n - 1), 2)
    checks <- c()
    for (k in seq_len(r$K)) {
      X <- r$X[[k]]; X0 <- X; diag(X0) <- 0
      m_r <- as.vector(X0 %*% one) / (n - 1); m_c <- as.vector(t(X0) %*% one) / (n - 1)
      M <- as.numeric(t(one) %*% X0 %*% one) / (n * (n - 1))
      ab <- (W %x% I) %*% c(m_r, m_c) - (n - 1) / (n - 2) * M * c(one, one)   # 9.3, 9.4 stacked
      a <- ab[1:n]; b <- ab[n + 1:n]
      G <- (X0 - a %*% t(one) - one %*% t(b) - M * J) * A                    # 9.5
      checks <- c(checks,
        setNames(max(abs(a - r$srm[[k]]$actor)), paste0("alpha (9.3), indicator ", k)),
        setNames(max(abs(b - r$srm[[k]]$partner)), paste0("beta (9.4), indicator ", k)),
        setNames(max(abs(G - r$G[[k]])), paste0("Gamma (9.5), indicator ", k)),
        setNames(max(abs(G %*% one)), paste0("Gamma 1 = 0, indicator ", k)),
        setNames(max(abs(t(G) %*% one)), paste0("Gamma' 1 = 0, indicator ", k)))
    }
    S <- r$S
    d_exp <- as.vector(S %*% one) / (n - 1); d_elic <- as.vector(t(S) %*% one) / (n - 1)
    dd <- as.vector((S - t(S)) %*% one) / (n - 1)
    checks <- c(checks,
      "C %*% C^-1 = I"            = max(abs(Cm %*% W - I[1:2, 1:2])),
      "d_exp (9.6/9.6a)"          = max(abs(d_exp - r$persons$D_exp)),
      "d_elic (9.8/9.8a)"         = max(abs(d_elic - r$persons$D_elic)),
      "D_Act = D_Prt (9.7, 9.9)"  = abs(mean(d_exp) - mean(d_elic)),
      "Delta d (9.10)"            = max(abs(dd - r$persons$Delta_D)),
      "1' Delta d = 0"            = abs(sum(dd)))
    if (r$K == 1 && n == 4)
      checks <- c(checks, "S2_Act = S2_Prt (n = 4)" = abs(var(d_exp) - var(d_elic)))
    ok <- all(checks < tol); ok_all <- ok_all && ok
    cat(sprintf("\nRound robin %s (n = %d, indicators = %d): %s\n", g, n, r$K,
                if (ok) "all checks pass" else "CHECK FAILED"))
    if (show_C && g == names(res$rr)[1]) {
      cat("  Mixing matrix C:\n"); print(round(Cm, 4))
      cat("  C^-1 (weights of Eqs. 9.3 and 9.4):\n"); print(round(W, 4))
      cat(sprintf("  Weight on M.. = %.4f\n", (n - 1) / (n - 2)))
    }
    print(data.frame(check = names(checks), max_abs_difference = signif(checks, 3),
                     row.names = NULL), row.names = FALSE)
  }
  invisible(ok_all)
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
                          partner = "partner", y = "listening", B = B, seed = 2026)
  print(res)
  role_test(res, example_roles, role = "KP")      # Delta D for key persons
  kp_contrast_test(res, example_roles, kp_role = "KP")
  verify_matrix(res)
  invisible(res)
}

# Run the example when the file is executed with Rscript (not when sourced):
if (sys.nframe() == 0) run_example()
