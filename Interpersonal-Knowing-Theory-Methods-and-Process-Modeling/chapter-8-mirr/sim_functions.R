# sim_functions.R
# Generating model for multi-informant, multi-round robin (MIRR) family data.
# Each family has four members (M, F, O, Y). Each member is a key person (source) who rates
# all 12 directed relationships among the four members on two indicators of positivity
# (POS1, POS2) and two of negativity (NEG1, NEG2), 1-7 scale.
#
# Rating by key person k of member i's behavior toward member j, construct c, indicator q:
#   X = mean_c + roleA_k[i] + roleP_k[j] + a_k[i] + b_k[j] + g_k[i,j] + e_q
# role means (fixed, Table 3) + random actor, partner, relationship effects + indicator error.
# Agreement across informants: each key person's random effects are a weighted mix of a
# family-level shared effect (the same for all four informants) and the key person's own view:
#   effect_k = sqrt(w) * shared + sqrt(1 - w) * own_k
# with w = w_self when the key person rates themself and w_other otherwise (actor and partner),
# and w_rel for relationship effects. The mix keeps each key person's variances, reciprocities,
# and positivity-negativity correlations at their population values.
suppressPackageStartupMessages(library(Matrix))

near_pd <- function(R) as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
is_pd <- function(R) min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) > 1e-8

# random-effect parameters per key person, after removing the variance that the fixed role
# means already contribute (TripleR treats role differences as actor/partner variance)
make_params <- function(w_other, w_self, w_rel) {
  P <- list()
  for (k in ROLES) {
    vA <- vP <- rAP <- numeric(0)
    for (c in CONSTRUCTS) {
      pv <- pop_var[pop_var$source == k & pop_var$construct == c, ]
      fa <- role_actor[[c]][k, ]; fp <- role_partner[[c]][k, ]
      vA[c] <- pv$varA - sum(fa^2) / 3
      vP[c] <- pv$varP - sum(fp^2) / 3
      if (vA[c] <= 0 || vP[c] <= 0) stop("role means exceed variance for ", k, " ", c)
      cov_tot <- pv$genrec * sqrt(pv$varA * pv$varP)
      rAP[c] <- max(min((cov_tot - sum(fa * fp) / 3) / sqrt(vA[c] * vP[c]), .95), -.95)
    }
    # order: aPOS, bPOS, aNEG, bNEG
    Rap <- diag(4)
    Rap[1, 2] <- Rap[2, 1] <- rAP[["POS"]]; Rap[3, 4] <- Rap[4, 3] <- rAP[["NEG"]]
    Rap[1, 3] <- Rap[3, 1] <- r_PN$actor;   Rap[2, 4] <- Rap[4, 2] <- r_PN$partner
    x <- r_PN$actor * sqrt(max(rAP[["POS"]], 0) * max(rAP[["NEG"]], 0))
    Rap[1, 4] <- Rap[4, 1] <- Rap[2, 3] <- Rap[3, 2] <- x
    if (!is_pd(Rap)) Rap <- near_pd(Rap)
    pvP <- pop_var[pop_var$source == k & pop_var$construct == "POS", ]
    pvN <- pop_var[pop_var$source == k & pop_var$construct == "NEG", ]
    # order: g_ij POS, g_ji POS, g_ij NEG, g_ji NEG
    Rr <- matrix(c(1, pvP$dyrec, r_PN$rel_same, r_PN$rel_cross,
                   pvP$dyrec, 1, r_PN$rel_cross, r_PN$rel_same,
                   r_PN$rel_same, r_PN$rel_cross, 1, pvN$dyrec,
                   r_PN$rel_cross, r_PN$rel_same, pvN$dyrec, 1), 4, byrow = TRUE)
    if (!is_pd(Rr)) Rr <- near_pd(Rr)
    P[[k]] <- list(sd_ap = sqrt(c(vA[["POS"]], vP[["POS"]], vA[["NEG"]], vP[["NEG"]])),
                   sd_r = sqrt(c(pvP$varR, pvP$varR, pvN$varR, pvN$varR)),
                   sd_e = sqrt(c(POS = pvP$varE, NEG = pvN$varE)),
                   Rap = Rap, Rr = Rr, rAP_random = rAP)
  }
  Rbar_ap <- Reduce(`+`, lapply(P, `[[`, "Rap")) / 4
  Rbar_r  <- Reduce(`+`, lapply(P, `[[`, "Rr")) / 4
  own <- function(R, w) {                 # covariance of the key person's own view
    U <- (R - w * Rbar) / (1 - w); if (!is_pd(U)) U <- near_pd(U); U
  }
  for (k in ROLES) {
    Rbar <- Rbar_ap
    P[[k]]$U_other <- own(P[[k]]$Rap, w_other); P[[k]]$U_self <- own(P[[k]]$Rap, w_self)
    Rbar <- Rbar_r
    P[[k]]$U_rel <- own(P[[k]]$Rr, w_rel)
  }
  list(keys = P, Rbar_ap = Rbar_ap, Rbar_r = Rbar_r,
       w = c(other = w_other, self = w_self, rel = w_rel))
}

rmvn <- function(n, S) matrix(rnorm(n * ncol(S)), n) %*% chol(S)

# simulate G families; returns long data: one row per family x source x actor x partner
simulate_families <- function(G, params, round_to_scale = TRUE, first_id = 1) {
  pairs <- t(combn(4, 2)); out <- vector("list", G * 4); r <- 0
  w <- params$w
  for (f in seq_len(G)) {
    S_ap <- rmvn(4, params$Rbar_ap)                 # shared actor/partner, members x 4
    S_r  <- rmvn(nrow(pairs), params$Rbar_r)        # shared relationship, dyads x 4
    for (k in ROLES) {
      Pk <- params$keys[[k]]; kk <- match(k, ROLES)
      V <- matrix(NA, 4, 4)
      for (i in 1:4) {
        ws <- if (i == kk) w[["self"]] else w[["other"]]
        U <- if (i == kk) Pk$U_self else Pk$U_other
        V[i, ] <- sqrt(ws) * S_ap[i, ] + sqrt(1 - ws) * rmvn(1, U)
      }
      V <- sweep(V, 2, Pk$sd_ap, `*`)              # columns: aPOS, bPOS, aNEG, bNEG
      Vr <- sqrt(w[["rel"]]) * S_r + sqrt(1 - w[["rel"]]) * rmvn(nrow(pairs), Pk$U_rel)
      Vr <- sweep(Vr, 2, Pk$sd_r, `*`)
      g <- list(POS = matrix(NA, 4, 4), NEG = matrix(NA, 4, 4))
      for (p in seq_len(nrow(pairs))) {
        i <- pairs[p, 1]; j <- pairs[p, 2]
        g$POS[i, j] <- Vr[p, 1]; g$POS[j, i] <- Vr[p, 2]
        g$NEG[i, j] <- Vr[p, 3]; g$NEG[j, i] <- Vr[p, 4]
      }
      d <- expand.grid(partner = 1:4, actor = 1:4)[, 2:1]; d <- d[d$actor != d$partner, ]
      rownames(d) <- NULL
      for (c in CONSTRUCTS) {
        a <- V[, if (c == "POS") 1 else 3]; b <- V[, if (c == "POS") 2 else 4]
        true <- construct_mean[[c]] + role_actor[[c]][k, d$actor] + role_partner[[c]][k, d$partner] +
          a[d$actor] + b[d$partner] + g[[c]][cbind(d$actor, d$partner)]
        # rounding to integers adds about 1/12 to the error variance (Sheppard), so the
        # generating error variance is reduced by 1/12 when ratings are rounded
        sd_e <- if (round_to_scale) sqrt(max(Pk$sd_e[[c]]^2 - 1 / 12, .01)) else Pk$sd_e[[c]]
        for (q in 1:2) {
          x <- true + rnorm(nrow(d), 0, sd_e)
          if (round_to_scale) x <- pmin(pmax(round(x), 1), 7)
          d[[paste0(c, q)]] <- x
        }
      }
      r <- r + 1
      out[[r]] <- data.frame(family = first_id + f - 1, source = k,
                             actor = ROLES[d$actor], partner = ROLES[d$partner],
                             d[, c("POS1", "POS2", "NEG1", "NEG2")])
    }
  }
  out <- do.call(rbind, out); rownames(out) <- NULL; out
}

# matrix for one family x source x indicator (rows actor, columns partner, order M F O Y)
rr_matrix <- function(dk, var) {
  X <- matrix(NA, 4, 4, dimnames = list(ROLES, ROLES))
  X[cbind(match(dk$actor, ROLES), match(dk$partner, ROLES))] <- dk[[var]]
  X
}
