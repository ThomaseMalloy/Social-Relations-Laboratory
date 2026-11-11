# ==============================================================================
# Chapter 13, Step 0 (optional): Example round robin data for ARRMA
#
# Simulates 40 round robins of n = 6 people. Each person rates every partner
# (judgment) and predicts how that partner rated them (metaperception). The
# data are generated from known individual and dyadic ARRMA parameters so the
# estimates in Step 2 can be checked against the truth.
#
# Output (data/):
#   example_round_robins.csv        raw long data, one row per directed pair
#                                    (group.id, actor.id, partner.id, judgment, meta)
#   arrma_person_effects.csv        SRM actor (.a) and partner (.p) effects
#   arrma_relationship_effects.csv  SRM relationship (.r) effects
#
# The two effects files are computed here in base R with the round robin
# equations of Chapter 4 (Equations 4.1-4.3; Warner, Kenny, & Stoto, 1979),
# which are the estimates TripleR reports. They have the same layout that
# 1_TripleR_export_effects.R writes, so Step 2 runs without TripleR.
#
# Run from the chapter-13-arrma folder. Base R only.
# ==============================================================================

if (!dir.exists("data")) stop("Set the working directory to the chapter-13-arrma folder first.")
set.seed(13)

G <- 40   # round robins
n <- 6    # people per round robin

# ---- True parameters --------------------------------------------------------
# Individual level (standardized SRM effects)
rho_A  <- .30   # generalized reciprocity: cor(actor judgment, partner judgment)
b_ind  <- .50   # unique generalized assumed reciprocity (b*)
c_ind  <- .25   # unique generalized meta-accuracy (c*)
# Dyadic level (standardized relationship effects)
rho_D  <- .40   # dyadic reciprocity: cor(gamma_ij, gamma_ji)
b_dy   <- .45   # unique dyadic assumed reciprocity (b*)
c_dy   <- .20   # unique dyadic meta-accuracy (c*)

# Residual SDs that make the metaperception effects standardized (variance 1)
res_sd <- function(b, c, r) sqrt(1 - (b^2 + c^2 + 2 * b * c * r))

# Variance shares of the raw scores (actor, partner, relationship, error)
sd_a <- sqrt(.30); sd_p <- sqrt(.15); sd_r <- sqrt(.45); sd_e <- sqrt(.10)

rows <- list()
for (g in seq_len(G)) {
  # Person effects for judgment (actor, partner) with reciprocity rho_A
  z1 <- rnorm(n); z2 <- rho_A * z1 + sqrt(1 - rho_A^2) * rnorm(n)
  a_j <- z1; p_j <- z2
  # Metaperception actor effect from ARRMA: a_mp = b* a_j + c* p_j + e
  a_m <- b_ind * a_j + c_ind * p_j + res_sd(b_ind, c_ind, rho_A) * rnorm(n)
  p_m <- rnorm(n)   # partner effect in metaperception (not part of ARRMA)

  # Relationship effects: g_ij and g_ji correlated rho_D
  gj <- matrix(NA, n, n); gm <- matrix(NA, n, n)
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    u <- rnorm(1); v <- rho_D * u + sqrt(1 - rho_D^2) * rnorm(1)
    gj[i, j] <- u; gj[j, i] <- v
  }
  # i's metaperception of j's response to i: gm[i, j] = b* g_ij + c* g_ji + e
  for (i in 1:n) for (j in 1:n) if (i != j)
    gm[i, j] <- b_dy * gj[i, j] + c_dy * gj[j, i] + res_sd(b_dy, c_dy, rho_D) * rnorm(1)

  for (i in 1:n) for (j in 1:n) if (i != j) {
    rows[[length(rows) + 1]] <- data.frame(
      group.id = g, actor.id = paste0(g, "_", i), partner.id = paste0(g, "_", j),
      judgment = 5 + sd_a * a_j[i] + sd_p * p_j[j] + sd_r * gj[i, j] + sd_e * rnorm(1),
      meta     = 5 + sd_a * a_m[i] + sd_p * p_m[j] + sd_r * gm[i, j] + sd_e * rnorm(1))
  }
}
raw <- do.call(rbind, rows)
raw$judgment <- round(raw$judgment, 3); raw$meta <- round(raw$meta, 3)
write.csv(raw, "data/example_round_robins.csv", row.names = FALSE)

# ---- First-stage SRM in base R (Chapter 4, Equations 4.1-4.3) ---------------
srm_round_robin <- function(d, var) {
  ids <- sort(unique(c(d$actor.id, d$partner.id))); n <- length(ids)
  X <- matrix(NA, n, n, dimnames = list(ids, ids))
  X[cbind(match(d$actor.id, ids), match(d$partner.id, ids))] <- d[[var]]
  Mi <- rowMeans(X, na.rm = TRUE); Mj <- colMeans(X, na.rm = TRUE); M <- mean(X, na.rm = TRUE)
  k1 <- (n - 1)^2 / (n * (n - 2)); k2 <- (n - 1) / (n * (n - 2)); k3 <- (n - 1) / (n - 2)
  a <- k1 * Mi + k2 * Mj - k3 * M          # actor effects   (4.1)
  p <- k1 * Mj + k2 * Mi - k3 * M          # partner effects (4.2)
  R <- X - outer(a, rep(1, n)) - outer(rep(1, n), p) - M   # relationship (4.3)
  list(person = data.frame(id = ids, a = a, p = p),
       rel = data.frame(actor.id = rep(ids, n), partner.id = rep(ids, each = n),
                        r = as.vector(R))[!is.na(as.vector(R)), ])
}

pe <- list(); re <- list()
for (g in unique(raw$group.id)) {
  d  <- raw[raw$group.id == g, ]
  sj <- srm_round_robin(d, "judgment"); sm <- srm_round_robin(d, "meta")
  pe[[length(pe) + 1]] <- data.frame(group.id = g, id = sj$person$id,
    judgment.a = sj$person$a, judgment.p = sj$person$p,
    meta.a = sm$person$a, meta.p = sm$person$p)
  r <- merge(setNames(sj$rel, c("actor.id", "partner.id", "judgment.r")),
             setNames(sm$rel, c("actor.id", "partner.id", "meta.r")))
  re[[length(re) + 1]] <- data.frame(group.id = g, r)
}
write.csv(do.call(rbind, pe), "data/arrma_person_effects.csv", row.names = FALSE)
write.csv(do.call(rbind, re), "data/arrma_relationship_effects.csv", row.names = FALSE)

cat("Wrote data/example_round_robins.csv (", nrow(raw), "directed ratings ),\n",
    "data/arrma_person_effects.csv and data/arrma_relationship_effects.csv\n")
cat("True values: individual b* =", b_ind, " c* =", c_ind, " reciprocity =", rho_A, "\n",
    "            dyadic     b* =", b_dy,  " c* =", c_dy,  " reciprocity =", rho_D, "\n")
