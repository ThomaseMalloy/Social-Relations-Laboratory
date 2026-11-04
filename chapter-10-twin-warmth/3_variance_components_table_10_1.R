# ==============================================================================
# Chapter 10, Table 10.1: SRM variance partitioning of warm and kind ratings
#
# Base R only (no packages). For each round robin the actor, partner, and
# relationship variances are estimated with the Warner, Kenny, and Stoto (1979)
# round robin formulas. Each variance is then averaged over the 50 round robins
# on a platform and tested against zero with a one-tailed t test,
# t(49) = mean / (SD / sqrt(50)).
#
# The script also recomputes the actor and partner effects from the raw
# ratings and checks them against data/ch10_twin_warmth_srm_effects.csv,
# so the effects file can be verified without TripleR.
#
# Run from the chapter-10-twin-warmth folder.
# ==============================================================================

if (!dir.exists("data")) stop("Set the working directory to the chapter-10-twin-warmth folder first, e.g. setwd(\"C:/.../chapter-10-twin-warmth\"), or open chapter-10-twin-warmth.Rproj.")

raw_files <- c(Kimi     = "data/kimi_distinct_round_robins.csv",
               DeepSeek = "data/deepseek_distinct_round_robins.csv",
               Nemotron = "data/nemotron_distinct_round_robins.csv")

# Round robin matrix for one group: rows = actors, columns = partners
rr_matrix <- function(g, y) {
  ids <- sort(unique(c(g$actor.id, g$partner.id))); n <- length(ids)
  M <- matrix(NA, n, n, dimnames = list(ids, ids))
  M[cbind(match(g$actor.id, ids), match(g$partner.id, ids))] <- g[[y]]
  M
}

# Actor and partner effects (Warner et al., 1979)
srm_effects <- function(M) {
  n <- nrow(M)
  ri <- rowMeans(M, na.rm = TRUE); cj <- colMeans(M, na.rm = TRUE); m <- mean(M, na.rm = TRUE)
  k1 <- (n - 1)^2 / (n * (n - 2)); k2 <- (n - 1) / (n * (n - 2)); k3 <- (n - 1) / (n - 2)
  list(a = k1 * ri + k2 * cj - k3 * m, b = k1 * cj + k2 * ri - k3 * m, m = m)
}

# Variance components for one round robin
srm_variances <- function(M) {
  n <- nrow(M); e <- srm_effects(M)
  G <- M - outer(e$a, rep(1, n)) - outer(rep(1, n), e$b) - e$m   # relationship effects
  up <- upper.tri(G)
  sums  <- (G + t(G))[up]; diffs <- (G - t(G))[up]
  MSB <- sum(sums^2 / 2)  / ((n - 1) * (n - 2) / 2 - 1)
  MSW <- sum(diffs^2 / 2) / ((n - 1) * (n - 2) / 2)
  c(actor        = sum(e$a^2) / (n - 1) - MSB / (2 * (n - 2)) - MSW / (2 * n),
    partner      = sum(e$b^2) / (n - 1) - MSB / (2 * (n - 2)) - MSW / (2 * n),
    relationship = (MSB + MSW) / 2)
}

effects_file <- read.csv("data/ch10_twin_warmth_srm_effects.csv")
table_rows <- list()

for (pl in names(raw_files)) {
  r <- read.csv(raw_files[[pl]])
  groups <- split(r, r$roundrobin)
  for (y in c("warm", "kind")) {
    V <- t(sapply(groups, function(g) srm_variances(rr_matrix(g, y))))
    mu <- colMeans(V); se <- apply(V, 2, sd) / sqrt(nrow(V)); tval <- mu / se
    pct <- 100 * mu / sum(mu)
    table_rows[[paste(pl, y)]] <- data.frame(
      platform = pl, item = y, M = mean(r[[y]]), SD = sd(r[[y]]),
      actor_pct = pct[["actor"]], partner_pct = pct[["partner"]], relationship_pct = pct[["relationship"]],
      actor_t = tval[["actor"]], partner_t = tval[["partner"]], relationship_t = tval[["relationship"]],
      actor_p = pt(tval[["actor"]], nrow(V) - 1, lower.tail = FALSE),
      partner_p = pt(tval[["partner"]], nrow(V) - 1, lower.tail = FALSE))

    # Check the effects file against effects recomputed from the raw ratings
    mine <- do.call(rbind, lapply(groups, function(g) {
      e <- srm_effects(rr_matrix(g, y))
      data.frame(id = as.integer(names(e$a)), a = e$a, b = e$b)
    }))
    ef <- effects_file[effects_file$platform == pl, ]
    chk <- merge(mine, ef, by = "id")
    stopifnot(nrow(chk) == 250,
              max(abs(chk$a - chk[[paste0(y, ".actor")]]))   < 1e-10,
              max(abs(chk$b - chk[[paste0(y, ".partner")]])) < 1e-10)
  }
}

table_10_1 <- do.call(rbind, table_rows); rownames(table_10_1) <- NULL
print(format(table_10_1, digits = 3), row.names = FALSE)
cat("\nEffects in data/ch10_twin_warmth_srm_effects.csv match the raw ratings.\n")
dir.create("output", showWarnings = FALSE)
write.csv(table_10_1, "output/table_10_1_variance_components.csv", row.names = FALSE)
