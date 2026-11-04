# ==============================================================================
# Chapter 10, Step 1: First-stage SRM analysis with TripleR
#
# Malloy & Dumais (2026) digital twins, distinct-persona condition.
# Three platforms (Kimi, DeepSeek, Nemotron), each with 50 round robins of
# n = 5 twins and 1,000 directed ratings of "warm" and "kind" (1-5).
#
# Input  (data/): kimi_, deepseek_, nemotron_distinct_round_robins.csv
#                 human_agreeableness.csv
# Output (data/): ch10_twin_warmth_srm_effects.csv, one row per twin with
#                 warm and kind actor and partner effects and the human's
#                 agreeableness score (the input to 2_lavaan_sem.R).
#
# Run from the chapter-10-twin-warmth folder.
# ==============================================================================

if (!dir.exists("data")) stop("Set the working directory to the chapter-10-twin-warmth folder first, e.g. setwd(\"C:/.../chapter-10-twin-warmth\"), or open chapter-10-twin-warmth.Rproj.")

library(TripleR)

raw_files <- c(Kimi     = "data/kimi_distinct_round_robins.csv",
               DeepSeek = "data/deepseek_distinct_round_robins.csv",
               Nemotron = "data/nemotron_distinct_round_robins.csv")
traits <- read.csv("data/human_agreeableness.csv")

effects <- list()
for (pl in names(raw_files)) {
  rr_data <- read.csv(raw_files[[pl]])

  # Univariate SRM analysis for each item; print() shows the variance
  # components and their significance tests across the 50 round robins.
  warm <- RR(warm ~ actor.id * partner.id | roundrobin, data = rr_data)
  kind <- RR(kind ~ actor.id * partner.id | roundrobin, data = rr_data)
  cat("\n==========", pl, "==========\n")
  print(warm)
  print(kind)

  # Actor and partner effects, one row per twin
  w <- warm$effects[, 1:4]; names(w) <- c("id", "group.id", "warm.actor", "warm.partner")
  k <- kind$effects[, 1:4]; names(k) <- c("id", "group.id", "kind.actor", "kind.partner")
  e <- merge(w, k, by = c("id", "group.id"))

  # Add the human trait
  h <- traits[traits$platform == pl, c("id", "group.id", "agreeableness")]
  e <- merge(e, h, by = c("id", "group.id"))
  e$platform <- pl

  # Checks: 250 twins, no missing values, effects sum to zero in each round robin
  stopifnot(nrow(e) == 250, !anyNA(e))
  stopifnot(all(abs(tapply(e$warm.actor, e$group.id, sum)) < 1e-8),
            all(abs(tapply(e$kind.partner, e$group.id, sum)) < 1e-8))

  effects[[pl]] <- e
}

out <- do.call(rbind, effects)
out <- out[order(out$platform, out$group.id, out$id),
           c("platform", "group.id", "id", "agreeableness",
             "warm.actor", "kind.actor", "warm.partner", "kind.partner")]
write.csv(out, "data/ch10_twin_warmth_srm_effects.csv", row.names = FALSE)
cat("\nWrote data/ch10_twin_warmth_srm_effects.csv (", nrow(out), "rows )\n")
