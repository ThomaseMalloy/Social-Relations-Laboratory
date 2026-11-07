# ==============================================================================
# Chapter 13, Step 1: First-stage SRM with TripleR, exported for ARRMA
#
# Input : a long round robin file, one row per directed pair, with
#           a group (round robin) id, an actor (perceiver) id, a partner
#           (target) id, the judgment of the partner, and the metaperception
#           (the actor's prediction of the partner's judgment of the actor).
# Output (data/):
#   arrma_person_effects.csv        group.id, id, judgment.a, judgment.p, meta.a, meta.p
#   arrma_relationship_effects.csv  group.id, actor.id, partner.id, judgment.r, meta.r
#
# These two files are the input to 2_arrma_estimate.R. If you already have
# TripleR output saved elsewhere, you can skip this step and point Step 2 at
# your own files (see the settings block at the top of Step 2).
#
# Run from the chapter-13-arrma folder. Needs the TripleR package.
# ==============================================================================

if (!dir.exists("data")) stop("Set the working directory to the chapter-13-arrma folder first.")
library(TripleR)

# ---- Settings: change these to match your file ------------------------------
raw_file   <- "data/example_round_robins.csv"
group_col  <- "group.id"
actor_col  <- "actor.id"
partner_col<- "partner.id"
judg_col   <- "judgment"   # i's judgment of j
meta_col   <- "meta"       # i's prediction of j's judgment of i
# -----------------------------------------------------------------------------

d <- read.csv(raw_file)
d <- data.frame(group.id = d[[group_col]], actor.id = d[[actor_col]],
                partner.id = d[[partner_col]],
                judgment = d[[judg_col]], meta = d[[meta_col]])

# Univariate SRM analyses; print() shows variance components and tests.
rr_j <- RR(judgment ~ actor.id * partner.id | group.id, data = d)
rr_m <- RR(meta     ~ actor.id * partner.id | group.id, data = d)
print(rr_j); print(rr_m)

# Person level: TripleR's $effects has id, group.id, <var>.a, <var>.p first
ej <- rr_j$effects[, 1:4]; names(ej) <- c("id", "group.id", "judgment.a", "judgment.p")
em <- rr_m$effects[, 1:4]; names(em) <- c("id", "group.id", "meta.a", "meta.p")
person <- merge(ej, em, by = c("id", "group.id"))

# Dyad level: TripleR's $effectsRel holds the relationship effects
pick_rel <- function(x, newname) {
  x <- as.data.frame(x)
  a  <- grep("actor|perceiver", names(x), value = TRUE)[1]
  p  <- grep("partner|target",  names(x), value = TRUE)[1]
  g  <- grep("group",           names(x), value = TRUE)[1]
  r  <- setdiff(names(x)[sapply(x, is.numeric)], c(a, p, g, "dyad"))
  r  <- r[length(r)]                     # the relationship effect column
  out <- data.frame(x[[g]], x[[a]], x[[p]], x[[r]])
  names(out) <- c("group.id", "actor.id", "partner.id", newname)
  out
}
rel <- merge(pick_rel(rr_j$effectsRel, "judgment.r"),
             pick_rel(rr_m$effectsRel, "meta.r"),
             by = c("group.id", "actor.id", "partner.id"))

stopifnot(!anyNA(person), !anyNA(rel))
write.csv(person[, c("group.id", "id", "judgment.a", "judgment.p", "meta.a", "meta.p")],
          "data/arrma_person_effects.csv", row.names = FALSE)
write.csv(rel, "data/arrma_relationship_effects.csv", row.names = FALSE)
cat("Wrote data/arrma_person_effects.csv (", nrow(person), "people ) and\n",
    "data/arrma_relationship_effects.csv (", nrow(rel), "directed dyads )\n")
