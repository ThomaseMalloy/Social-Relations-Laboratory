# 01_simulate.R
# Simulate multi-informant, multi-round robin (MIRR) data for G families of four
# (mother, father, older child, younger child). SIMULATED DATA, not results.
# Output: data/mirr_raw_long.(csv|sav|rds)  one row per family x source x actor x partner
#         data/mirr_raw_wide.(csv|sav|rds)  one row per family; SOURCE_ITEM_ACTORPARTNER
for (f in c("00_population_values.R", "sim_functions.R", "param_functions.R")) source(f)
G <- 50; set.seed(2026)
params <- make_params(W_OTHER, W_SELF, W_REL)
raw <- simulate_families(G, params, round_to_scale = TRUE)
dir.create("data", showWarnings = FALSE)
item_lab <- c(POS1 = "positivity item 1", POS2 = "positivity item 2",
              NEG1 = "negativity item 1", NEG2 = "negativity item 2")
write_all_formats(raw, "data/mirr_raw_long",
  labels = c(family = "family (round robin group)", source = "key person who made the rating",
             actor = "member whose behavior is rated", partner = "member the behavior is directed to",
             item_lab))
wide <- raw_to_wide(raw)
lab <- setNames(sapply(names(wide)[-1], function(n) {
  p <- strsplit(n, "_")[[1]]
  paste0(ROLE_LABELS[p[1]], "'s rating: ", item_lab[p[2]], ", ",
         ROLE_LABELS[substr(p[3], 1, 1)], " toward ", ROLE_LABELS[substr(p[3], 2, 2)])
}), names(wide)[-1])
write_all_formats(wide, "data/mirr_raw_wide", labels = c(family = "family", lab))
cat("Simulated", G, "families:", nrow(raw), "ratings rows;", ncol(wide) - 1, "wide rating variables\n")
print(head(raw, 12))
