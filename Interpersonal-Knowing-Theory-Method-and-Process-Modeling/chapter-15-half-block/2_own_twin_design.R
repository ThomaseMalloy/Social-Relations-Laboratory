# Table 15.1 design: raters (rows) rate their own digital twin (the diagonal) and the other
# twins in the block (columns). The same code works for the reverse design (twins rating
# their own human and other humans): put the twins in the actor column.
source("half_block_functions.R")
d <- read.csv("data/example_own_twin.csv")

# 1. Own twin vs. other twins, side by side
cmp <- do.call(rbind, lapply(c("warm", "kind"), function(v) hb_own_compare(d, v)))
print(cmp, digits = 3)
write.csv(cmp, "output/own_twin_comparison.csv", row.names = FALSE)

# 2. Rater-level values for second-stage modeling (e.g., predicting the own-twin effect)
write.csv(hb_own_twin(d, "warm"), "data/own_twin_effects_warm.csv", row.names = FALSE)

# 3. Multilevel version: own twin as a fixed effect, so it is not counted as relationship variance
if (requireNamespace("lme4", quietly = TRUE)) {
  d$actor <- paste(d$block, d$actor); d$partner <- paste(d$block, d$partner)
  m <- lme4::lmer(warm ~ own + (1 | actor) + (1 | partner), data = d)
  print(summary(m)$coefficients); print(lme4::VarCorr(m))
}
