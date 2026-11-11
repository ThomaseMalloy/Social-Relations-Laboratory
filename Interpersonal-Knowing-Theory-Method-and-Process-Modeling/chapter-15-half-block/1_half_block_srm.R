# Half-block SRM analysis: block estimates pooled across blocks (Blocko-style), effects
# for second-stage modeling, and the lme4 alternative.
source("half_block_functions.R")
d <- read.csv("data/example_half_block.csv")          # replace with your data
constructs <- list(warmth = c("warm", "kind"), competence = c("capable", "effective"))

res <- hb_all(d, constructs)
s <- res$summary
s[, -1] <- round(s[, -1], 3)
print(s)
write.csv(res$summary, "output/half_block_results.csv", row.names = FALSE)
write.csv(res$blocks, "data/half_block_block_estimates.csv", row.names = FALSE)
write.csv(hb_effects(d, unlist(constructs)), "data/half_block_effects.csv", row.names = FALSE)

if (requireNamespace("lme4", quietly = TRUE)) {
  for (cn in names(constructs)) {
    cat("\nlme4,", cn, "\n"); print(hb_lme4(d, constructs[[cn]])$components, digits = 3)
  }
}
sink("output/half_block_results.txt"); print(s); sink()
