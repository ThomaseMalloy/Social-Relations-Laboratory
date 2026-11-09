# ABD SRM analysis of the Chapter 15 social-status study (Claude; 50 blocks of 3 high-
# and 3 low-status agents; zero acquaintance). Replace the file name and constructs
# for your own data.
source("abd_functions.R")
d <- read.csv("data/claude_status_abd.csv")
constructs <- list(warmth = c("warm", "kind"), competence = c("capable", "effective"))
groups <- c("High", "Low")

res <- abd_all(d, constructs, groups)
S <- res$summary
S[, -1] <- round(S[, -1], 3)
print(S, row.names = FALSE)
cat("\nPercent stable construct variance\n"); print(res$stable, digits = 3)

write.csv(res$summary, "output/abd_results.csv", row.names = FALSE)
write.csv(res$blocks, "data/abd_block_estimates.csv", row.names = FALSE)
write.csv(abd_effects(d, unlist(constructs), groups), "data/abd_effects.csv", row.names = FALSE)
sink("output/abd_results.txt"); print(S, row.names = FALSE); print(res$stable, digits = 3); sink()

# For conditions measured on the same blocks (e.g., 0, 2, 4, 6 turns), run abd_all() on each
# condition and compare block estimates with abd_compare(), for example:
#   r0 <- abd_all(d0, constructs, groups); r2 <- abd_all(d2, constructs, groups)
#   abd_compare(r0$blocks, r2$blocks, "High_Low.warmth.partner")
