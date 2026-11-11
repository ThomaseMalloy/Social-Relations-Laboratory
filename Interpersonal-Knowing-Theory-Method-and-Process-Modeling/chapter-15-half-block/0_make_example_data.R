# Simulate half-block example data with known values (100 blocks of 4 actors x 4 partners).
# Two constructs (warmth: warm, kind; competence: capable, effective).
source("half_block_functions.R")
set.seed(15)
simulate_half_block <- function(K = 100, na = 4, np = 4,
                                s2 = list(actor = c(.10, .05), partner = c(.60, .40), relationship = c(.15, .10)),
                                r  = c(actor = .50, partner = -.40, relationship = .20),
                                error = .20, own_effect = 0) {
  mk <- function(n, comp) {          # bivariate (warmth, competence) construct-level effects
    R <- matrix(c(1, r[[comp]], r[[comp]], 1), 2)
    (matrix(rnorm(n * 2), n) %*% chol(R)) %*% diag(sqrt(s2[[comp]]))
  }
  out <- vector("list", K)
  for (k in 1:K) {
    A <- mk(na, "actor"); P <- mk(np, "partner")
    g <- expand.grid(a = 1:na, p = 1:np); G <- mk(nrow(g), "relationship")
    own <- as.integer(g$a == g$p)    # used only by the own-twin example
    y <- function(c) 3 + A[g$a, c] + P[g$p, c] + G[, c] + own_effect * own + rnorm(nrow(g), 0, sqrt(error))
    out[[k]] <- data.frame(block = k, actor = paste0("A", k, "_", g$a), partner = paste0("P", k, "_", g$p),
                           own = own, warm = y(1), kind = y(1), capable = y(2), effective = y(2))
  }
  do.call(rbind, out)
}
d <- simulate_half_block()
write.csv(d, "data/example_half_block.csv", row.names = FALSE)
d2 <- simulate_half_block(own_effect = .5)
d2$partner <- sub("^P", "T", d2$partner)  # twins
write.csv(d2, "data/example_own_twin.csv", row.names = FALSE)
cat("True construct values: actor .10/.05, partner .60/.40, relationship .15/.10, error .20;",
    "\nactor r = .50, partner r = -.40, relationship r = .20; own-twin effect = .50\n",
    file = "output/true_values.txt")
