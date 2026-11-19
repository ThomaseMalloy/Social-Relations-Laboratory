# 04_check_recovery.R
# Does the generating model reproduce the population values? Simulates a LARGE number of
# families (default 2,000), with and without rounding to the 1-7 scale, runs the built-in
# first-stage SRM, and compares second-stage means with the population values:
# variances (Table 2), role means (Table 3), agreement (Table 4), reciprocities, and
# positivity-negativity correlations. Writes output/recovery_check.csv.
for (f in c("00_population_values.R", "srm_functions.R", "sim_functions.R", "param_functions.R")) source(f)
G <- if (length(commandArgs(TRUE))) as.integer(commandArgs(TRUE)[1]) else 2000
params <- make_params(W_OTHER, W_SELF, W_REL)
res <- list()
for (rounded in c(FALSE, TRUE)) {
  set.seed(99)
  long <- srm_mirr_builtin(simulate_families(G, params, round_to_scale = rounded))
  m <- aggregate(estimate ~ source + variable + target, data = long, FUN = mean)
  pick <- function(k, v, t = "") m$estimate[m$source == k & m$variable == v & m$target == t]
  for (k in ROLES) for (c in CONSTRUCTS) {
    pv <- pop_var[pop_var$source == k & pop_var$construct == c, ]
    vA <- pick(k, paste0(c, "varA")); vP <- pick(k, paste0(c, "varP")); vR <- pick(k, paste0(c, "varR"))
    add <- function(what, pop, est) res[[length(res) + 1]] <<- data.frame(
      rounded = rounded, source = k, construct = c, quantity = what, population = pop, recovered = est)
    add("actor variance", pv$varA, vA); add("partner variance", pv$varP, vP)
    add("relationship variance", pv$varR, vR); add("error variance", pv$varE, pick(k, paste0(c, "varE")))
    add("generalized reciprocity r", pv$genrec, pick(k, paste0(c, "covAP")) / sqrt(vA * vP))
    add("dyadic reciprocity r", pv$dyrec, pick(k, paste0(c, "covR")) / vR)
    for (t in ROLES) {
      add(paste("role mean, actor", t), role_actor[[c]][k, t], pick(k, paste0(c, "act"), t))
      add(paste("role mean, partner", t), role_partner[[c]][k, t], pick(k, paste0(c, "prt"), t))
    }
  }
  ag <- agreement(long)
  res[[length(res) + 1]] <- data.frame(rounded = rounded, source = "all", construct = "both",
    quantity = c("other-other agreement", "self-other agreement", "relationship agreement"),
    population = unname(agree_target[c("other_other", "self_other", "relationship")]),
    recovered = c(mean(ag$actor_partner$other_other), mean(ag$actor_partner$self_other), mean(ag$relationship)))
}
res <- do.call(rbind, res)
res$difference <- res$recovered - res$population
dir.create("output", showWarnings = FALSE)
write.csv(res, "output/recovery_check.csv", row.names = FALSE)
options(width = 120)
for (r in c(FALSE, TRUE)) {
  cat("\n=== ", G, " families, ", if (r) "ROUNDED to 1-7 integers" else "continuous (not rounded)", " ===\n", sep = "")
  x <- res[res$rounded == r & !grepl("role mean", res$quantity), ]
  print(transform(x[, -1], population = round(population, 3), recovered = round(recovered, 3),
                  difference = round(difference, 3)), row.names = FALSE)
  y <- res[res$rounded == r & grepl("role mean", res$quantity), ]
  cat("Role means (Table 3): largest absolute difference", round(max(abs(y$difference)), 3),
      "; mean absolute difference", round(mean(abs(y$difference)), 3), "\n")
}
