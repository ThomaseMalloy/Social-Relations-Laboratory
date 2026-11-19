# 05_compare_engines.R
# Run on a machine with TripleR installed, after 02_srm_analysis.R has run with TripleR.
# Recomputes every first-stage parameter with the built-in estimator and compares it with
# the TripleR values in output/srm_params_long.rds. Differences should be ~0 (rounding).
for (f in c("00_population_values.R", "srm_functions.R", "sim_functions.R", "param_functions.R")) source(f)
if (!requireNamespace("TripleR", quietly = TRUE)) stop("TripleR is not installed: install.packages('TripleR')")
raw <- readRDS("data/mirr_raw_long.rds")
tr <- readRDS("output/srm_params_long.rds")
bi <- srm_mirr_builtin(raw)
m <- merge(tr[, c("family", "name", "estimate")], bi[, c("family", "name", "estimate")],
           by = c("family", "name"), suffixes = c(".TripleR", ".builtin"))
m$diff <- m$estimate.TripleR - m$estimate.builtin
m$variable <- sub("^[MFOY]_", "", sub("_[MFOY]{1,2}$", "", m$name))
s <- aggregate(diff ~ variable, data = m, FUN = function(x) max(abs(x), na.rm = TRUE))
names(s)[2] <- "max_abs_difference"
print(s[order(-s$max_abs_difference), ], row.names = FALSE)
cat("\nCompared", nrow(m), "values; largest absolute difference:", max(abs(m$diff), na.rm = TRUE), "\n")
write.csv(m, "output/compare_TripleR_builtin.csv", row.names = FALSE)
