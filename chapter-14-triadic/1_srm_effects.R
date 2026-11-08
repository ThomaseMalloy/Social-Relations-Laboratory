# ==============================================================================
# Chapter 14, Step 1: First-stage SRM of triadic idiographic round robins
#
# Reads the KP's triadic assessments (Table 14.2), runs the SRM in every round
# robin (one round robin = one KP in one group), and writes the grand means,
# effect estimates, and variance components as data sets for second-stage
# modeling (Bond & Malloy, 2023).
#
# Input: either
#   long  - one row per KP x group x perceiver x target, columns
#           kp.id, group, perceiver, target, <construct 1>, <construct 2>, ...
#   wide  - one row per KP x group, cell columns <construct>_<perceiver>_<target>
#           (e.g. satisfaction_KP_N1 = KP's own view of N1's satisfaction;
#            satisfaction_N1_KP = KP's view of how N1 sees the KP's satisfaction)
# The key person must be labelled "KP". Nominees can have any labels.
#
# Output (data/):
#   triadic_person_effects.csv        kp.id, group, focal, role, <c>.a (perceiver), <c>.b (target)
#   triadic_relationship_effects.csv  kp.id, group, perceiver, target, <c>.g
#   triadic_kp_level.csv              one row per KP x group: grand mean, KP and mean
#                                     nominee effects, variance components, covariances
#   triadic_kp_wide.csv               one row per KP, columns suffixed _friends, _coworkers
#
# With use_TripleR = TRUE the effects are estimated with TripleR and checked
# against the base R equations; the script stops if they differ. With FALSE
# (or if TripleR is not installed) the base R equations are used. They are the
# same Warner, Kenny, & Stoto (1979) estimators.
#
# Run from the chapter-14-triadic folder.
# ==============================================================================

if (!dir.exists("data")) stop("Set the working directory to the chapter-14-triadic folder first.")
source("srm_triadic_functions.R")

# ---- Settings: change these to match your data ------------------------------
input_format <- "long"                              # "long" or "wide"
input_file   <- "data/example_triadic_long.csv"     # or "data/example_triadic_wide.csv"
kp_col       <- "kp.id"
group_col    <- "group"
perceiver_col<- "perceiver"                         # long format only
target_col   <- "target"                            # long format only
constructs   <- c("listening", "satisfaction")      # one or more
use_TripleR  <- TRUE
# -----------------------------------------------------------------------------

raw <- read.csv(input_file, stringsAsFactors = FALSE)
if (input_format == "wide") {
  d <- wide_to_long(raw, kp_col, group_col, constructs)
} else {
  d <- data.frame(kp.id = raw[[kp_col]], group = raw[[group_col]],
                  perceiver = raw[[perceiver_col]], target = raw[[target_col]],
                  raw[constructs], stringsAsFactors = FALSE)
}
check_triadic_data(d, constructs)
cat("Read", nrow(d), "triadic assessments from", length(unique(d$kp.id)), "key people in",
    length(unique(d$group)), "groups:", paste(unique(d$group), collapse = ", "), "\n")

# ---- Base R first-stage SRM -------------------------------------------------
res <- srm_all(d, constructs)

# ---- TripleR (optional) -----------------------------------------------------
if (use_TripleR && requireNamespace("TripleR", quietly = TRUE)) {
  library(TripleR)
  d$rr.id <- paste(d$kp.id, d$group, sep = ":")
  d$p.id  <- paste(d$rr.id, d$perceiver, sep = ":")   # person ids unique across round robins
  d$t.id  <- paste(d$rr.id, d$target, sep = ":")
  for (v in constructs) {
    f  <- as.formula(paste(v, "~ p.id * t.id | rr.id"))
    rr <- RR(f, data = d, na.rm = FALSE)
    cat("\n==== TripleR:", v, "(pooled over round robins) ====\n"); print(rr)
    e  <- as.data.frame(rr$effects)
    e  <- data.frame(key = as.character(e[[1]]), a = e[[paste0(v, ".a")]], b = e[[paste0(v, ".p")]])
    mine <- data.frame(key = paste(res$person$kp.id, res$person$group, res$person$focal, sep = ":"),
                       a0 = res$person[[paste0(v, ".a")]], b0 = res$person[[paste0(v, ".b")]])
    chk <- merge(mine, e, by = "key")
    if (nrow(chk) != nrow(mine)) stop("TripleR and base R person ids do not match for ", v)
    dif <- max(abs(chk$a - chk$a0), abs(chk$b - chk$b0))
    if (dif > 1e-6) stop("TripleR and base R effects differ for ", v, " (max ", signif(dif, 3), ")")
    cat("TripleR and base R perceiver and target effects agree for", v, "(max difference", signif(dif, 2), ")\n")
  }
} else if (use_TripleR) {
  cat("TripleR is not installed; using the base R round robin equations (same estimators).\n",
      "To install it: install.packages(\"TripleR\")\n")
}

# ---- Write the second-stage data sets ---------------------------------------
write.csv(res$person,       "data/triadic_person_effects.csv",       row.names = FALSE)
write.csv(res$relationship, "data/triadic_relationship_effects.csv", row.names = FALSE)
write.csv(res$kp,           "data/triadic_kp_level.csv",             row.names = FALSE)
write.csv(kp_wide(res$kp),  "data/triadic_kp_wide.csv",              row.names = FALSE)

# Self-checks: effects sum to zero, and mean + a + b + g reproduces every score
for (v in constructs) {
  s <- aggregate(res$person[, paste0(v, c(".a", ".b"))],
                 res$person[c("kp.id", "group")], sum)
  stopifnot(max(abs(as.matrix(s[, 3:4]))) < 1e-8)
  m  <- merge(d[c("kp.id", "group", "perceiver", "target", v)], res$relationship,
              by = c("kp.id", "group", "perceiver", "target"))
  m  <- merge(m, setNames(res$person[c("kp.id", "group", "focal", paste0(v, ".a"))],
                          c("kp.id", "group", "perceiver", "a")))
  m  <- merge(m, setNames(res$person[c("kp.id", "group", "focal", paste0(v, ".b"))],
                          c("kp.id", "group", "target", "b")))
  m  <- merge(m, res$kp[c("kp.id", "group", paste0(v, ".mean"))])
  rec <- m[[paste0(v, ".mean")]] + m$a + m$b + m[[paste0(v, ".g")]]
  stopifnot(max(abs(rec - m[[v]])) < 1e-8)
}
cat("Checks passed: effects sum to zero in every round robin; grand mean + perceiver +",
    "target + relationship effects reproduce every score.\n")
cat("Wrote data/triadic_person_effects.csv (", nrow(res$person), "rows ),",
    "data/triadic_relationship_effects.csv (", nrow(res$relationship), "rows ),\n",
    "data/triadic_kp_level.csv (", nrow(res$kp), "rows ), data/triadic_kp_wide.csv (",
    length(unique(res$kp$kp.id)), "rows )\n")
