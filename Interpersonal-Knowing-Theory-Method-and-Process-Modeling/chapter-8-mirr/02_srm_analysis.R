# 02_srm_analysis.R
# First-stage SRM of every idiographic round robin: one TripleR analysis per family x key
# person (4 key persons x 50 families = 200 round robins), latent constructs with two
# indicators (POS1/POS2, NEG1/NEG2) plus the positivity x negativity bivariate analysis.
# Then: TripleR's pooled (multi-group) analysis for each key-person role, printed.
#
# Output (folder output/):
#   tripler_output_by_family.txt   printed TripleR output, every family x key person
#   tripler_output_pooled.txt      printed TripleR multi-group output, one per key-person role
#   srm_params_long.(csv|sav|rds)  one row per family x source x parameter (x target)
#   srm_params_wide.(csv|sav|rds)  one row per family; all parameters of M, F, O, Y in that order
#   variable_dictionary.csv        every wide variable name with its meaning
#
# ENGINE: "TripleR" when the package is installed (install.packages("TripleR")), otherwise
# the built-in estimator in srm_functions.R, which follows the TripleR formulas. Run
# 04_compare_engines.R to confirm the two agree on your machine.
for (f in c("00_population_values.R", "srm_functions.R", "sim_functions.R", "param_functions.R")) source(f)
ENGINE <- if (requireNamespace("TripleR", quietly = TRUE)) "TripleR" else "builtin"
cat("Engine:", ENGINE, "\n")
raw <- readRDS("data/mirr_raw_long.rds")
dir.create("output", showWarnings = FALSE)

# ---- TripleR route ---------------------------------------------------------------------------
# Extract one family x key person into the structure srm_to_long() expects.
tripler_one <- function(dk) {
  dk$actor <- factor(dk$actor, levels = ROLES); dk$partner <- factor(dk$partner, levels = ROLES)
  est <- list(); printed <- character(0)
  for (c in CONSTRUCTS) {
    fml <- as.formula(sprintf("%s1/%s2 ~ actor*partner", c, c))
    fit <- TripleR::RR(fml, data = dk, varname = c, na.rm = TRUE)
    printed <- c(printed, capture.output(print(fit)))
    vc <- fit$varComp
    getv <- function(type) vc$estimate[match(type, vc$type)]
    vA <- getv("actor variance"); vP <- getv("partner variance"); vR <- getv("relationship variance")
    vE <- getv("error variance"); cAP <- getv("actor-partner covariance"); cR <- getv("relationship covariance")
    tot <- vA + vP + vR + vE
    eff <- fit$effects
    ida <- grep("\\.a$", names(eff), value = TRUE)[1]; idp <- grep("\\.p$", names(eff), value = TRUE)[1]
    idcol <- intersect(c("id", "actor.id"), names(eff))[1]
    a <- eff[[ida]][match(ROLES, as.character(eff[[idcol]]))]
    b <- eff[[idp]][match(ROLES, as.character(eff[[idcol]]))]
    er <- fit$effectsRel
    rcol <- grep("(\\.r|_rel|\\.rel)$", names(er), value = TRUE)[1]
    g <- matrix(NA, 4, 4)
    g[cbind(match(as.character(er$actor.id), ROLES), match(as.character(er$partner.id), ROLES))] <- er[[rcol]]
    if (anyNA(a) || anyNA(b) || sum(!is.na(g)) != 12)
      stop("Could not read TripleR effects; check names(fit$effects) and names(fit$effectsRel)")
    X <- (rr_matrix(dk, paste0(c, 1)) + rr_matrix(dk, paste0(c, 2))) / 2
    est[[c]] <- list(mean = mean(X, na.rm = TRUE), a = a, b = b, g = g,
      comp = c(varA = vA, varP = vP, varR = vR, varE = vE, covAP = cAP, covR = cR,
               stdA = vA / tot, stdP = vP / tot, stdR = vR / tot, stdE = vE / tot,
               genrec = safe_r(cAP, vA, vP), dyrec = safe_r(cR, vR, vR)))
  }
  fitb <- TripleR::RR(POS1/POS2 + NEG1/NEG2 ~ actor*partner, data = dk, na.rm = TRUE)
  printed <- c(printed, capture.output(print(fitb)))
  bv <- fitb$bivariate
  pick <- function(pattern) bv$estimate[grep(pattern, bv$type, ignore.case = TRUE)[1]]
  cx <- est$POS$comp; cy <- est$NEG$comp
  covAA <- pick("^actor-actor"); covPP <- pick("^partner-partner")
  covAP <- pick("^actor-partner"); covPA <- pick("^partner-actor")
  covRi <- pick("intra"); covRe <- pick("inter")
  est$biv <- c(covAA = covAA, covPP = covPP, covAP = covAP, covPA = covPA,
               covRintra = covRi, covRinter = covRe,
               rAA = safe_r(covAA, cx[["varA"]], cy[["varA"]]), rPP = safe_r(covPP, cx[["varP"]], cy[["varP"]]),
               rAP = safe_r(covAP, cx[["varA"]], cy[["varP"]]), rPA = safe_r(covPA, cx[["varP"]], cy[["varA"]]),
               rRintra = safe_r(covRi, cx[["varR"]], cy[["varR"]]), rRinter = safe_r(covRe, cx[["varR"]], cy[["varR"]]))
  list(est = est, printed = printed)
}

builtin_one <- function(dk) {
  est <- list()
  for (c in CONSTRUCTS) est[[c]] <- srm_latent(rr_matrix(dk, paste0(c, 1)), rr_matrix(dk, paste0(c, 2)))
  est$biv <- srm_latent_bivariate(est$POS, est$NEG)
  printed <- c(sprintf("%-10s %10s %10s", "", "POS", "NEG"),
               sprintf("%-10s %10.3f %10.3f", names(est$POS$comp), est$POS$comp, est$NEG$comp))
  list(est = est, printed = printed)
}

# ---- first stage: every family x key person -----------------------------------------------
long <- list(); txt <- character(0)
for (f in sort(unique(raw$family))) for (k in ROLES) {
  dk <- raw[raw$family == f & raw$source == k, ]
  r <- if (ENGINE == "TripleR") tripler_one(dk) else builtin_one(dk)
  txt <- c(txt, "", strrep("=", 78),
           sprintf("Family %d, key person %s (%s)  [%s]", f, k, ROLE_LABELS[k], ENGINE),
           strrep("=", 78), r$printed)
  long[[length(long) + 1]] <- srm_to_long(f, k, r$est)
}
writeLines(txt, "output/tripler_output_by_family.txt")
long <- do.call(rbind, long)

# ---- pooled multi-group analysis per key-person role (TripleR prints this) ------------------
if (ENGINE == "TripleR") {
  ptxt <- character(0)
  for (k in ROLES) {
    dk <- raw[raw$source == k, ]
    dk$actor <- paste(dk$family, dk$actor, sep = ":"); dk$partner <- paste(dk$family, dk$partner, sep = ":")
    for (c in CONSTRUCTS) {
      fit <- TripleR::RR(as.formula(sprintf("%s1/%s2 ~ actor*partner | family", c, c)),
                         data = dk, varname = c, na.rm = TRUE)
      ptxt <- c(ptxt, "", strrep("=", 78), sprintf("Key person %s (%s), %s, %d families",
                k, ROLE_LABELS[k], c, length(unique(dk$family))), strrep("=", 78),
                capture.output(print(fit)))
    }
  }
  writeLines(ptxt, "output/tripler_output_pooled.txt")
}

# ---- write long and wide parameter files ------------------------------------------------------
lab_long <- c(family = "family", source = "key person (source of the ratings)",
              variable = "construct and SRM parameter", target = "member or directed dyad (blank = whole round robin)",
              estimate = "first-stage SRM estimate", label = "meaning", name = "wide-format variable name")
write_all_formats(long, "output/srm_params_long", labels = lab_long)
wide <- long_to_wide(long)
dict <- unique(long[, c("name", "source", "variable", "target", "label")])
dict$label <- paste0(ROLE_LABELS[dict$source], "'s view: ", dict$label)
write.csv(dict, "output/variable_dictionary.csv", row.names = FALSE)
write_all_formats(wide, "output/srm_params_wide",
                  labels = c(family = "family", setNames(dict$label, dict$name)))
cat("First-stage SRM done:", length(unique(long$family)), "families x 4 key persons;",
    ncol(wide) - 1, "parameters per family (wide)\n")
