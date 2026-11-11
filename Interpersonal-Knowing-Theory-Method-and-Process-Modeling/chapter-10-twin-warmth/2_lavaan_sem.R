# ==============================================================================
# Chapter 10, Step 2: Second-stage multiple-group SEM with lavaan
#
# Input: data/ch10_twin_warmth_srm_effects.csv (from 1_TripleR_srm_effects.R),
#        one row per twin: platform, group.id (round robin), id,
#        agreeableness (human), warm/kind actor and partner effects.
#
# Model: human agreeableness -> latent twin warmth actor effect (WA)
#        -> latent twin warmth partner effect (WP), plus a direct path
#        agreeableness -> WP (accuracy). Platform is the group.
# Standard errors are cluster-robust by round robin.
# ==============================================================================

if (!dir.exists("data")) stop("Set the working directory to the chapter-10-twin-warmth folder first, e.g. setwd(\"C:/.../chapter-10-twin-warmth\"), or open chapter-10-twin-warmth.Rproj.")

library(lavaan)

d <- read.csv("data/ch10_twin_warmth_srm_effects.csv")
d$platform <- factor(d$platform, levels = c("Kimi", "DeepSeek", "Nemotron"))
d$rr <- paste(d$platform, d$group.id)          # unique round robin id

fit_mg <- function(model, ...) {
  sem(model, data = d, group = "platform", group.label = levels(d$platform),
      cluster = "rr", fixed.x = FALSE, ...)
}

measurement <- '
  WA =~ warm.actor + kind.actor
  WP =~ warm.partner + kind.partner
  warm.actor ~~ warm.partner'

# Model 1: mediation only
m1 <- paste(measurement, '
  WA ~ agreeableness
  WP ~ WA')

# Model 2: adds the accuracy path; all parameters free across platforms
m2 <- paste(measurement, '
  WA ~ agreeableness
  WP ~ WA + agreeableness')

# Model 5 (final): accuracy path and residual covariance equal across platforms
m5 <- '
  WA =~ warm.actor + kind.actor
  WP =~ warm.partner + kind.partner
  warm.actor ~~ c(r, r, r) * warm.partner
  WA ~ agreeableness
  WP ~ WA + c(acc, acc, acc) * agreeableness'

fits <- list(
  M1_mediation_only    = fit_mg(m1),
  M2_accuracy_free     = fit_mg(m2),
  M3_equal_loadings    = fit_mg(m2, group.equal = "loadings"),
  M4_equal_paths       = fit_mg(m2, group.equal = "regressions"),
  M5_final             = fit_mg(m5))

# Fit of all models (Table 10.5)
fm <- c("chisq.scaled", "df", "pvalue.scaled", "rmsea.robust",
        "rmsea.ci.lower.robust", "rmsea.ci.upper.robust",
        "cfi.robust", "tli.robust", "srmr")
round(t(sapply(fits, fitMeasures, fit.measures = fm)), 3)

# Scaled chi-square difference tests against Model 2
lavTestLRT(fits$M1_mediation_only, fits$M2_accuracy_free)
lavTestLRT(fits$M2_accuracy_free, fits$M3_equal_loadings)
lavTestLRT(fits$M2_accuracy_free, fits$M4_equal_paths)
lavTestLRT(fits$M2_accuracy_free, fits$M5_final)

# Final model estimates (Table 10.6, Figure 10.5)
summary(fits$M5_final, standardized = TRUE, rsquare = TRUE)
