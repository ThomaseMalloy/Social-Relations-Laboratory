# ==============================================================================
# Chapter 14, Step 4 (optional): Second-stage models with the exported effects
#
# Examples of second-stage modeling (Bond & Malloy, 2023) with the data sets
# from Step 1. Replace them with the models your theory calls for.
#   A. Multilevel model of target effects (person file): is the KP's target
#      effect higher than the nominees' (egocentric bias), in both groups?
#   B. Path model (KP wide file, lavaan): does the KP's egocentric bias among
#      friends predict it among coworkers, controlling for the KP's perceiver
#      effect (the multivariate model sketched in "Multivariate Triadic Modeling")?
# Needs lme4 and lavaan.
# ==============================================================================

if (!dir.exists("data")) stop("Set the working directory to the chapter-14-triadic folder first.")
person <- read.csv("data/triadic_person_effects.csv", stringsAsFactors = FALSE)
w      <- read.csv("data/triadic_kp_wide.csv",        stringsAsFactors = FALSE)
sink("output/second_stage_example.txt", split = TRUE)

# ---- A. Multilevel model of satisfaction target effects ---------------------
if (requireNamespace("lme4", quietly = TRUE)) {
  person$group <- factor(person$group, levels = c("friends", "coworkers"))
  # Target effects sum to zero in each round robin, so a KP random intercept
  # would be zero; the model compares the KP with the nominees within KP.
  mA <- lme4::lmer(satisfaction.b ~ 0 + group + group:role + (1 | kp.id:group), data = person,
                   control = lme4::lmerControl(check.conv.singular = "ignore"))
  cat("A. Satisfaction target effects by role (KP vs nominee) and group\n")
  print(summary(mA)$coefficients, digits = 3)
  cat("   The model SEs treat the 4 effects in a round robin as independent; they are not\n",
      "  (they sum to zero), so use the paired t tests with KPs as units for inference:\n")
  for (g in levels(person$group)) {
    d <- person[person$group == g, ]
    kp <- d$satisfaction.b[d$role == "KP"]
    nm <- tapply(d$satisfaction.b[d$role != "KP"], d$kp.id[d$role != "KP"], mean)
    tt <- t.test(kp, nm[d$kp.id[d$role == "KP"]], paired = TRUE, alternative = "greater")
    cat(sprintf("   %-9s beta_KP - mean beta_N = %.3f, t(%d) = %.2f, one-sided p %s\n",
                g, tt$estimate, as.integer(tt$parameter), tt$statistic,
                ifelse(tt$p.value < .001, "< .001", sprintf("= %.3f", tt$p.value))))
  }
}

# ---- B. Path model on the KP-level wide file --------------------------------
if (requireNamespace("lavaan", quietly = TRUE)) {
  w$ego_friends   <- w$satisfaction.b.kp_friends   - w$satisfaction.b.nom_friends
  w$ego_coworkers <- w$satisfaction.b.kp_coworkers - w$satisfaction.b.nom_coworkers
  model <- '
    ego_coworkers ~ ego_friends + satisfaction.a.kp_coworkers
    ego_friends   ~~ satisfaction.a.kp_coworkers
  '
  fit <- lavaan::sem(model, data = w)
  cat("\nB. Path model: KP egocentric bias (satisfaction) among coworkers\n")
  print(lavaan::parameterEstimates(fit, standardized = TRUE)[, c("lhs", "op", "rhs", "est", "se", "z", "pvalue", "std.all")],
        digits = 3)
}
sink()
cat("Wrote output/second_stage_example.txt\n")
