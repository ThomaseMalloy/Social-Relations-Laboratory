# 03_second_stage.R
# Second-stage analyses of the first-stage SRM parameters (family is the unit of analysis),
# in the form of Tables 2-4 of the family views study:
#   Table 2: mean variance components per key person and construct, t tests across families
#   Table 3: mean actor and partner effects of each member in each key person's view
#   Table 4: self-other and other-other agreement in actor and partner effects
# Reads output/srm_params_long.rds (written by 02_srm_analysis.R).
for (f in c("00_population_values.R", "param_functions.R")) source(f)
long <- readRDS("output/srm_params_long.rds")
G <- length(unique(long$family))

# ---- Table 2: variance components ------------------------------------------------------------
t2 <- list()
for (k in ROLES) for (c in CONSTRUCTS) for (p in c("varA", "varP", "varR", "varE", "covAP", "covR")) {
  x <- long$estimate[long$source == k & long$variable == paste0(c, p)]
  m <- mean(x); se <- sd(x) / sqrt(length(x)); t <- m / se
  pval <- if (startsWith(p, "var")) pt(t, length(x) - 1, lower.tail = FALSE) else 2 * pt(-abs(t), length(x) - 1)
  t2[[length(t2) + 1]] <- data.frame(source = k, construct = c, parameter = p, estimate = m, se = se,
                                     t = t, df = length(x) - 1, p = pval)
}
t2 <- do.call(rbind, t2)
# standardized: variance / total; reciprocities as correlations of the pooled components
std <- do.call(rbind, lapply(split(t2, list(t2$source, t2$construct)), function(d) {
  v <- setNames(d$estimate, d$parameter); tot <- sum(v[c("varA", "varP", "varR", "varE")])
  d$standardized <- ifelse(startsWith(d$parameter, "var"), d$estimate / tot,
    ifelse(d$parameter == "covAP", v[["covAP"]] / sqrt(v[["varA"]] * v[["varP"]]), v[["covR"]] / v[["varR"]]))
  d }))
t2 <- merge(t2, std[, c("source", "construct", "parameter", "standardized")])
t2 <- t2[order(match(t2$source, ROLES), match(t2$construct, CONSTRUCTS),
               match(t2$parameter, c("varA", "varP", "varR", "varE", "covAP", "covR"))), ]
# population value used in the simulation (covariances shown as correlations)
t2$population <- mapply(function(k, c, p) {
  pv <- pop_var[pop_var$source == k & pop_var$construct == c, ]
  switch(p, varA = pv$varA, varP = pv$varP, varR = pv$varR, varE = pv$varE,
         covAP = pv$genrec, covR = pv$dyrec) }, t2$source, t2$construct, t2$parameter)
write.csv(t2, "output/table2_variance_components.csv", row.names = FALSE)

# ---- Table 3: role means of actor and partner effects ----------------------------------------
t3 <- aggregate(estimate ~ source + variable + target,
                data = long[long$variable %in% c("POSact", "NEGact", "POSprt", "NEGprt"), ], FUN = mean)
t3w <- reshape(t3, idvar = c("source", "variable"), timevar = "target", direction = "wide")
names(t3w) <- sub("estimate.", "", names(t3w))
t3w <- t3w[order(match(t3w$variable, c("POSact", "NEGact", "POSprt", "NEGprt")), match(t3w$source, ROLES)),
           c("variable", "source", ROLES)]
# population (Table 3) for comparison
t3w$pop <- apply(t3w, 1, function(r) {
  c <- substr(r[["variable"]], 1, 3); tab <- if (grepl("act", r[["variable"]])) role_actor[[c]] else role_partner[[c]]
  paste(sprintf("%.3f", tab[r[["source"]], ]), collapse = " ") })
write.csv(t3w, "output/table3_role_means.csv", row.names = FALSE)

# ---- Table 4: agreement --------------------------------------------------------------------
ag <- agreement(long)
write.csv(ag$actor_partner, "output/table4_agreement.csv", row.names = FALSE)

options(width = 120)
cat("Second stage,", G, "families\n\nTable 2 (variance components)\n")
print(transform(t2, estimate = round(estimate, 3), se = round(se, 3), t = round(t, 2),
                p = signif(p, 2), standardized = round(standardized, 3)), row.names = FALSE)
cat("\nTable 3 (mean actor and partner effects; columns = member; pop = population values)\n")
print(transform(t3w, M = round(M, 3), F = round(F, 3), O = round(O, 3), Y = round(Y, 3)), row.names = FALSE)
cat("\nTable 4 (agreement)\n"); print(transform(ag$actor_partner, self_other = round(self_other, 2),
                                                other_other = round(other_other, 2)), row.names = FALSE)
cat("\nRelationship-effect agreement between informants:", round(ag$relationship, 2), "\n")
