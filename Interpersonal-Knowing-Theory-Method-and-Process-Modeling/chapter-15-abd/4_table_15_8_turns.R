# Chapter 15, Tables 15.8 and 15.9: partner and relationship construct variance for five LLM
# platforms (Haiku 5.5, Sonnet 5.5, DeepSeek, Kimi K3, Nemotron 3) after 0, 2, 4 and 6
# conversational turns. Each file is one ABD study: 50 blocks of 3 high- and 3 low-status
# agents, 900 directed ratings of warm, kind, capable and effective.
source("abd_functions.R")
constructs <- list(warmth = c("warm", "kind"), competence = c("capable", "effective"))
platforms <- c(haiku = "Haiku 5.5", sonnet = "Sonnet 5.5", deepseek = "DeepSeek",
               kimi = "Kimi K3", nemotron = "Nemotron 3")
stars <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", ""))))
fmt <- function(x) sub("^(-?)0\\.", "\\1.", sprintf("%.2f", x))

out <- NULL; blocks <- list()
for (pl in names(platforms)) for (t in c(0, 2, 4, 6)) {
  d <- read.csv(sprintf("data/turn_study/%s_turns%d.csv", pl, t))
  names(d) <- c("block", "actor", "actor_group", "partner", "partner_group",
                "warm", "kind", "capable", "effective")
  res <- abd_all(d, constructs, c("High", "Low"))
  blocks[[paste(pl, t)]] <- res$blocks
  for (dir in c("High_Low", "Low_High")) for (k in names(constructs))
    for (cmp in c("actor", "partner", "relationship")) {
      r <- res$summary[res$summary$term == paste(dir, k, cmp, sep = "."), ]
      out <- rbind(out, data.frame(platform = platforms[[pl]], turns = t, direction = dir,
        construct = k, component = cmp, variance = r$estimate, relative = r$relative,
        t = r$t, df = r$df, p_one_tailed = r$p))
    }
}
write.csv(out, "output/table_15_8_long.csv", row.names = FALSE)

# Table 15.8 layout: relative partner and relationship variance with one-tailed stars
out$cell <- with(out, paste0(fmt(relative), stars(p_one_tailed)))
w <- reshape(out[out$component != "actor", c("platform", "turns", "direction", "construct", "component", "cell")],
             idvar = c("platform", "turns", "direction", "construct"), timevar = "component", direction = "wide")
w <- reshape(w, idvar = c("platform", "turns", "direction"), timevar = "construct", direction = "wide")
w <- w[order(w$direction, match(w$platform, platforms), w$turns), ]
names(w) <- sub("cell\\.", "", names(w))
write.csv(w, "output/table_15_8.csv", row.names = FALSE)
for (dir in c("High_Low", "Low_High")) {
  cat("\nTable 15.8,", if (dir == "High_Low") "Panel A: High -> Low" else "Panel B: Low -> High", "\n")
  print(w[w$direction == dir, c("platform", "turns", "partner.warmth", "relationship.warmth",
                                "partner.competence", "relationship.competence")], row.names = FALSE)
}
cat("\nActor variance (omitted from Table 15.8), largest relative value:",
    fmt(max(out$relative[out$component == "actor"])), "\n")

# Table 15.9: medians across the five platforms
m <- aggregate(relative ~ construct + turns + direction + component,
               data = out[out$component != "actor", ], FUN = median)
m$relative <- fmt(m$relative)
m2 <- reshape(m, idvar = c("construct", "turns", "direction"), timevar = "component", direction = "wide")
m2 <- reshape(m2, idvar = c("construct", "turns"), timevar = "direction", direction = "wide")
m2 <- m2[order(-xtfrm(m2$construct), m2$turns), ]
write.csv(m2, "output/table_15_9.csv", row.names = FALSE)
cat("\nTable 15.9 (medians across platforms)\n"); print(m2, row.names = FALSE)

# Changes across turns are tested on the same 50 blocks with paired t tests, for example
# the drop in Haiku High -> Low warmth partner variance from 0 to 2 turns:
print(abd_compare(blocks[["haiku 0"]], blocks[["haiku 2"]], "High_Low.warmth.partner"))
