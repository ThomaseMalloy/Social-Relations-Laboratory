# param_functions.R
# Collect first-stage SRM parameters into the standard naming format
#   SOURCE_VARIABLE_TARGET
# SOURCE = key person who made the ratings (M, F, O, Y)
# VARIABLE = construct + parameter, e.g. POSact, POSprt, POSrel, POSvarA, POSgenrec
# TARGET = member the parameter is about (M, F, O, Y) or the directed dyad (FO = father toward
#          older child); omitted for parameters of the whole round robin (means, variances).
# Order everywhere: sources M, F, O, Y; constructs POS, NEG; members M, F, O, Y;
# dyads MF, MO, MY, FM, FO, FY, OM, OF, OY, YM, YF, YO.

DYADS <- c("MF", "MO", "MY", "FM", "FO", "FY", "OM", "OF", "OY", "YM", "YF", "YO")

# one row-set per family x source, from a list with elements mean, a, b, g, comp, biv
# (both the built-in estimator and the TripleR extractor return this structure)
srm_to_long <- function(family, source, est) {
  rows <- list()
  add <- function(variable, target, value, label)
    rows[[length(rows) + 1]] <<- data.frame(family = family, source = source,
      variable = variable, target = target, estimate = unname(value), label = label)
  for (c in CONSTRUCTS) {
    e <- est[[c]]; cl <- if (c == "POS") "positivity" else "negativity"
    add(paste0(c, "mean"), "", e$mean, paste("within-person grand mean,", cl))
    comp_lab <- c(varA = "actor variance", varP = "partner variance",
                  varR = "relationship variance", varE = "error variance",
                  covAP = "actor-partner covariance", covR = "relationship covariance",
                  stdA = "actor variance, proportion of total",
                  stdP = "partner variance, proportion of total",
                  stdR = "relationship variance, proportion of total",
                  stdE = "error variance, proportion of total",
                  genrec = "generalized reciprocity (actor-partner r)",
                  dyrec = "dyadic reciprocity (relationship r)")
    for (p in names(comp_lab)) add(paste0(c, p), "", e$comp[[p]], paste0(comp_lab[[p]], ", ", cl))
    for (i in seq_along(ROLES)) add(paste0(c, "act"), ROLES[i], e$a[i],
                                    paste0("actor effect of ", ROLE_LABELS[ROLES[i]], ", ", cl))
    for (i in seq_along(ROLES)) add(paste0(c, "prt"), ROLES[i], e$b[i],
                                    paste0("partner effect of ", ROLE_LABELS[ROLES[i]], ", ", cl))
    for (dy in DYADS) {
      i <- match(substr(dy, 1, 1), ROLES); j <- match(substr(dy, 2, 2), ROLES)
      add(paste0(c, "rel"), dy, e$g[i, j], paste0("relationship effect, ", ROLE_LABELS[ROLES[i]],
                                                 " toward ", ROLE_LABELS[ROLES[j]], ", ", cl))
    }
  }
  biv_lab <- c(covAA = "actor-actor covariance", covPP = "partner-partner covariance",
               covAP = "POS actor - NEG partner covariance", covPA = "POS partner - NEG actor covariance",
               covRintra = "relationship covariance, same dyad direction",
               covRinter = "relationship covariance, reversed dyad direction",
               rAA = "actor-actor r", rPP = "partner-partner r",
               rAP = "POS actor - NEG partner r", rPA = "POS partner - NEG actor r",
               rRintra = "relationship r, same direction", rRinter = "relationship r, reversed direction")
  for (p in names(biv_lab)) add(paste0("PN", p), "", est$biv[[p]], paste0(biv_lab[[p]], ", positivity with negativity"))
  out <- do.call(rbind, rows)
  out$name <- ifelse(out$target == "", paste(out$source, out$variable, sep = "_"),
                     paste(out$source, out$variable, out$target, sep = "_"))
  out
}

# long -> wide: one row per family, columns in the fixed order of the long file
long_to_wide <- function(long) {
  ord <- unique(long$name)                       # order of first family = standard order
  w <- reshape(long[, c("family", "name", "estimate")], idvar = "family", timevar = "name",
               direction = "wide")
  names(w) <- sub("^estimate\\.", "", names(w))
  w <- w[order(w$family), c("family", ord)]
  rownames(w) <- NULL
  w
}

# raw ratings long -> wide (one row per family; SOURCE_ITEM_ACTORPARTNER)
raw_to_wide <- function(raw) {
  l <- do.call(rbind, lapply(c("POS1", "POS2", "NEG1", "NEG2"), function(v)
    data.frame(family = raw$family, name = paste0(raw$source, "_", v, "_", raw$actor, raw$partner),
               value = raw[[v]])))
  ord <- as.vector(outer(c("POS1", "POS2", "NEG1", "NEG2"), DYADS, function(v, d) paste0(v, "_", d)))
  ord <- as.vector(t(outer(ROLES, ord, paste, sep = "_")))
  w <- reshape(l, idvar = "family", timevar = "name", direction = "wide")
  names(w) <- sub("^value\\.", "", names(w))
  w <- w[order(w$family), c("family", ord)]; rownames(w) <- NULL
  w
}

# write CSV + SPSS (.sav, with variable labels) + R (.rds)
write_all_formats <- function(df, stem, labels = NULL) {
  write.csv(df, paste0(stem, ".csv"), row.names = FALSE, na = "")
  saveRDS(df, paste0(stem, ".rds"))
  if (requireNamespace("haven", quietly = TRUE)) {
    d2 <- df
    if (!is.null(labels)) for (v in intersect(names(labels), names(d2))) attr(d2[[v]], "label") <- labels[[v]]
    haven::write_sav(d2, paste0(stem, ".sav"))
  } else warning("haven not installed: SPSS file not written")
}

# ---- agreement across informants (Table 4 logic) -------------------------------------------
agreement <- function(long) {
  res <- list()
  for (c in CONSTRUCTS) for (eff in c("act", "prt")) for (t in ROLES) {
    v <- long[long$variable == paste0(c, eff) & long$target == t, ]
    W <- reshape(v[, c("family", "source", "estimate")], idvar = "family", timevar = "source",
                 direction = "wide"); names(W) <- sub("estimate.", "", names(W))
    others <- setdiff(ROLES, t)
    so <- cor(W[[t]], rowMeans(W[, others]))
    cc <- combn(others, 2); oo <- tanh(mean(atanh(apply(cc, 2, function(p) cor(W[[p[1]]], W[[p[2]]])))))
    res[[length(res) + 1]] <- data.frame(construct = c, effect = eff, target = t,
                                         self_other = so, other_other = oo)
  }
  ap <- do.call(rbind, res)
  rr <- list()
  for (c in CONSTRUCTS) {
    v <- long[long$variable == paste0(c, "rel"), ]
    rs <- c()
    for (dy in DYADS) {
      W <- reshape(v[v$target == dy, c("family", "source", "estimate")], idvar = "family",
                   timevar = "source", direction = "wide")
      cc <- combn(ROLES, 2)
      rs <- c(rs, apply(cc, 2, function(p) cor(W[[paste0("estimate.", p[1])]], W[[paste0("estimate.", p[2])]])))
    }
    rr[[c]] <- mean(rs)
  }
  list(actor_partner = ap, relationship = unlist(rr))
}
