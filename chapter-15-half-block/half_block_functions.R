# Half-block SRM functions (Chapter 15)
# A half-block has actors (rows) who rate partners (columns); partners do not rate actors.
# Each block is analyzed separately and the block estimates are pooled (Kenny's approach,
# as in Blocko). Base R only. lme4 is used only by hb_lme4().
#
# Data: long format, one row per directed rating, with columns
#   block, actor, partner, and one column per indicator (e.g., warm, kind).

# Actor x partner matrix for one block and one variable
hb_matrix <- function(x, var) {
  a <- sort(unique(x$actor)); p <- sort(unique(x$partner))
  X <- matrix(NA_real_, length(a), length(p), dimnames = list(a, p))
  X[cbind(match(x$actor, a), match(x$partner, p))] <- x[[var]]
  X
}

# Random-effects ANOVA estimates for one complete block (Chapter 6, Eqs. 6.1-6.3)
#   actor variance        = (MS_actor   - MS_AxP) / n_partners
#   partner variance      = (MS_partner - MS_AxP) / n_actors
#   relationship + error  =  MS_AxP
hb_one <- function(X) {
  if (anyNA(X)) stop("hb_one() needs a complete block; use hb_lme4() for missing data.")
  na <- nrow(X); np <- ncol(X)
  M <- mean(X); a <- rowMeans(X) - M; b <- colMeans(X) - M
  G <- X - outer(a, rep(1, np)) - outer(rep(1, na), b) - M
  MSA <- np * sum(a^2) / (na - 1)
  MSP <- na * sum(b^2) / (np - 1)
  MSE <- sum(G^2) / ((na - 1) * (np - 1))
  list(mean = M, a = a, b = b, G = G,
       var = c(actor = (MSA - MSE) / np, partner = (MSP - MSE) / na, relationship = MSE))
}

# Covariances between two variables in the same block, by polarization:
# Cov(x, y) = [V(x + y) - V(x - y)] / 4 (exact, because the estimators are quadratic forms)
hb_cov <- function(X, Y) (hb_one(X + Y)$var - hb_one(X - Y)$var) / 4

# Estimates for one block.
# constructs: named list; each element is a vector of 1 or 2 indicator names, e.g.
#   list(warmth = c("warm", "kind"), competence = c("capable", "effective"))
# With 2 indicators, the construct variance of each component is the covariance of the
# component across the two indicators, which removes error (Kenny, 2020):
#   relationship = Cov(g_ind1, g_ind2);  error = mean(MS_AxP) - relationship.
# With 1 indicator, relationship and error are confounded ("relationship_error").
hb_block <- function(x, constructs) {
  out <- c()
  for (cn in names(constructs)) {
    ind <- constructs[[cn]]
    for (v in ind) {                                   # indicator-level (univariate)
      e <- hb_one(hb_matrix(x, v))$var
      out[paste0(v, c(".actor", ".partner", ".relationship_error"))] <- e
    }
    if (length(ind) == 2) {
      c12 <- hb_cov(hb_matrix(x, ind[1]), hb_matrix(x, ind[2]))
      v1 <- hb_one(hb_matrix(x, ind[1]))$var; v2 <- hb_one(hb_matrix(x, ind[2]))$var
      out[paste0(cn, ".actor")]        <- c12[["actor"]]
      out[paste0(cn, ".partner")]      <- c12[["partner"]]
      out[paste0(cn, ".relationship")] <- c12[["relationship"]]
      out[paste0(cn, ".error")]        <- mean(c(v1[["relationship"]], v2[["relationship"]])) - c12[["relationship"]]
    }
  }
  # Covariances between constructs: average over the cross-construct indicator pairs
  cns <- names(constructs)
  if (length(cns) > 1) for (i in 1:(length(cns) - 1)) for (j in (i + 1):length(cns)) {
    cc <- sapply(constructs[[cns[i]]], function(u) sapply(constructs[[cns[j]]], function(w)
      hb_cov(hb_matrix(x, u), hb_matrix(x, w))), simplify = "array")
    cc <- matrix(cc, nrow = 3)   # rows: actor, partner, relationship; columns: indicator pairs
    nm <- paste0(cns[i], "_", cns[j])
    out[paste0(nm, ".actor_cov")]        <- mean(cc[1, ])
    out[paste0(nm, ".partner_cov")]      <- mean(cc[2, ])
    out[paste0(nm, ".relationship_cov")] <- mean(cc[3, ])
  }
  out
}

# All blocks: block-level estimates, pooled estimates, tests, standardized values,
# and correlations. Variances are tested one-tailed and covariances two-tailed
# with one-sample t tests across blocks (df = blocks - 1).
hb_all <- function(d, constructs) {
  chk <- split(d, d$block)
  B <- do.call(rbind, lapply(chk, hb_block, constructs = constructs))
  K <- nrow(B)
  est <- colMeans(B); se <- apply(B, 2, sd) / sqrt(K); t <- est / se
  is_cov <- grepl("_cov$", names(est))
  p <- ifelse(is_cov, 2 * pt(-abs(t), K - 1), pt(t, K - 1, lower.tail = FALSE))
  testable <- !grepl("relationship_error$|\\.error$", names(est))
  res <- data.frame(term = names(est), estimate = est, SE = se,
                    t = ifelse(testable, t, NA), df = K - 1, p = ifelse(testable, p, NA),
                    row.names = NULL)
  # Standardized (proportion of total) for each construct and indicator
  res$std <- NA_real_
  for (cn in c(names(constructs), unlist(constructs))) {
    rows <- grep(paste0("^", cn, "\\.(actor|partner|relationship|error|relationship_error)$"), res$term)
    res$std[rows] <- res$estimate[rows] / sum(res$estimate[rows])
  }
  # Correlations between constructs
  res$r <- NA_real_
  for (k in which(is_cov)) {
    nm <- sub("\\..*$", "", res$term[k]); comp <- sub("^.*\\.(\\w+)_cov$", "\\1", res$term[k])
    cn <- strsplit(nm, "_")[[1]]
    v1 <- est[paste0(cn[1], ".", comp)]; v2 <- est[paste0(cn[2], ".", comp)]
    res$r[k] <- if (v1 > 0 && v2 > 0) est[k] / sqrt(v1 * v2) else NA  # undefined if a variance is <= 0
  }
  list(blocks = data.frame(block = names(chk), B, row.names = NULL), summary = res)
}

# Actor and partner effects for second-stage modeling (one row per person and block)
hb_effects <- function(d, vars) {
  do.call(rbind, lapply(split(d, d$block), function(x) {
    A <- data.frame(block = x$block[1], id = sort(unique(x$actor)), role = "actor")
    P <- data.frame(block = x$block[1], id = sort(unique(x$partner)), role = "partner")
    for (v in vars) { h <- hb_one(hb_matrix(x, v)); A[[v]] <- h$a; P[[v]] <- h$b }
    rbind(A, P)
  }))
}

# Multilevel alternative (handles missing ratings and blocks of different sizes).
# With two indicators stacked, the actor:partner (dyad) term separates relationship
# from error, and the indicator-specific actor and partner terms keep that variance
# out of the error term.
hb_lme4 <- function(d, indicators) {
  if (!requireNamespace("lme4", quietly = TRUE)) stop("Install lme4.")
  long <- do.call(rbind, lapply(indicators, function(v)
    data.frame(block = d$block, actor = paste(d$block, d$actor), partner = paste(d$block, d$partner),
               indicator = v, y = d[[v]])))
  long$dyad <- paste(long$actor, long$partner)
  f <- if (length(indicators) == 2)
    y ~ indicator + (1 | actor) + (1 | partner) + (1 | dyad) + (1 | actor:indicator) + (1 | partner:indicator)
  else y ~ 1 + (1 | actor) + (1 | partner)
  m <- lme4::lmer(f, data = long, REML = TRUE)
  vc <- as.data.frame(lme4::VarCorr(m))[, c("grp", "vcov")]
  vc$grp <- c(actor = "actor", partner = "partner", dyad = "relationship",
              `actor:indicator` = "actor (indicator-specific)",
              `partner:indicator` = "partner (indicator-specific)",
              Residual = if (length(indicators) == 2) "error" else "relationship + error")[vc$grp]
  core <- vc$grp %in% c("actor", "partner", "relationship", "error", "relationship + error")
  vc$std <- ifelse(core, vc$vcov / sum(vc$vcov[core]), NA)
  list(model = m, components = vc)
}

# Own-twin (self-referenced) half-block, Table 15.1: actor i rates twin i' and the other
# twins in the block. 'own' is a 0/1 column marking the human's own twin.
# The own-twin effect for human i is gamma_ii' minus the mean of gamma_ij for the
# other twins. Because relationship effects sum to zero across a row, this equals
# gamma_ii' * n / (n - 1). Test the block means of these effects against zero across
# blocks (humans in the same block are not independent).
hb_own_twin <- function(d, var) {
  do.call(rbind, lapply(split(d, d$block), function(x) {
    X <- hb_matrix(x, var); G <- hb_one(X)$G; n <- ncol(G)
    own <- x[x$own == 1, c("actor", "partner")]
    r <- match(own$actor, rownames(G)); c <- match(own$partner, colnames(G))
    g_own <- G[cbind(r, c)]
    g_oth <- sapply(seq_along(r), function(k) mean(G[r[k], -c[k]]))
    x_own <- X[cbind(r, c)]
    x_oth <- sapply(seq_along(r), function(k) mean(X[r[k], -c[k]]))
    data.frame(block = x$block[1], actor = own$actor,
               rating_own = x_own, rating_others = x_oth,          # raw ratings
               gamma_own = g_own, gamma_others = g_oth,            # relationship effects
               own_effect = g_own - g_oth)                         # = gamma_own * n / (n - 1)
  }))
}

# Side-by-side comparison of own-twin and other-twin ratings and relationship effects,
# tested across blocks (block means). Note: rater-level own-twin effects also contain
# relationship and error variance, so their raw spread is not a test of rater differences.
hb_own_compare <- function(d, var) {
  o <- hb_own_twin(d, var)
  bm <- aggregate(o[, c("rating_own", "rating_others", "gamma_own", "gamma_others", "own_effect")],
                  list(block = o$block), mean)
  tt <- t.test(bm$own_effect, alternative = "greater")
  data.frame(variable = var,
             mean_rating_own = mean(o$rating_own), mean_rating_others = mean(o$rating_others),
             mean_gamma_own = mean(o$gamma_own), mean_gamma_others = mean(o$gamma_others),
             own_effect = mean(o$own_effect), t = unname(tt$statistic), df = unname(tt$parameter),
             p_one_tailed = tt$p.value)
}
