# Asymmetric block design (ABD) SRM functions (Chapter 15)
# An ABD block has two groups, A and B. Members of A rate members of B (A->B) and members
# of B rate members of A (B->A); no one rates members of their own group.
# Each direction is a half-block. Each block is analyzed separately and the block
# estimates are pooled across blocks, as in Blocko (Kenny & Xuan). Base R only.
#
# Data: long format, one row per directed rating, with columns
#   block, actor, actor_group, partner, and one column per indicator.

# ---- half-block building blocks --------------------------------------------------------
abd_matrix <- function(x, var, rows = sort(unique(x$actor)), cols = sort(unique(x$partner))) {
  X <- matrix(NA_real_, length(rows), length(cols), dimnames = list(rows, cols))
  X[cbind(match(x$actor, rows), match(x$partner, cols))] <- x[[var]]
  X
}
# Random-effects ANOVA for one complete direction (Chapter 6, Eqs. 6.1-6.3)
hb_one <- function(X) {
  if (anyNA(X)) stop("A block has missing ratings. Complete blocks are required (as in Blocko).")
  na <- nrow(X); np <- ncol(X)
  M <- mean(X); a <- rowMeans(X) - M; b <- colMeans(X) - M
  G <- X - outer(a, rep(1, np)) - outer(rep(1, na), b) - M
  MSA <- np * sum(a^2) / (na - 1); MSP <- na * sum(b^2) / (np - 1)
  MSE <- sum(G^2) / ((na - 1) * (np - 1))
  list(a = a, b = b, G = G,
       var = c(actor = (MSA - MSE) / np, partner = (MSP - MSE) / na, relationship = MSE))
}
# Covariance of two variables in the same direction, by polarization
hb_cov <- function(X, Y) (hb_one(X + Y)$var - hb_one(X - Y)$var) / 4

# Reciprocity between the two directions for one block.
# X = ratings A->B (rows A, cols B); Y = ratings B->A (rows B, cols A).
#   relationship (dyadic) covariance: Cov(g_ij, g_ji) = sum(gij * gji) / ((nA-1)(nB-1))
#   A members' actor-partner covariance: sum(a_i * b_i) / (nA - 1) - Cov(g, g') / nB
#   B members' actor-partner covariance: sum(a_j * b_j) / (nB - 1) - Cov(g, g') / nA
abd_recip <- function(X, Y) {
  hx <- hb_one(X); hy <- hb_one(Y); nA <- nrow(X); nB <- ncol(X)
  Gy <- t(hy$G)[rownames(X), colnames(X)]
  cgg <- sum(hx$G * Gy) / ((nA - 1) * (nB - 1))
  c(A = sum(hx$a * hy$b[rownames(X)]) / (nA - 1) - cgg / nB,
    B = sum(hy$a[colnames(X)] * hx$b) / (nB - 1) - cgg / nA,
    relationship = cgg)
}

# ---- one block ----------------------------------------------------------------------------
# constructs: named list, each 1 or 2 indicator names.
# groups: c(A, B), the two values of actor_group (e.g., c("High", "Low")).
abd_block <- function(x, constructs, groups) {
  A <- groups[1]; B <- groups[2]
  xa <- x[x$actor_group == A, ]; xb <- x[x$actor_group == B, ]
  rA <- sort(unique(xa$actor)); rB <- sort(unique(xb$actor))
  M <- function(d, v, r, c) abd_matrix(d, v, r, c)
  out <- c()
  dirs <- list(list(nm = paste0(A, "_", B), d = xa, r = rA, c = rB),
               list(nm = paste0(B, "_", A), d = xb, r = rB, c = rA))
  for (cn in names(constructs)) {
    ind <- constructs[[cn]]
    for (dr in dirs) {
      v <- lapply(ind, function(i) hb_one(M(dr$d, i, dr$r, dr$c))$var)
      for (k in seq_along(ind))
        out[paste(dr$nm, ind[k], c("actor", "partner", "relationship"), sep = ".")] <- v[[k]]
      if (length(ind) == 2) {
        s <- hb_cov(M(dr$d, ind[1], dr$r, dr$c), M(dr$d, ind[2], dr$r, dr$c))   # stable
        u <- (v[[1]] + v[[2]]) / 2 - s                                           # unstable
        out[paste(dr$nm, cn, c("actor", "partner", "relationship"), sep = ".")] <- s
        out[paste(dr$nm, cn, c("actor_unstable", "partner_unstable", "relationship_unstable"), sep = ".")] <- u
      }
    }
    # reciprocity: each indicator, and the construct (average of the two cross-indicator pairs)
    XA <- function(i) M(xa, i, rA, rB); YB <- function(i) M(xb, i, rB, rA)
    for (i in ind) out[paste("recip", i, c(A, B, "relationship"), sep = ".")] <- abd_recip(XA(i), YB(i))
    if (length(ind) == 2)
      out[paste("recip", cn, c(A, B, "relationship"), sep = ".")] <-
        (abd_recip(XA(ind[1]), YB(ind[2])) + abd_recip(XA(ind[2]), YB(ind[1]))) / 2
  }
  out
}

# ---- all blocks ---------------------------------------------------------------------------
# Returns block-level estimates and a summary table: pooled estimates, SE, t, p
# (variances one-tailed, covariances two-tailed, df = blocks - 1), relative variances,
# percent stable variance, and reciprocity correlations.
abd_all <- function(d, constructs, groups) {
  chk <- split(d, d$block)
  Bk <- do.call(rbind, lapply(chk, abd_block, constructs = constructs, groups = groups))
  K <- nrow(Bk); est <- colMeans(Bk); se <- apply(Bk, 2, sd) / sqrt(K); t <- est / se
  is_rec <- grepl("^recip\\.", names(est)); is_uns <- grepl("_unstable$", names(est))
  p <- ifelse(is_rec, 2 * pt(-abs(t), K - 1), pt(t, K - 1, lower.tail = FALSE))
  S <- data.frame(term = names(est), estimate = est, SE = se, t = t, df = K - 1, p = p, row.names = NULL)
  S$t[is_uns] <- NA; S$p[is_uns] <- NA
  # Indicator relationship includes error, so it is not tested
  ind_rel <- grepl("\\.relationship$", S$term) & !is_rec &
    sapply(strsplit(S$term, "\\."), `[`, 2) %in% unlist(constructs)
  S$t[ind_rel] <- NA; S$p[ind_rel] <- NA
  # Relative variance (Blocko's definition): share of the total, where the total for a
  # construct is stable + unstable variance (the mean indicator total). As in Blocko,
  # negative variance estimates are set to zero for this calculation.
  S$relative <- NA_real_
  dirs <- c(paste0(groups[1], "_", groups[2]), paste0(groups[2], "_", groups[1]))
  for (dr in dirs) for (v in c(names(constructs), unlist(constructs))) {
    rows <- grep(paste0("^", dr, "\\.", v, "\\.(actor|partner|relationship)(_unstable)?$"), S$term)
    if (length(rows)) { e <- pmax(S$estimate[rows], 0); S$relative[rows] <- e / sum(e) }
  }
  # Reciprocity correlations
  S$r <- NA_real_
  g <- function(tm) { x <- S$estimate[S$term == tm]; if (length(x)) x else NA }
  for (k in which(is_rec)) {
    pr <- strsplit(S$term[k], "\\.")[[1]]; v <- pr[2]; who <- pr[3]
    if (who == "relationship") {
      v1 <- g(paste(dirs[1], v, "relationship", sep = ".")); v2 <- g(paste(dirs[2], v, "relationship", sep = "."))
    } else {
      own <- if (who == groups[1]) dirs[1] else dirs[2]; oth <- setdiff(dirs, own)
      v1 <- g(paste(own, v, "actor", sep = ".")); v2 <- g(paste(oth, v, "partner", sep = "."))
    }
    S$r[k] <- if (isTRUE(v1 > 0 && v2 > 0)) S$estimate[k] / sqrt(v1 * v2) else NA
  }
  # Percent stable construct variance
  stable <- do.call(rbind, lapply(dirs, function(dr) do.call(rbind, lapply(names(constructs)[lengths(constructs) == 2], function(cn) {
    s <- sum(pmax(est[paste(dr, cn, c("actor", "partner", "relationship"), sep = ".")], 0))
    u <- sum(pmax(est[paste(dr, cn, c("actor_unstable", "partner_unstable", "relationship_unstable"), sep = ".")], 0))
    data.frame(direction = dr, construct = cn, stable = s / (s + u))
  }))))
  list(blocks = data.frame(block = names(chk), Bk, row.names = NULL), summary = S, stable = stable)
}

# Actor and partner effects for second-stage modeling (one row per person, block, direction)
abd_effects <- function(d, vars, groups) {
  do.call(rbind, lapply(split(d, d$block), function(x) do.call(rbind, lapply(1:2, function(k) {
    g <- groups[k]; xd <- x[x$actor_group == g, ]
    rows <- sort(unique(xd$actor)); cols <- sort(unique(xd$partner))
    A <- data.frame(block = x$block[1], id = rows, group = g, role = "actor", direction = paste0(g, "_", groups[3 - k]))
    P <- data.frame(block = x$block[1], id = cols, group = groups[3 - k], role = "partner", direction = A$direction[1])
    for (v in vars) { h <- hb_one(abd_matrix(xd, v, rows, cols)); A[[v]] <- h$a; P[[v]] <- h$b }
    rbind(A, P)
  }))))
}

# Paired comparison of block estimates between two conditions (e.g., 0 vs. 2 turns)
abd_compare <- function(blocks1, blocks2, term) {
  m <- merge(blocks1[, c("block", term)], blocks2[, c("block", term)], by = "block")
  tt <- t.test(m[[3]], m[[2]], paired = TRUE)
  data.frame(term = term, mean_1 = mean(m[[2]]), mean_2 = mean(m[[3]]), diff = mean(m[[3]] - m[[2]]),
             t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value)
}
