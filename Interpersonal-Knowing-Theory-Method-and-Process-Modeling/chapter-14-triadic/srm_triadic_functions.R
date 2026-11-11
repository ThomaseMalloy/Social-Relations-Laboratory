# ==============================================================================
# Chapter 14: Base R functions for the triadic idiographic key person design
#
# Sourced by the numbered scripts. Each round robin is one key person (KP) in
# one social group: the KP's triadic assessments of every focal person's
# (KP, N1, N2, N3, ...) view of every other focal person's relationship with
# that person's partner (Table 14.2). Rows are perceivers (as assessed by the
# KP), columns are targets.
#
# Estimators are the round robin equations of Chapter 4 (Warner, Kenny, &
# Stoto, 1979), the same estimates TripleR reports (Schönbrodt et al., 2012).
# ==============================================================================

# ---- Reading data -----------------------------------------------------------
# Long format: one row per KP x group x perceiver x target, with one column per
# construct. Wide format: one row per KP x group, with one column per cell
# named <construct><sep><perceiver><sep><target>, e.g. sat_KP_N1, sat_N1_KP.
wide_to_long <- function(w, kp_col, group_col, constructs, sep = "_") {
  out <- list()
  for (v in constructs) {
    cols <- grep(paste0("^", v, sep), names(w), value = TRUE)
    if (!length(cols)) stop("No columns start with '", v, sep, "'")
    parts <- strsplit(sub(paste0("^", v, sep), "", cols), sep, fixed = TRUE)
    if (any(lengths(parts) != 2)) stop("Column names must be ", v, sep, "<perceiver>", sep, "<target>")
    for (k in seq_along(cols)) {
      out[[length(out) + 1]] <- data.frame(
        kp.id = w[[kp_col]], group = w[[group_col]],
        perceiver = parts[[k]][1], target = parts[[k]][2],
        construct = v, value = w[[cols[k]]])
    }
  }
  d <- do.call(rbind, out)
  d <- reshape(d, idvar = c("kp.id", "group", "perceiver", "target"),
               timevar = "construct", direction = "wide")
  names(d) <- sub("^value\\.", "", names(d))
  rownames(d) <- NULL
  d
}

check_triadic_data <- function(d, constructs) {
  need <- c("kp.id", "group", "perceiver", "target", constructs)
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  if (any(d$perceiver == d$target)) stop("Rows with perceiver == target: the diagonal is not collected (Table 14.2).")
  key <- paste(d$kp.id, d$group, d$perceiver, d$target)
  if (anyDuplicated(key)) stop("Duplicate KP x group x perceiver x target rows.")
  rr <- split(d, list(d$kp.id, d$group), drop = TRUE)
  for (nm in names(rr)) {
    x <- rr[[nm]]; ids <- unique(c(x$perceiver, x$target)); n <- length(ids)
    if (n < 4) stop("Round robin ", nm, " has ", n, " focal people; at least 4 (KP + 3 nominees) are needed.")
    if (nrow(x) != n * (n - 1)) stop("Round robin ", nm, " is incomplete: ", nrow(x), " of ", n * (n - 1), " cells.")
    if (!"KP" %in% ids) stop("Round robin ", nm, " has no focal person labelled 'KP'.")
    if (anyNA(x[constructs])) stop("Round robin ", nm, " has missing values.")
  }
  invisible(TRUE)
}

# ---- One round robin --------------------------------------------------------
# X: n x n matrix, NA on the diagonal. Returns grand mean, effects, variances.
rr_matrix <- function(x, var) {
  ids <- unique(c("KP", sort(setdiff(unique(c(x$perceiver, x$target)), "KP"))))
  n <- length(ids)
  X <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
  X[cbind(match(x$perceiver, ids), match(x$target, ids))] <- x[[var]]
  X
}

srm_one <- function(X) {
  n  <- nrow(X); ids <- rownames(X)
  Mi <- rowMeans(X, na.rm = TRUE); Mj <- colMeans(X, na.rm = TRUE)
  M  <- mean(X, na.rm = TRUE)
  k1 <- (n - 1)^2 / (n * (n - 2)); k2 <- (n - 1) / (n * (n - 2)); k3 <- (n - 1) / (n - 2)
  a  <- k1 * Mi + k2 * Mj - k3 * M                                   # perceiver (4.1)
  b  <- k1 * Mj + k2 * Mi - k3 * M                                   # target    (4.2)
  G  <- X - outer(a, rep(1, n)) - outer(rep(1, n), b) - M            # relationship (4.3)
  # Variance components (Warner et al., 1979)
  up <- which(upper.tri(G), arr.ind = TRUE)
  gij <- G[up]; gji <- t(G)[up]
  dfB <- (n - 1) * (n - 2) / 2 - 1
  dfW <- (n - 1) * (n - 2) / 2
  MSB <- sum((gij + gji)^2 / 2) / dfB
  MSW <- sum((gij - gji)^2 / 2) / dfW
  var_g  <- (MSB + MSW) / 2
  cov_gg <- (MSB - MSW) / 2
  var_a  <- sum(a^2) / (n - 1) - MSB / (2 * (n - 2)) - MSW / (2 * n)
  var_b  <- sum(b^2) / (n - 1) - MSB / (2 * (n - 2)) - MSW / (2 * n)
  cov_ab <- sum(a * b) / (n - 1) - MSB / (2 * (n - 2)) + MSW / (2 * n)
  list(n = n, ids = ids, mean = M, a = setNames(a, ids), b = setNames(b, ids), G = G,
       var = c(perceiver = var_a, target = var_b, relationship = var_g,
               perceiver_target = cov_ab, dyadic = cov_gg))
}

# Bivariate SRM covariances from the univariate estimator (polarization):
# Cov(x, y) = [V(x + y) - V(x - y)] / 4. Because the estimators are quadratic
# forms in the data, this is exactly the bivariate round robin estimator.
srm_cov <- function(X, Y) {
  s <- srm_one(X + Y)$var; d <- srm_one(X - Y)$var
  c(perceiver_perceiver = (s[["perceiver"]] - d[["perceiver"]]) / 4,
    target_target       = (s[["target"]] - d[["target"]]) / 4,
    relationship_intra  = (s[["relationship"]] - d[["relationship"]]) / 4,  # Cov(g_ij x, g_ij y)
    relationship_inter  = (s[["dyadic"]] - d[["dyadic"]]) / 4)              # Cov(g_ij x, g_ji y)
}

# Sampling variance of a round robin grand mean given its SRM components:
# Var(M) = [(n-1)(s2a + s2b + 2 sab) + s2g + sgg] / (n(n-1))
grand_mean_sampling_var <- function(v, n)
  ((n - 1) * (v[["perceiver"]] + v[["target"]] + 2 * v[["perceiver_target"]]) +
     v[["relationship"]] + v[["dyadic"]]) / (n * (n - 1))

# ---- All round robins -------------------------------------------------------
# Returns the three second-stage data sets: person, relationship, kp.
srm_all <- function(d, constructs) {
  person <- list(); rel <- list(); kp <- list()
  rr <- split(d, list(d$kp.id, d$group), drop = TRUE)
  for (nm in names(rr)) {
    x <- rr[[nm]]; kp_id <- x$kp.id[1]; grp <- x$group[1]
    fits <- lapply(constructs, function(v) srm_one(rr_matrix(x, v))); names(fits) <- constructs
    ids <- fits[[1]]$ids; n <- length(ids)
    p <- data.frame(kp.id = kp_id, group = grp, focal = ids,
                    role = ifelse(ids == "KP", "KP", "nominee"))
    for (v in constructs) { p[[paste0(v, ".a")]] <- fits[[v]]$a; p[[paste0(v, ".b")]] <- fits[[v]]$b }
    person[[nm]] <- p
    off <- which(!diag(n) == 1, arr.ind = TRUE)
    r <- data.frame(kp.id = kp_id, group = grp, perceiver = ids[off[, 1]], target = ids[off[, 2]])
    for (v in constructs) r[[paste0(v, ".g")]] <- fits[[v]]$G[off]
    rel[[nm]] <- r
    k <- data.frame(kp.id = kp_id, group = grp, n = n)
    for (v in constructs) {
      f <- fits[[v]]; nom <- ids != "KP"
      k[[paste0(v, ".mean")]]       <- f$mean
      k[[paste0(v, ".a.kp")]]       <- f$a[["KP"]]
      k[[paste0(v, ".a.nom")]]      <- mean(f$a[nom])
      k[[paste0(v, ".b.kp")]]       <- f$b[["KP"]]
      k[[paste0(v, ".b.nom")]]      <- mean(f$b[nom])
      k[[paste0(v, ".var.a")]]      <- f$var[["perceiver"]]
      k[[paste0(v, ".var.b")]]      <- f$var[["target"]]
      k[[paste0(v, ".var.g")]]      <- f$var[["relationship"]]
      k[[paste0(v, ".var.total")]]  <- f$var[["perceiver"]] + f$var[["target"]] + f$var[["relationship"]]
      k[[paste0(v, ".cov.ab")]]     <- f$var[["perceiver_target"]]
      k[[paste0(v, ".cov.gg")]]     <- f$var[["dyadic"]]
      k[[paste0(v, ".var.mean")]]   <- grand_mean_sampling_var(f$var, n)
    }
    if (length(constructs) >= 2) {
      for (i in 1:(length(constructs) - 1)) for (j in (i + 1):length(constructs)) {
        cv <- srm_cov(rr_matrix(x, constructs[i]), rr_matrix(x, constructs[j]))
        tag <- paste0(constructs[i], ".", constructs[j])
        k[[paste0(tag, ".cov.aa")]]     <- cv[["perceiver_perceiver"]]
        k[[paste0(tag, ".cov.bb")]]     <- cv[["target_target"]]
        k[[paste0(tag, ".cov.gg.intra")]] <- cv[["relationship_intra"]]
        k[[paste0(tag, ".cov.gg.inter")]] <- cv[["relationship_inter"]]
      }
    }
    kp[[nm]] <- k
  }
  bind <- function(z) { z <- do.call(rbind, z); rownames(z) <- NULL; z[order(z$kp.id, z$group), ] }
  list(person = bind(person), relationship = bind(rel), kp = bind(kp))
}

# KP-level file in wide form: one row per KP, columns suffixed by group
kp_wide <- function(kp) {
  w <- reshape(kp, idvar = "kp.id", timevar = "group", direction = "wide", sep = "_")
  rownames(w) <- NULL
  w
}
