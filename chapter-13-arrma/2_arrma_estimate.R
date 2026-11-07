# ==============================================================================
# Chapter 13, Step 2: Estimate individual and dyadic ARRMA
#
# Reads SRM effects exported from TripleR (Step 1) and estimates the minimal
# ARRMA models two ways:
#   (a) the standard (scalar) equations: 13.10 and 13.11 (individual),
#       13.17 and 13.18 (dyadic)
#   (b) the matrix equations: theta-hat = R_zz^-1 r_zy, with
#       R_zz = Z'Z/(N - 1) and r_zy = Z'y/(N - 1) for standardized variables
# and checks that both agree with each other and with lm().
#
#   Individual ARRMA (Eq. 13.5):  alpha_mp = b* alpha_r + c* beta_r + e
#   Dyadic ARRMA     (Eq. 13.12): gamma_imp,ji = b* gamma_ij + c* gamma_ji + e
#     A = reciprocity, B = b* (unique assumed reciprocity),
#     C = c* (unique meta-accuracy)
#
# Input (data/):
#   arrma_person_effects.csv        one row per person
#   arrma_relationship_effects.csv  one row per directed dyad (actor i, partner j)
#     meta.r for actor i and partner j is i's prediction of j's response to i,
#     that is, gamma_imp,ji in the chapter's notation.
# Output (output/): arrma_results.csv and arrma_results.txt
#
# Run from the chapter-13-arrma folder. Base R only.
# ==============================================================================

if (!dir.exists("data")) stop("Set the working directory to the chapter-13-arrma folder first.")
dir.create("output", showWarnings = FALSE)

# ---- Settings: change these to match your TripleR export --------------------
person_file  <- "data/arrma_person_effects.csv"
rel_file     <- "data/arrma_relationship_effects.csv"
# person file columns
p_group <- "group.id"; p_id <- "id"
alpha_r_col  <- "judgment.a"   # perceiver (actor) effect in judgments
beta_r_col   <- "judgment.p"   # target (partner) effect in judgments
alpha_mp_col <- "meta.a"       # perceiver (actor) effect in metaperceptions
# relationship file columns
r_group <- "group.id"; r_actor <- "actor.id"; r_partner <- "partner.id"
gamma_r_col  <- "judgment.r"   # relationship effect in judgments
gamma_mp_col <- "meta.r"       # relationship effect in metaperceptions
n_boot <- 2000                 # bootstrap resamples of round robins (0 = skip)
set.seed(2026)
# -----------------------------------------------------------------------------

# ---- Read the effects -------------------------------------------------------
pe <- read.csv(person_file)
ind <- data.frame(group = pe[[p_group]], id = pe[[p_id]],
                  alpha_r = pe[[alpha_r_col]], beta_r = pe[[beta_r_col]],
                  alpha_mp = pe[[alpha_mp_col]])
ind <- ind[complete.cases(ind), ]

re <- read.csv(rel_file)
rel <- data.frame(group = re[[r_group]], i = as.character(re[[r_actor]]),
                  j = as.character(re[[r_partner]]),
                  g_r = re[[gamma_r_col]], g_mp = re[[gamma_mp_col]])
# Pair each directed dyad (i, j) with its reverse (j, i)
rev <- rel[, c("group", "i", "j", "g_r")]
names(rev) <- c("group", "j", "i", "g_ji")
dy <- merge(rel, rev, by = c("group", "i", "j"))
dy <- data.frame(group = dy$group, i = dy$i, j = dy$j,
                 g_ij = dy$g_r,        # i's unique response to j
                 g_ji = dy$g_ji,       # j's unique response to i
                 g_impji = dy$g_mp)    # i's prediction of j's unique response to i
dy <- dy[complete.cases(dy), ]

# ---- (a) Standard (scalar) equations ----------------------------------------
# x1 = the actor's own response, x2 = the partner's response, y = metaperception
arrma_scalar <- function(x1, x2, y) {
  rA <- cor(x1, x2)   # reciprocity                         (A)
  rB <- cor(x1, y)    # raw assumed reciprocity correlation
  rC <- cor(x2, y)    # raw meta-accuracy correlation
  b  <- (rB - rA * rC) / (1 - rA^2)   # Eq. 13.10 / 13.17
  c  <- (rC - rA * rB) / (1 - rA^2)   # Eq. 13.11 / 13.18
  c(reciprocity = rA, raw_AR = rB, raw_MA = rC, b_star = b, c_star = c,
    R2 = b * rB + c * rC)
}

# ---- (b) Matrix equations ---------------------------------------------------
arrma_matrix <- function(x1, x2, y) {
  N   <- length(y)
  Z   <- cbind(scale(x1), scale(x2))       # n x 2 standardized design matrix
  yz  <- scale(y)
  Rzz <- crossprod(Z) / (N - 1)            # predictor correlation matrix
  rzy <- crossprod(Z, yz) / (N - 1)        # predictor-criterion correlations
  rho <- Rzz[1, 2]
  Rinv_closed <- (1 / (1 - rho^2)) * matrix(c(1, -rho, -rho, 1), 2)  # closed-form 2x2 inverse
  stopifnot(max(abs(Rinv_closed - solve(Rzz))) < 1e-10)
  theta <- solve(Rzz, rzy)                 # theta-hat = R_zz^-1 r_zy
  list(theta = setNames(as.vector(theta), c("b_star", "c_star")),
       Rzz = Rzz, rzy = rzy, Rinv = Rinv_closed,
       R2 = as.numeric(crossprod(theta, rzy)))
}

# ---- Inference --------------------------------------------------------------
# Individual level: SRM effects sum to zero within each round robin, so one
# degree of freedom is lost per round robin: df = N - G - k (k = 2 predictors).
t_tests <- function(theta, Rinv, R2, df) {
  se <- sqrt((1 - R2) / df * diag(Rinv))
  t  <- theta / se
  data.frame(estimate = theta, se = se, t = t, df = df, p = 2 * pt(-abs(t), df))
}
# Both levels: bootstrap by resampling whole round robins (people and dyads
# within a round robin are not independent).
boot_ci <- function(dat, f, B) {
  if (B == 0) return(NULL)
  groups <- unique(dat$group)
  est <- replicate(B, {
    g <- sample(groups, replace = TRUE)
    d <- do.call(rbind, lapply(g, function(x) dat[dat$group == x, ]))
    f(d)
  })
  t(apply(est, 1, quantile, c(.025, .975)))
}

# ---- Individual ARRMA -------------------------------------------------------
s_ind <- arrma_scalar(ind$alpha_r, ind$beta_r, ind$alpha_mp)
m_ind <- arrma_matrix(ind$alpha_r, ind$beta_r, ind$alpha_mp)
l_ind <- coef(lm(scale(alpha_mp) ~ scale(alpha_r) + scale(beta_r), data = ind))[-1]
stopifnot(max(abs(s_ind[c("b_star", "c_star")] - m_ind$theta)) < 1e-10,
          max(abs(l_ind - m_ind$theta)) < 1e-10)
N <- nrow(ind); G <- length(unique(ind$group))
tt_ind <- t_tests(m_ind$theta, m_ind$Rinv, m_ind$R2, N - G - 2)
ci_ind <- boot_ci(ind, function(d) arrma_scalar(d$alpha_r, d$beta_r, d$alpha_mp), n_boot)

# ---- Dyadic ARRMA -----------------------------------------------------------
s_dy <- arrma_scalar(dy$g_ij, dy$g_ji, dy$g_impji)
m_dy <- arrma_matrix(dy$g_ij, dy$g_ji, dy$g_impji)
l_dy <- coef(lm(scale(g_impji) ~ scale(g_ij) + scale(g_ji), data = dy))[-1]
stopifnot(max(abs(s_dy[c("b_star", "c_star")] - m_dy$theta)) < 1e-10,
          max(abs(l_dy - m_dy$theta)) < 1e-10)
ci_dy <- boot_ci(dy, function(d) arrma_scalar(d$g_ij, d$g_ji, d$g_impji), n_boot)

# ---- Report -----------------------------------------------------------------
fmt <- function(x) sprintf("%.3f", x)
report <- function(label, s, m, ci, tt = NULL, nobs) {
  cat("\n================", label, "================\n")
  cat("Observations:", nobs, "\n\n")
  cat("Raw correlations\n")
  cat("  A  Reciprocity                     ", fmt(s["reciprocity"]), "\n")
  cat("     Raw assumed reciprocity (x1, y)  ", fmt(s["raw_AR"]), "\n")
  cat("     Raw meta-accuracy (x2, y)        ", fmt(s["raw_MA"]), "\n\n")
  cat("Standard equations (13.10/13.11 or 13.17/13.18)\n")
  cat("  B  b* unique assumed reciprocity    ", fmt(s["b_star"]), "\n")
  cat("  C  c* unique meta-accuracy          ", fmt(s["c_star"]), "\n\n")
  cat("Matrix equations: theta-hat = R_zz^-1 r_zy\n")
  cat("  R_zz =\n"); print(round(m$Rzz, 3))
  cat("  r_zy =", fmt(m$rzy), "\n")
  cat("  R_zz^-1 = 1/(1 - rho^2) [1, -rho; -rho, 1] =\n"); print(round(m$Rinv, 3))
  cat("  theta-hat = [b*, c*] =", fmt(m$theta), "  (identical to the standard equations)\n")
  cat("  R-squared =", fmt(m$R2), "\n")
  if (!is.null(tt)) {
    cat("\nTests (df = N - G - 2)\n"); print(round(tt, 4))
  }
  if (!is.null(ci)) {
    cat("\n95% bootstrap CIs, resampling round robins\n")
    print(round(ci[c("reciprocity", "b_star", "c_star", "raw_MA"), ], 3))
  }
}
sink("output/arrma_results.txt", split = TRUE)
cat("ARRMA estimates from", person_file, "and", rel_file, "\n")
report("Individual ARRMA (Eq. 13.5)", s_ind, m_ind, ci_ind, tt_ind, N)
cat("  x1 = alpha_r (perceiver effect, judgment), x2 = beta_r (target effect,\n",
    " judgment), y = alpha_mp (perceiver effect, metaperception)\n")
report("Dyadic ARRMA (Eq. 13.12)", s_dy, m_dy, ci_dy, NULL, nrow(dy))
cat("  x1 = gamma_ij, x2 = gamma_ji, y = gamma_imp,ji. Each dyad enters twice\n",
    " (as i,j and as j,i), so use the bootstrap CIs, not lm() standard errors.\n")
sink()

res <- rbind(data.frame(level = "individual", t(s_ind)), data.frame(level = "dyadic", t(s_dy)))
write.csv(res, "output/arrma_results.csv", row.names = FALSE)
cat("\nWrote output/arrma_results.txt and output/arrma_results.csv\n")
