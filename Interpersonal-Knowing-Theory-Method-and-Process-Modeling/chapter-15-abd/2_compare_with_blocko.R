# Check: this R code against Blocko output for the same data (ABD50N6STATUS.out).
# Blocko values are copied from that output (absolute variances, relative stable
# construct variances, percent stable variance, and reciprocity correlations).
source("abd_functions.R")
d <- read.csv("data/claude_status_abd.csv")
res <- abd_all(d, list(warmth = c("warm", "kind"), competence = c("capable", "effective")), c("High", "Low"))
S <- res$summary
v <- function(tm, col = "estimate") S[[col]][S$term == tm]
blocko <- rbind(
  # direction, variable, actor, partner, relationship (absolute variances)
  c("High_Low", "warm", .176, .118, .424), c("High_Low", "kind", .219, .126, .403),
  c("High_Low", "capable", .161, .128, .352), c("High_Low", "effective", .140, .111, .358),
  c("Low_High", "warm", .204, .147, .527), c("Low_High", "kind", .214, .112, .499),
  c("Low_High", "capable", .110, .179, .368), c("Low_High", "effective", .110, .161, .339),
  c("High_Low", "warmth", .196, .121, .350), c("High_Low", "competence", .151, .109, .269),
  c("Low_High", "warmth", .203, .122, .439), c("Low_High", "competence", .113, .169, .263))
cmp <- do.call(rbind, lapply(1:nrow(blocko), function(i) {
  b <- blocko[i, ]; comp <- c("actor", "partner", "relationship")
  r <- sapply(comp, function(cc) v(paste(b[1], b[2], cc, sep = ".")))
  data.frame(direction = b[1], variable = b[2], component = comp,
             blocko = as.numeric(b[3:5]), R = round(r, 3))
}))
rel <- data.frame(direction = rep(c("High_Low", "Low_High"), each = 6),
                  variable = rep(rep(c("warmth", "competence"), each = 3), 2),
                  component = rep(c("actor", "partner", "relationship"), 4),
                  blocko = c(.267, .165, .478, .241, .175, .431, .238, .143, .516, .178, .265, .414))
rel$R <- round(mapply(function(a, b, c) v(paste(a, b, c, sep = "."), "relative"), rel$direction, rel$variable, rel$component), 3)
rec <- data.frame(term = paste("recip", rep(c("warm", "kind", "capable", "effective", "warmth", "competence"), each = 3),
                               c("High", "Low", "relationship"), sep = "."),
                  blocko = c(-.014, .007, -.073, .106, .176, -.102, -.046, .000, -.028, .037, .171, -.019,
                             .056, .122, -.100, .010, .070, -.019))
rec$R <- round(sapply(rec$term, function(tm) v(tm, "r")), 3)
stab <- data.frame(res$stable[, 1:2], blocko = c(.910, .847, .897, .857), R = round(res$stable$stable, 3))
cat("Absolute variances\n"); print(cmp, row.names = FALSE)
cat("\nRelative stable construct variances\n"); print(rel, row.names = FALSE)
cat("\nReciprocity correlations\n"); print(rec, row.names = FALSE)
cat("\nPercent stable variance\n"); print(stab, row.names = FALSE)
all_diff <- c(cmp$R - cmp$blocko, rel$R - rel$blocko, rec$R - rec$blocko, stab$R - stab$blocko)
cat(sprintf("\n%d values compared; largest absolute difference = %.3f\n", length(all_diff), max(abs(all_diff))))
